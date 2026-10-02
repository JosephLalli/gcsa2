# External-memory GCSA2 construction

Status: the production disk-first route is implemented through prefix doubling,
bounded pruning, final component packing, streaming LCP, framed temporary-file
compression, and durable predecessor retirement. Whole-process input decoding
and global disk admission remain open work. The public `.gcsa` and `.lcp`
formats remain unchanged. This document records both the measured baseline and
the invariants of the disk-first route; the final section tracks which parts
are actually implemented.

## Baseline audit

The legacy builder uses temporary files, but a temporary file is generally a
serialization point between whole-file in-memory operations rather than a
bounded external-memory data structure.

| Location | Resident structures and behavior |
| --- | --- |
| `InputGraph::read()` (`src/files.cpp`) | Reserves and reads every `KMer` in the requested logical input. |
| `InputGraph::readKeys()` | Accumulates one global `vector<key_type>`, globally sorts it, then merges equal labels. Each physical input is also read into a whole-file `vector<KMer>`. |
| `InputGraph::readFrom()` | Accumulates one global `vector<node_type>`, then sorts and deduplicates it. Each input is again read as a whole-file `vector<KMer>`. |
| `PathGraph::PathGraph(const InputGraph&, ...)` | Reads one complete logical input, sorts all `KMer` records by label, converts label values to ranks, and writes one path/rank pair. |
| `PathGraph::extend()` | Reads one complete path file and rank sidecar, sorts all paths by `from`, and builds an in-memory `ValueIndex`. It reserves 1,048,576 `PathNode` records and `(2^step + 1) * 1,048,576` ranks **per OpenMP thread** before generating joins. |
| `PathGraphBuilder::write()` | The memory-limit branch checks whether the eventual output file could later be reloaded. It does not constrain the source vectors, `ValueIndex`, sort workspace, or concurrent thread-local buffers. |
| `PathGraphBuilder::sort()` | Closes a generated path/rank pair, reloads both entire files, sorts all paths by label, and rewrites them. |
| `PathGraph::prune()` | Merges label-sorted inputs, but `BufferWindow<PriorityNode>` retains a complete equal-label/extended range. A repetitive range can therefore be unbounded. |
| `PathGraphMerger` | Opens a `ReadBuffer<PathNode>` and `ReadBuffer<rank_type>` per input. Each buffer is sized as one million **elements**, not bytes. |
| `MergedGraph::MergedGraph()` | Streams most path data, but `SameFromSet::nodes` can grow with a merged group and output remains an ephemeral PID-named file set. |
| `GCSA::GCSA()` | Simultaneously retains one full BWT bitvector per alphabet character, outdegree/occurrence/redundancy counters, sampled-position and sample-boundary bitvectors, `sample_buffer`, and dense `prev_occ`, plus reader buffers. |
| `LCPArray::LCPArray()` | Allocates the leaf sequence and every range-minimum level in one `int_vector<0>`, initializes it at width 8, and only then bit-compresses it. |

The up-front estimate in `GCSA::GCSA()` and the future-file check in
`PathGraphBuilder::write()` are guards. They can terminate before an allocation,
but they do not provide an external execution route.

All ten suspected peaks in the implementation request were present in the
checked-out GCSA2 revision. Two details differed slightly from the initial
description: `PathGraph::extend()` allocated the large path/rank staging
vectors per OpenMP thread rather than as one shared pair, and
`MergedGraph::MergedGraph()` was already mostly streaming but retained an
unbounded `SameFromSet::nodes` group. The table describes the preserved legacy
route. The workspace-selected route now replaces key/start/initial-path
preprocessing, prefix-doubling joins and label sorting, final event/component
construction, dense previous-occurrence state, and simultaneous raw LCP levels.
It now also spills pruning equal-label/range state, externally reduces
`MergedGraph` same-from sets, bounds input/output shard descriptors, and derives
merged-graph output buffers from bytes rather than element counts. Remaining
resident limitations are called out in the implementation ledger below.

### Production fixture and measured legacy baseline

The chromosome-scale acceptance fixture is the packed chr20 graph pruned with
`vg prune -k 32` from the panSC indexing work. Its identities are fixed so a
benchmark cannot silently substitute a smaller or differently pruned graph:

| Input | Bytes | SHA-256 |
| --- | ---: | --- |
| `chr20.exact_walk_k32_m0.pruned.vg` | 440,223,843 | `04e566fb4810e273fcaad1c4796b4dc785a336796b604089bbae1df6007bc525` |
| `chr20.exact_walk_k32_m0.mapping` | 64,428,032 | `04843ec74911550700b051fed5579e3a7c29748309a7fb5b55c99dba17cd4a05` |

The immutable files are under
`/mnt/ssd/lalli/hprc_v2_vg_rna/smoke_chr20/exact_walk_k32_loss_diagnostic_v1_20260821T052743-0500/index/`.
No benchmark step regenerates or reprunes them.

The preferred second acceptance fixture is the existing chr21 k32-pruned retry
under
`/mnt/ssd/lalli/hprc_v2_vg_rna/smoke_chr21/exact_walk_k_sweep_v1_20260815T115236-0500/k32_gcsa_retry_cap400g/index/`:

| Input | Bytes | SHA-256 |
| --- | ---: | --- |
| `chr21.exact_walk_k32_m0.pruned.vg` | 246,263,907 | `1c75a17d98ceae0f1b0d7cf67649f91e1e4136f1929be8bf54b884759519129e` |
| `chr21.exact_walk_k32_m0.mapping` | 25,928,080 | `85e0763a8b1040569694272492259f69cc87168b918d04a0576811c030ca2a99` |

Its successful legacy `vg index -k16 -X4 -Z2048 -t32` run took 2:07:13
and 105,438,580 KiB maximum RSS (100.55 GiB).

The external-memory acceptance build completed on the same immutable fixture
with a 23 GiB GCSA2 working-set goal, a hard 25 GiB cgroup limit, swap disabled,
four process workers, and a 1 TiB disk limit. The initial process was stopped
only after step 4 had committed; `--gcsa-resume` restored that frontier and the
successor process completed normally. Active construction time across the two
processes was approximately 15:45 (16:15 elapsed including the planned pause),
favoring feasibility over the 2:07 legacy runtime. The resume process used
22,375,576 KiB maximum RSS (21.34 GiB), 78.8% below the legacy peak, while the
cgroup reached its exact 25 GiB ceiling through reclaimable file cache without
an OOM, OOM kill, or swap. Peak live run-directory bytes were 432,556,537,957
(402.85 GiB), and the completed retained run occupied about 262 GiB.

Step 4 generated 1,109,279,137 joined paths and streamed 1,208,172,107 records
directly into label runs. Grouped left-context records avoided 10.8 GiB of
transient payload, while direct join-to-label streaming avoided 97.3 GiB of
intermediate path/rank I/O in each direction. The final graph contained
399,778,113 paths, 413,060,024 edges, 631,824,030 pointers, and 177,682,537
samples. `vg index -V` queried 54,709,753 patterns and reported `Index
verification complete`.

| External output | Bytes | SHA-256 |
| --- | ---: | --- |
| GCSA | 1,357,870,701 | `2b1021eb926596de29c2cdd280e35b7f0c69d4d8a47ab9152e72fab81614b05f` |
| LCP | 406,123,913 | `7a74ef3b67cb0243cb9433ee240ee4fb7cf225b8f7775aea1a013b0e4eb11e7f` |

The acceptance binary used GCSA2 `30b4886`. Later component-streaming,
phase-lifetime, aggregate-budget sizing, and incremental-checksum commits pass
the fixture tests but are not included in the chromosome-scale measurements
above.

Only preexisting pruned fixtures are eligible for acceptance runs. The chr19
production prune exited successfully after 49:45:06 at 154,737,660 KiB maximum
RSS, but its wrapper did not publish the 748,778,080-byte `.partial` graph as a
committed fixture. No competing chr19 index build was launched. The remaining
whole-physical-graph input load, allocations outside the shared byte budget,
and global disk-admission gaps below must be closed or deployment-bounded before
that race is considered ready.

The successful legacy command used `vg index -k 16 -X 4 -Z 700 -t 32 -V`.
It produced a 682,869,673-byte GCSA and a 273,490,577-byte LCP in 1:55:35
wall time. `/usr/bin/time -v` measured 63,554,776 KiB maximum RSS (60.61 GiB),
while GCSA2 reported 223.140 GB read and 193.262 GB written. The often-quoted
“about 100 GB” was a conservative recollection rather than the recorded peak
for this exact fixture.

The external completion run resumed after a committed step-2 prune with a
96 GiB GCSA2 budget and a hard 128 GiB cgroup ceiling. It completed and passed
`vg index` verification in 2:51:15. Maximum process RSS was 33,936,664 KiB
(32.36 GiB), 46.6% below the legacy peak. Sampled cgroup peak was
82,382,069,760 bytes (76.73 GiB); live inspection showed that the difference
was inactive clean filesystem cache accumulated by the still-legacy final
merge, not anonymous heap. Swap was disabled and no memory-limit or OOM event
occurred. Peak live run-directory usage was 337,799,719,952 bytes (314.60 GiB),
below the explicit 1 TiB disk limit, and the retained completed workspace plus
outputs occupied about 177 GiB.

The completed GCSA and LCP are byte-identical to the legacy outputs:

| Output | Bytes | SHA-256 |
| --- | ---: | --- |
| GCSA | 682,869,673 | `f4d0f89a6f188802aeb5ad19c3cc4cd6d6ac9cbee892cef5776b8c3eaee45eca` |
| LCP | 273,490,577 | `2a98d96f0a19c1a169f408b7b93f7cfadcf0f2176b86da6c82c5cb6cf0476940` |

This was deliberately a feasibility-first completion run, not a clean-wall-
time comparison: it reused the committed frontier created by earlier forced-
spill attempts. The command itself read 451,411,742,720 bytes and wrote
541,563,432,960 bytes according to `/usr/bin/time`; GCSA2's logical counters
reported 497.497 GB read and 417.810 GB written. A separate 23 GiB internal /
25 GiB cgroup run forced multiple immutable join spills, reached only
10,824,118,272 bytes (10.08 GiB) sampled cgroup peak, and was deliberately
terminated before its incomplete join phase committed. Resumption retained its
predecessor and ignored the incomplete successor exactly as designed.

**Do not read the peak run-directory figures above as a workspace-to-input
multiplier.** They are the right numbers for what they measure — the high-water
mark a given run needed on disk — but three things make the ratio
`peak / pruned-input-bytes` non-transferable, and a downstream consumer has
already been misled by it. First, both peaks were reached by resumed,
deliberately forced-spill runs that carry a predecessor's committed frontier
alongside the successor's live state; the completed retained workspaces are
much smaller (about 262 GiB against the 402.85 GiB peak, and about 177 GiB
against 314.60 GiB). Second, spill volume is driven by the configured memory
budget, so the tighter the budget the larger the workspace: the 23 GiB run
peaked higher in absolute terms than the 96 GiB one despite a smaller input.
Third, the ratio is not even monotone in input size. Dividing peak by input
gives 767x for the chr20 fixture (337,799,719,952 / 440,223,843) and 1,756x for
chr21 (432,556,537,957 / 246,263,907), while a 96 GiB-budget build recording no
resume or forced spill, over a 4,231,235,421-byte pruned corpus, measured
60,717,348,762 bytes —
**14.35x** — settling to 4.66x after cleanup (`hprc_v2_vg_rna`,
`notes/chr21_or_vs_exact_dedup_downstream.md`, 2026-09-14). Size a disk budget from
a measurement at the intended budget, input shape and resume posture, not by
scaling any figure on this page.

The chromosome completion binary predates the later parallel run-sort,
direct-join, external-preprocessing, final-event, streaming-LCP, and dirty-cache
changes described below. Those changes pass unit and small integration tests,
but their chromosome-scale runtime and RSS have not been remeasured. In
particular, the newer 75/25 run allocation intentionally uses more of an
explicitly generous budget to eliminate an otherwise unnecessary full merge
pass; users seeking the smallest RSS should configure a smaller sort-run budget
and retain the hard cgroup ceiling.

### `vg` integration audit

`vg index` previously constructed `ConstructionParameters` with only doubling
steps and the legacy `-Z` size limit. Generated de Bruijn files were anonymous
temporary files and were deleted after construction. `vg autoindex` similarly
exposed a GCSA temporary-size setting but had no durable construction frontier.
The convenience helper in `src/build_index.cpp` still constructs an in-memory
index for small programmatic/test callers; it is not used by the production
`vg index` or autoindex recipes.

The external route now gives both production entry points a workspace, resume
flag, byte-valued memory limit, and byte-valued disk limit. `vg` records source
graph size and checksum, atomically commits its generated de Bruijn inputs into
`WORK/inputs`, and restores them only when their semantic manifest matches.
The GCSA and LCP outputs are each written to a partial path and renamed only
after successful serialization. Publishing the two-file pair through one
atomic marker remains future work.

### Logical inputs are not spill shards

`PriorityNode::file` currently denotes the original `PathGraph` file index.
`SameFromFile` requires both the same `from` node and the same file before it
allows a range to be merged. Equal paths from different source files take the
explicitly redundant branch in `PathGraph::prune()`. Consequently, splitting
one input chromosome into several apparent GCSA2 inputs changes pruning
semantics and index size.

The disk-first route uses distinct strong types:

```cpp
struct logical_file_id_t { std::uint32_t value; };
struct physical_shard_id_t { std::uint64_t value; };
```

All semantic comparisons use `logical_file_id_t`. A logical file may have any
number of physical shards. A shard identifier is only an artifact locator and
must not be implicitly convertible to a logical identifier.

## Resource invariants

One process-wide `MemoryBudget` grants byte reservations through move-only RAII
tokens. Buffers, decoded records, priority queues, hash tables, thread-local
output, compression workspaces, and bounded queues must reserve before
allocation. A task whose reservation does not fit waits; a single request larger
than the usable ceiling is an error in the algorithm or configuration and must
be subdivided. The configured ceiling includes an explicit safety margin for
allocator, OpenMP, SDSL, and library overhead.

I/O buffer settings are bytes. The element capacity is computed as
`max(1, bytes / sizeof(T))`. Compute concurrency and I/O concurrency are
scheduled separately, but both consume tokens from the same budget.

Linux cgroups charge clean filesystem page cache to `MemoryMax`. The external
preprocessor, join, label sorter, final-event streams, spillable per-path
start-node sets, disk arrays, and durable checkpoint writer/validator/restorer
therefore use bounded byte buffers and,
where implemented, sequential access advice plus rolling cache eviction.
Sequential join and final-event writers sync and evict completed prefixes every
512 MiB while retaining a 64 MiB tail; readers evict consumed prefixes. A
backwards seek for a pathological join key resets the reader watermark so every
replay remains bounded. External final-scan `ReadBuffer` instances use their
existing descriptor for `pread()` and retire consumed prefixes every 32 MiB
while retaining an 8 MiB tail. This adds no hidden descriptor per stream, and a
backward seek simply rereads an evicted page. `posix_fadvise()` is best-effort
and non-semantic;
`fdatasync()` is the durability boundary. Pruning and `MergedGraph` now use
bounded spill structures and byte-sized buffers. Every active framed-stream
encoder reserves its input block, worst-case output, zstd context, and, when the
block is large enough to occupy more than one zstd job, the workspace a
multi-threaded context really allocates: one compression context per engaged
worker, a round buffer of `workers + 3` jobs, and a job-output pool of
`2 * workers + 3` buffers, all sized from the job rather than the block. A block
shorter than one job engages one worker whatever was requested, so the writer
configures and the estimate charges only the usable workers; charging the
request instead reserved a job pool the encoder never allocated at large blocks
while under-charging the pool it did allocate at small ones. Framed readers
reserve encoded and decoded blocks plus decoder context before construction. Final
mapping/input allocations and some older pruning/merged-graph reservations are
not all admitted through the process-wide `MemoryBudget`. The in-process token
budget is therefore not yet an end-to-end RSS ceiling; a hard cgroup remains
the acceptance boundary.

`--max-open-files` is now enforced for the prune, final merge, and final-scan
run sets.
Pruning reserves its two possible spill descriptors before dividing the
remainder between two-descriptor input and output cache entries. Final merge
reserves fourteen descriptors for four sequential outputs, two possible spill
files, and the worst overlapping two-way `SameFromSet` external-sort state;
each remaining input-cache entry consumes two. The final scan additionally
accounts for all merged-graph readers, event writers, two mutable disk arrays,
one retained start-set reader, and the external sorter's `2 * fan-in + 4`
descriptors. Its fan-in is reduced to fit the ceiling. The supported global
minimum is 64; the isolated prune/merge algorithms can still operate with 16.
This ceiling covers GCSA2 construction files, not descriptors already inherited
from the embedding process.

### Compute and I/O concurrency

Large initial join and label runs are sorted in place with the existing OpenMP
parallel quicksort. The vector is already covered by the phase reservation, so
using the construction thread team does not add one run-sized allocation per
worker. Small runs remain sequential to avoid team-startup overhead. Run
generation uses three quarters of the available phase budget; the remaining
quarter covers allocator/OpenMP slack and bounded decoding and output-cache
state. A sole sorted run flows directly to final path/rank encoding instead of
being copied through a degenerate one-way merge. The 75/25 split is
intentionally chosen so the chr20 step-4 label set (about 67.1 GiB of 96-byte
records) is one run with a 96 GiB phase budget rather than narrowly exceeding a
64 GiB run and paying for a full merge pass.

Merge emitters remain single-writer, while prefix-doubling joins can now use
deterministic range partitions of `(logical_file_id, join_key)`. Worker tasks
reserve byte tokens, consume disjoint immutable input ranges, and publish
separate physical output shards. Those shards retain their shared logical ID
and are compacted by a byte-budgeted, descriptor-bounded k-way merge before
the next prune. Because every input already has the complete label order, the
compactor never decodes the frontier into a new general sort. It performs as
many merge levels as a small fan-in requires, publishes deterministic physical
shard ids, and may retain several sorted streams when one monolith is neither
needed nor desirable;
pruning itself projects its globally label-sorted merge directly to one output
shard per logical input. A heavy key is recursively split on both left and
right record ranges; splitting the left range suppresses duplicate sorted
bypasses, while splitting the right range gives each task a disjoint bypass.
Separate admission for compute and writers remains future work for filesystems
that cannot sustain all admitted writers concurrently.

The intended bounded pipeline is therefore read-ahead for partition N+1,
parallel join/sort for partition N, and ordered writeback for partition N-1.
Queues, decompression state, and prefetched pages count against the same global
budget. Independent merge groups can use the same scheduler. Work stealing may
change task completion order, but semantic range names and ordered publication
keep committed artifacts byte-deterministic.

Range tasks use a versioned multiprocessing ABI. The coordinator uses
`posix_spawn()` to `exec` workers with immutable input paths, exact key ranges, unique
partial output names, and explicit memory/disk reservations. Workers finalize
artifacts, but only the coordinator validates them and commits task markers.
This isolates allocator fragmentation and worker crashes and permits per-worker
cgroup controls. It never `fork`s a live OpenMP process. The parent admits the
sum of worker reservations under one byte-token budget with a 1/16 safety
margin, divides OpenMP threads among admitted children, validates exact output
counts and lengths, and kills sibling workers after a failure. Each completed
semantic range is committed independently, so a resumed generation restores
valid ranges and only respawns missing work. Cgroup-charged filesystem cache
still needs a deployment-level ceiling in addition to allocator byte tokens.

### Framed temporary streams and small-file policy

Large immutable path, rank, fixed-record join-run, and final-event streams can
use a version-2 framed format below their record API. Its 32-byte header identifies the codec and
logical block size. Each independently decodable block has a 40-byte header,
logical record count, raw/stored lengths, and checksum. A 16-byte-per-block
physical/logical offset index and a 56-byte footer provide counts, total logical
bytes, and a whole-stream checksum. The writer spools that index to one
construction-only sidecar and appends it at commit, so neither the number of
blocks nor random-seek metadata grows in RAM or remains as thousands of small
files.

All `PathNode` pointers, rank offsets, and join partition boundaries remain
**logical uncompressed offsets or record ordinals**. Readers binary-search the
on-disk footer index, decode one admitted block, and can mix legacy raw shards
or join runs with framed streams in the same construction. A framed join run
contains the byte-identical version-1 logical header, records, and footer, so
compression does not change its payload checksum or worker range ABI. `auto`
uses zstd level 1 when its workspace fits and stores an
incompressible block raw; under a very small memory budget it falls back to a
raw stream. A codec's workspace includes the decoder's, not just the encoder's:
generated path/rank shards are bounded so the next `prune()` (and, after the
last step, `MergedGraph`) can hold every pair it may open inside its own merge
budget, and join runs are bounded so `planJoinPartitions()` can admit the run
and sidecar readers it will need. Both bounds are applied before a stream is
written, because a block that only the producer can afford commits a generation
its consumer cannot read. A generation committed before those bounds existed
keeps its own block: the merge buffer is then raised to hold one decoded pair of
it, and a workspace whose block needs more than a sixteenth of `--memory-limit`
is refused with the limit it does need rather than silently over-allocated. Explicit `zstd` instead reports the minimum codec workspace it
cannot admit. Compression mode, block size, level, and worker count are
operational settings and may change on resume. Join child processes clamp
zstd's internal worker count to their assigned thread share, preventing
process-level and codec-level parallelism from multiplying unnoticed.
Join group-summary and one-byte detail sidecars use the same random-readable
framing. Their inner version-2 headers declare the record count up front and
leave terminal checksum/group fields for the authenticated footer, because a
streamed compressed first block cannot be patched in place. Readers retain raw
version-1 compatibility and continue to address summaries and details by
logical ordinal. Primary-run and sidecar stored/logical ratios are reported
separately so nearly-unique join-key workloads remain visible in telemetry.

A retained chr21 workspace inventory explains why the implementation batches
large streams instead of combining task markers into a general small-file
container: 534 files below 1 MiB occupied only 18.1 MiB, and 498 files below
64 KiB occupied 1.25 MiB, in an approximately 260 GiB workspace. Completion
records remain separate because their atomic rename is the task commit point.
Logical join bins are already ranges in shared pack files, while compressed
data blocks and their offset index travel in one artifact. Preliminary
`zstd -1` samples of the old workspace showed roughly 4--5x reduction for path
streams, 14--77x for generation rank streams, about 32x for join-partition
ranks, and 6--12x for final event streams. Those ratios motivate the codec but
are not a replacement for a new end-to-end chromosome benchmark.

An interrupted chr21 k32 diagnostic after the fixed-record codec landed made
the next bottleneck measurable. At 7,864,288 distributed records, 849,343,104
logical primary bytes occupied 114,797,355 stored bytes (0.135 stored/logical),
while 7,200,595 nearly-unique key groups produced 353,493,032 raw sidecar bytes.
The process was waiting in `balance_dirty_pages` and used only a small fraction
of one CPU during the final sample. The run was deliberately stopped before a
checkpoint; these measurements motivated framing the sidecars and are not an
end-to-end performance result.

The largest remaining I/O opportunity is further record compaction, not extra
writer processes. Framing now compresses fixed-width join runs and their exact
group/detail sidecars, but logical join records still repeat complete path
labels and may be larger than the variable-width source path/rank pair. A
versioned grouped codec writes one left context for consecutive expansions in
each label-sort run, followed by compact right-dependent records that reference
it. This removes
14 bytes of repeated `PathNode` context per grouped reference in addition to
label-prefix compression. A future join-run reference format could also avoid
repeating the complete rank payload before label sorting. Asynchronous
writeback (`sync_file_range` where available, followed by `fdatasync` only at a
commit boundary) may later overlap compute and flushing, but is not yet used:
the portable implementation currently syncs every 512 MiB before evicting that
prefix.

The chr20 trace makes the scheduling priorities concrete. Join-run encoding
occasionally used one CPU at about 39 MiB/s while the device was not saturated;
parallel block encoding/checksumming can help there. The label merge sustained
roughly 240--280 MiB/s and was storage-bound, so additional writers would not.
The final `MergedGraph` writer sustained roughly 100--120 MiB/s but retained
about 18.3 GB of clean staging data in cache; rolling cache eviction is the
first fix. The final index scan consumed approximately one CPU and 128 MiB/s of
cached logical input with no physical reads. It cannot be naively split because
`prev_occ`, the suffix-tree stack, and redundant-pointer aggregation are
ordered state;
the event/component split is what exposes safe parallel work around that
ordered core.

Generated paths and sorted bypass paths now feed the label-run builder directly
during the final join scan. The previous route materialized an approximately
61.5 GB step-4 path/rank pair and then decoded it into 67.1 GiB of fixed-width
sort records. The fused route removes one complete write/read cycle and reports
the exact avoided bytes for each generation. The v2 label-run codec additionally
groups adjacent expansions that share `from`, predecessor mask, and order; its
counters report both grouped references and payload bytes avoided. Rank payloads
in the earlier join runs remain the next compaction target. That compaction
remains more valuable than running two encoders against the same disk.

The current disk guard is exact for generated path/label sinks and deterministic
join-partition output: it counts committed and pending generation bytes against
`--disk-limit` and never compares those bytes with RAM. Those bytes are physical
throughout. A generation records the installed `st_size` total it was charged,
and the budget handed to its successor subtracts that value rather than the
logical `PathGraph::bytes()` payload; the subtraction saturates, so a frontier
that exceeds the limit produces an explicit refusal instead of an unsigned wrap
that would leave every later check vacuous. Framed output is admitted
using a conservative physical peak that includes per-block metadata and the
temporary index copy present during commit; completed and restored shards then
replace that reservation with their exact `st_size` totals. The `DiskBudget`
primitive can also inspect recursively live workspace bytes and `statvfs()`
free space, but it is not yet shared by every writer. In particular, final-event
copies, LCP raw levels, transient join-distribution runs/sidecars, and some
checkpoint copies currently stop on a failed write/ENOSPC rather than reserving
their complete output first. Consequently
`--disk-limit` is presently a spill-generation budget, not a strict global
workspace quota. A generation estimate remains diagnostic; RAM size is never a
disk-limit proxy.

## Durable workspace

The target persistent semantic organization is:

```text
build.json
inputs/
keys/
initial/
step-01/{prune,extend,join,label-sort}/
step-02/{prune,extend,join,label-sort}/
merged/
gcsa-components/
lcp-levels/
tasks/
logs/
```

The current slice commits `inputs/` plus semantic phase artifacts and
`*.complete` task records directly under the workspace root. Join and label
sort runs are immutable while in use, but are still PID-named transient files
inside a phase. Completed join-range outputs are durable semantic tasks keyed
by logical input and exact left/right half-open ranges, enabling
partition-granularity resume. Moving transient distribution and label runs into
the directory hierarchy remains format work.

`build.json` contains the workspace format, GCSA2 version, final format
versions, input identities/sizes/checksums, mapping identity, k-mer length,
doubling steps, sample period, LCP branching, alphabet, and logical input IDs.
Those fields form the semantic fingerprint. Memory, threads, run size,
partition target, merge fan-in, open-file limit, and I/O buffer size are stored
for provenance but excluded from the fingerprint. Temporary compression mode,
block size, level, worker count, and obsolete-artifact cleanup are likewise
operational: a resumed build may use a different codec or RAM/disk tradeoff for
new artifacts without invalidating already committed raw or framed inputs.

Durable phase artifacts have explicit little-endian headers and footers with
magic, format version, artifact kind, logical ID, physical shard ID, record and
payload byte counts, checksum, sort order, and key range. Join runs also use an
explicit versioned encoding. Transient label-sort runs use an explicit
little-endian, versioned grouped and prefix-compressed format. A group header
stores the left-path context once; following records store a reference flag and
only right-dependent node fields. The first label is complete and following
labels store LCP length plus suffix. They remain transient rather than workspace
tasks, but no longer dump an unstable C++ object layout. The final public index
keeps its unchanged legacy format.

### Commit and recovery protocol

1. Write and checksum an artifact incrementally at a unique `.partial` path.
2. Rewrite/finalize its header and append its footer.
3. Flush, `fdatasync()`, close, rename to the semantic final name, and sync its
   parent directory.
4. After every output artifact for a task exists, validate their metadata and
   atomically commit its semantic `*.complete` record using the same sync
   protocol. A later layout may place those records under `tasks/`.
5. Only a task record makes its outputs part of the restart frontier.

On resume, partial files and artifacts without a committed task are removed or
quarantined; immutable inputs and committed predecessors remain. Task records
are replayed in dependency order. Headers, lengths, footers, and task metadata
are always checked. Join-run payload checksums are selected by
`--verify-workspace`; normal construction trusts the checksum computed while a
newly synced run was written and avoids immediately rereading it. A newly
closed phase path/rank writer passes its durable physical checksum to checkpoint
adoption with the immutable fixed-record payload. Adoption therefore hard-links
the raw or framed file without a checksum reread when source and workspace share
a filesystem, and otherwise validates while performing one bounded copy.
Version-3 checkpoint
metadata records both logical record bytes and physical stored bytes and still
restores version-2 raw checkpoints. Normal same-filesystem restore validates
task identity, record count, and byte length and then hard-links the committed
inode into the next `PathGraph`; `--verify-workspace` additionally rereads the
complete checksum. Join-range worker outputs and final-event streams use the
same immutable opaque-payload protocol, avoiding both checkpoint rewrite and
normal resume copy.
The final-event writer incrementally hashes BWT masks, edge destinations,
samples, and occurrence events while producing them. Raw same-filesystem
adoption reuses those logical-stream digests. A framed stream is checksummed as
stored during adoption because its physical header/index/payload differs from
the logical records. The externally sorted redundancy stream remains raw and
retains a full checksum pass; cross-filesystem fallback always validates while
copying.
Final events retain a checked header/footer compatibility restore for older
workspaces. Cross-filesystem restore copies and checksum-validates in one pass.
A task with a missing or truncated output is invalidated and rerun along with
its dependants. `--verify-workspace` additionally discovers same-length payload
corruption by replaying the committed checksum; ordinary resume deliberately
trusts the checksum recorded by the closed writer. Cleanup first validates
predecessor and successor markers, atomically publishes a `*.retired`
deletion-intent journal, and only then unlinks predecessor artifacts not
referenced by another task. Recovery replays journals idempotently, including a
crash after any individual unlink. Once those deletions and their directory are
synced, the journal is atomically replaced by a `complete=1` audit record;
recovery does not revalidate an old immediate successor that was itself retired
later. A completed extend checkpoint also retires all of its fine-grained
join-partition and MSD-plan families. Enabling cleanup only on resume catches up
all generations older than the restored frontier. Cleanup never removes the
newest committed predecessor of an incomplete task.

Crash tests use deterministic process termination, without stack unwinding,
immediately before and after artifact rename, immediately before and after task
marker rename, and after one cleanup removal. A fresh process then reconstructs
the frontier. Renamed-but-unmarked artifacts are discarded; a marker that has
crossed its rename boundary keeps every referenced immutable artifact.

## External algorithms

### Keys and start nodes

Key extraction writes bounded sorted runs. Equal labels are combined inside
each run and again during a deterministic k-way merge, OR-ing predecessor and
successor masks exactly as `Key::merge()` does. The final merge assigns label
ranks while streaming and emits label/rank, last-character, k-mer-LCP leaf, and
de Bruijn count streams.

Start nodes use the same run/merge framework after optional node mapping. The
final deduplicated stream is consumed to build sparse/rank support without a
global node vector.

For `N` records and memory for `M` records, comparison-run construction costs
`O(N log M)` CPU and `O(N * (1 + merge_passes))` sequential I/O. Fixed-width
keys may instead use stable MSD radix distribution for linear work per visited
digit.

### Run-set PathGraph and label ordering

A `PathGraph` is an in-process run set of immutable shards carrying distinct
logical and physical IDs plus path/rank counts. Workspace checkpoints add
versioned artifact lengths and checksums. Initial extraction and join workers
may create many physical shards for one logical input; a bounded post-join sort
compacts those runs, and pruning can directly emit one label-sorted shard per
logical input because it already observes one global label-order stream.
The global label merger keeps only an LRU-bounded set of path/rank descriptor
pairs. Each open stream receives a 4--64 KiB sequential read window whose total
capacity is derived from, and capped at one quarter of, the merge group budget.
This preserves the descriptor and RAM ceilings while amortizing the 24-byte
path records and short rank payloads that would otherwise require one `pread()`
each. Budgets too small for one page deliberately fall back to direct reads.

Generated paths enter a byte-bounded buffer and are written as deterministic
label-sorted runs. The primary comparison is the existing first-label/rank
sequence order. Equivalent primary keys use semantic tie-breakers (`logical`,
`from`, `to`, predecessor mask, full interval, and input ordinal), making output
independent of thread scheduling. Multi-pass k-way merge costs sequential
`O(G * (1 + ceil(log_F R)))` I/O for `G` generated bytes, `R` runs, and fan-in
`F`.

### Prefix-doubling join

Each doubling generation performs a self-join within a logical input:

```text
left key  = (logical_file, path.to)
right key = (logical_file, path.from)
```

The planner takes a deterministic systematic sample of at most 4096 records
from each side, estimates both input bytes and matching-key fanout, and
recursively refines oversized ranges by four most-significant key bits. The
result is a canonical set of explicit prefix packs with estimated input/output
bytes plus suspected heavy keys. Pack memory is capped from the configured RAM
budget; if that cap is reached, a coarser pack is retained and the exact pass
continues splitting it. The plan is bound to both join-run payload checksums,
committed as a checksummed workspace artifact, and restored independently of
partition outputs. Each sorted join run emits an exact compact sidecar with one
fixed-width summary per join key and one byte of order/sorted detail per source
record. Planning and pathological-key subdivision read those sidecars instead
of rescanning the wide variable-rank `JoinRecord` payload. The plan is bound to
the run and sidecar checksums; completed partition task names carry the same
identities, so equal-length replacement runs cannot reuse stale output. The
already-sorted left/right join files act as shared pack files:
thousands of logical radix bins are half-open offsets in two files instead of
thousands of open descriptors. Sidecar summaries remain authoritative for every
range, count, and byte estimate, so a sampling miss can change task balance but
cannot change join semantics.

For each partition, the scheduler reserves a declared working set. If one side
fits, it is indexed and the other is streamed. Otherwise the partition is
recursively split on more key bits. A single unsplittable heavy key is processed
as bounded left blocks crossed with bounded right blocks; generated records are
immediately sent to label-run construction. This uses
`O(|L_key| * |R_key|)` generation work, as required by the output, but only
`O(block_L + block_R + output_buffer)` RAM.

Sorted paths bypass the join and enter label-run generation unchanged. After a
bounded chunk is label-sorted, consecutive records with the same left context
are encoded as one group header followed by compact references. The reader
reconstructs full `PathNode` values only as the merge consumes them. Complete
rank payloads are still generated before label sorting, preserving the exact
legacy comparator and keeping this optimization local to the transient format.

The production route uses an external sort-merge implementation. For every
logical input it writes explicit, versioned left-by-`to` and right-by-`from`
run records. Run creation uses leveled compaction, so file-name metadata is
`O(F log_F R)` for merge fan-in `F` and `R` initial runs. A matching key is
processed without materializing either key group: it loads bounded blocks from
the smaller side and replays disk-backed blocks/ranges from the other side.
This already handles an individual key larger than RAM and preserves joins
across physical shards. Generated records flow directly to grouped,
prefix-compressed bounded label-sort runs. The run format avoids repeated left
`PathNode` context; the complete generated rank label remains present as a
prefix-compressed sequence because it is the primary external-sort key.

### Pruning

`SpillableGroup` and `SpillableDeque` keep only bounded resident windows for
equal-label records and extended ranges, spilling excess records to sequential
files. `SameFromLogicalFile` computes the required all-from/all-logical summary
through bounded record access. Output branches replay spilled records only when
individual paths must be emitted. `MergedGraph` likewise externally
sorts/deduplicates its potentially huge same-from sets. Physical shard IDs never
participate in `SameFromFile` semantics.

### Final components and LCP

The implemented final merged-graph scan emits immutable BWT masks,
per-character source-path edge ranks, sampled path positions, sample IDs and
boundaries, nonzero occurrence pairs, and redundancy positions. Redundancy
positions are externally sorted. Dense previous-occurrence state and the
suffix-tree traversal stack are block-cached disk arrays. Mapped start nodes for
one path use two reusable `SpillableNodeSet` instances: small sets stay in RAM,
while oversized sets are externally sorted/deduplicated and rewound from disk
for occurrence, continuation, and sample passes. The asynchronous merged-graph
readers use byte-sized foreground/prefetch windows and explicitly release
consumed kernel-cache ranges, so widely separated cursors do not rely on cgroup
reclaim to bound clean cache. All streams are synced and committed together as
one workspace task; a crash restarts only this ordered scan. Mid-scan resume and
the assignment-log protocol remain future work.

Final raw components are built serially, one at a time, from those streams.
The library retains the resident-object constructor for compatibility, while
`GCSA::buildAndStore()` writes each fast/sparse BWT, rank support, edge, sample,
occurrence, and redundancy component immediately in the exact existing
`GCSA::load()` order and releases it before constructing the next component.
Fast BWT, edge-boundary, and sampled-path `bit_vector_il<512>` members are
serialized directly from their event streams, including SDSL's interleaved
cumulative counts and bounded select samples. They never materialize the dense
source bitvector or the full interleaved destination; empty `rank_support_il`
payloads are emitted without rebuilding the corresponding vector. Packed
sample IDs are likewise serialized word by word in the existing
`sdsl::int_vector<0>` layout instead of allocating one entry per sample.
Sparse BWT components replay the compact path-mask stream four times to write
their Elias--Fano low bits, unary high bits, and two high-bit select supports
directly. Their `rank_support_sd` objects have no serialized payload and bind
to the vectors at load time, so they require neither another scan nor another
resident sparse vector.
Sample boundaries replay the monotone endpoint stream twice to write the dense
bit payload and `select_support_mcl` payload directly. The latter is spooled in
4096-one blocks, matching SDSL's slow and fast initialization layouts without
retaining either a bit per sample ID or all select blocks in RAM.
Occurrence pointers replay their sorted `(path, extra-count)` events to write
both Elias--Fano `sd_vector` members directly: packed low bits, unary high bits,
and high-bit select supports. Select directories and block payloads spill to
the workspace, so no builder or directory grows with the occurrence count.
This currently uses eight sequential event passes in exchange for the bounded
RSS invariant; a fused multi-output encoder is a later bandwidth optimization.
Sorted redundancy events likewise stream the exact ordinary `sdsl::bit_vector`
and `select_support_mcl` payload used by `SadaCount`; its select payload is
spooled in 4096-one blocks, so direct packing retains no dense unary vector.
The partial file is synced and atomically published. Standalone `build_gcsa`,
`vg index`, and `vg autoindex` use this direct packer whenever a workspace is
configured; they reload the index only for explicit verification. This removes
the cumulative completed-index peak and the final component-sized SDSL builder
floor from the production file-building route.

LCP level `i+1` is produced by streaming groups of the configured branching
factor from level `i`; each raw internal level is independently checkpointed
and consumed. `LCPArray::buildAndStore()` then writes the packed hierarchy
directly in legacy leaf-to-root order, including SDSL's historical partial-word
padding, and atomically publishes the final `.lcp`. The production standalone
and `vg` routes therefore do not allocate `LCPArray::data`; they reload it only
when explicit verification needs the query object. The resident constructor is
retained for API compatibility and small callers.

Index verification also has a disk-first route. It emits fixed-width
`(label, mapped-from-node)` records from bounded input blocks, externally sorts
and deduplicates the expected occurrences, streams raw `locate()` occurrences
through a callback, externally reduces those records, and compares the two
sorted streams. The verification workspace defaults to at most 64 MiB even
when construction has a larger budget. It is currently serial and transient;
checkpointing its sort runs and partitioning disjoint label ranges are
throughput improvements rather than verification-correctness requirements.

## Implementation ledger

The ledger is intentionally conservative: a feature moves to **complete** only
when its production call path and forced-spill/recovery tests pass.

### Historical parallel-finalization checkpoint (2026-09-29)

At this checkpoint, items 1, 2, and 4 were
implemented and validated on chr21. The final candidate completed in **38:06
versus 41:07** (7.34% shorter), with byte-identical GCSA/LCP outputs. Final merge
was 424.4 versus 468.4 s; final events 578.1 versus 666.3 s; component encoding
71.7 versus 155.9 s. Average cores increased 3.14 to 3.76; time below two cores
fell 74.7% to 59.7%; peak sampled RSS remained about 14.6 GiB. Candidate total
cgroup peak was 27.498 GiB, with no memory-limit/OOM events. This is a single
shared-host fixture comparison, not a whole-genome performance prediction.

All 15 GCSA2 targets, all 35 vg integration checks, and the expanded 4,097-pattern
exact-output/resume fixture pass. Regression tests cover descriptor-limited
cache thrashing and backward buffer refills discovered in two stopped attempts.
The sealed candidate is `vg-gcsa-final-streams-20260929`, SHA256
`3a8fffc210c9c0c541efa73d6b9b294ba622167e4808ad05550323bef155e91c`, based on
vg `fd2ac2ba4` / GCSA2 `88d0621` plus the retained source patch. Ordinary `bin/vg`
is unchanged. Source edits are uncommitted. At that checkpoint, the next
recommended gate was a GCSA2-only comparison on retained chr2+chr18 inputs;
it was not launched.

Detailed implementation, receipts, CPU plot, resource tradeoffs, and preserved
attempt history: [final parallel result](/mnt/ssd/lalli/hprc_v2_vg_rna/gcsa2_distribution_concurrency_20260928/final_parallel/README.md).

### Completed serial-pass chr21 gate (2026-09-30)

Fresh `run_1` exited 0 and reproduced final-parallel `run_4`'s public MD5s,
`297f5fe840f7ec48337c97eacdc00c0f` (GCSA) and
`6fc46b008b54f34cc30dcbb812807a9a` (LCP). All 16 GCSA2 suite targets and all
35 vg integration checks pass, including the regression that keeps parallel
prune off the legacy construction route. The sealed candidate is
`bin/vg-gcsa-serial-passes-v2-20260930`, SHA256
`0d5016aacc95a48705b798b25c3144f3caffa481d0387212562ad5d4210b6a92`.

The three bounded changes exercised their intended paths. Four pruning rounds
used seven workers and seven root partitions each, with no parallel fallback;
their reported wall times totaled 157.2 s. Parallel prune is restricted to the
external-memory route because its extend implementation preserves logical input
identity across physical shards; the legacy extend implementation joins only
within each physical file and therefore retains serial prune. All nine main
initial/prune/extend path/rank checkpoints performed zero adoption checksum
rescan and reused about 31.7 GiB of checksum provenance from their closed
writers.

Final preparation processed 399,778,113 paths and 631,824,030 ranks with seven
workers and no spill fallback. Admission reserved 128 MiB for the bounded
node/rank arenas plus 64 MiB for predecessor batches, 192 MiB combined. Workers
prepare independent node/rank inputs; the previous-occurrence state, suffix-tree
stack, and event serialization remain ordered. The complete final scan therefore
retains serial state transitions.

The candidate took 2,363.3 s (39:23), versus 2,285.5 s (38:06) for `run_4`, a
3.4% wall-time increase. Candidate versus reference phase times were
preprocessing 369.731 versus 151.342 s, doubling 953.1 versus 1,021.3 s, merge
409.1 versus 424.4 s, events 500.3 versus 578.1 s, and components 76.5 versus
71.7 s. Average cores were 3.77 versus 3.76; time below two cores was 47.5%
versus 59.7%; sampled maximum RSS was 14.46 versus 14.57 GiB; and cgroup peak
was 28.125 versus 27.498 GiB, with no memory-limit or OOM event.

The main-thread sampled CPU profile supports the targeted checksum improvement:
the raw checkpoint checksum function fell from 97.4 to 10.6 sampled CPU seconds,
while `ExternalInputPreprocessor::buildLCP` was 76.6 to 74.8 s and `prepare`
73.6 to 71.8 s. These values are inclusive CPU samples rather than wall times;
overlapping values must not be subtracted, and they do not establish the cause
of the preprocessing wall increase. This single shared-host A/B establishes
exactness and the targeted mechanism behavior. Total wall time did not improve.

The isolated triplet-state prototypes are now implemented. All 16 retained
synthetic end-to-end runs have byte-identical public outputs; triplet pruning
used 24 workers across 65 certified partitions. Two large-fixture repetitions
observed 6.09 versus 10.27 s mean pruning time for triplet-only versus default
(40.66% lower), not an established whole-build or chr21 speedup. Parallel final
state is exact, but its 24-worker execution replayed 98.13M LCP updates and has
no established speed benefit. The production GCSA2 implementation and ordinary
`bin/vg` are unchanged, with no promotion or chr21/chr2+chr18/whole-genome
launch. All 16 GCSA2 suite targets and final medium default/combined query
smokes passed. The earlier tie-order assertion and interrupted suite are
retained diagnostic history: a test-only repair canonicalizes full records
within identical-label groups and compares downstream streams; no core algorithm
change was required. The recommended real-chr21 triplet-only measurement remains
unstarted; the state-chunk result below is completed historical evidence. The evidence navigator is
`/mnt/ssd/lalli/hprc_v2_vg_rna/gcsa2_distribution_concurrency_20260928/triplet_state_20260930/README.md`.

The state-chunk trial is complete:
`/mnt/ssd/lalli/hprc_v2_vg_rna/gcsa2_distribution_concurrency_20260928/state_chunks_20260930/README.md`.
It replaces repeated per-worker LCP replay with contiguous path chunks, using
ordered-pass canonical stack seeds (at most 257 frames), rank-owned
previous-occurrence updates, chunk-local queries, and coordinator emission. The
sealed binary SHA256 is
`21b90225129697d8785dbbe78f2e8c8454f541912d7fd266faa4c254375aa2a4`.
All 16 GCSA2 targets and focused 1/2/24-worker, tie, spill, mapping-dedup,
state-arena/predecessor-fallback, and resume checks passed. All 11 builds exited
0 and byte-matched the archived public pair; five medium builds independently
query-verified. On large input, worker replay fell 98,131,296 to 4,194,305
(95.7258%) for 4,194,305 planned paths in 249 batches. Fresh old/new/new/old
24-worker events averaged 8.75 s in both arms; build means were 86.1497 versus
90.08535 s, so no wall-speed gain is established. New-24 chunk queries averaged
0.0201025 s, but serial emission averaged 3.0748785 s. This is a work-efficiency
success, not a performance promotion. **Completed emission profiling** used
sealed instrumentation binary SHA256
`032045ae8efc7592e0e15b63c82123a18b32add572001eb52c64a82fb081caa4`;
all 10 runs passed, five medium runs query-verified, seven runs profiled, and
the full 16-target suite passed. On the two 24-worker runs, emission was
3.735008921/2.910755412 s and `writeAll()` was 3.697841370/2.873472212 s
(about 98.88% on average); CPU outside full-buffer flush was only
0.037077590/0.037206828 s. Close `fdatasync()` was 0.001495699/0.002553545 s,
with no periodic sync. This bounds the dominant scope to full-buffer writes,
but does not establish a kernel writeback or filesystem-contention cause and
does not establish speedup. Aggregate telemetry makes native CPU sampling
unnecessary. The next recommendation is a bounded write-stall/overlap
investigation, not yet implemented or authorized. Production implementation and
ordinary `vg` remain unchanged; no relink, promotion, chromosome, or
whole-genome launch is authorized. Evidence is at
`/mnt/ssd/lalli/hprc_v2_vg_rna/gcsa2_distribution_concurrency_20260928/emission_profile_20260930/README.md`.

| Slice | State |
| --- | --- |
| Shared byte-token `MemoryBudget` and statvfs-aware `DiskBudget` | memory primitives, major external-phase byte caps, and framed encoder/decoder/zstd-worker reservations implemented; pruning/merged-graph buffers are byte-bounded but not globally admitted; final-scan readers retire cgroup-charged cache without extra descriptors; mapping/input/library allocations remain outside the shared budget; disk guard uses conservative framed peaks and exact installed sizes for generated path/join volume but not every final/LCP/checkpoint writer |
| Versioned workspace artifacts, atomic publication, manifest compatibility, recovery tests | phase and join-range checkpoints plus abrupt-exit commit-boundary tests implemented; newly closed path/rank writers provide durable physical checksum provenance so normal same-filesystem adoption does not reread their payloads; the chr21 gate reused about 31.7 GiB across nine main checkpoints with zero adoption rescan; raw/framed mixed restore and journaled idempotent predecessor/family retirement implemented; transient distribution-run task records pending |
| Human-readable operational construction parameters | core budget/run, temporary-compression, process-worker, and safe-cleanup controls implemented and parser-tested; checkpoint cadence remains standalone-only |
| Framed temporary-file compression | versioned independently checksummed zstd/raw blocks, disk-spooled footer index, logical random access, raw fallback, mixed raw/framed join reading, exact primary/sidecar size telemetry, shared-budget admission, and multithreaded zstd contexts implemented for path/rank, fixed-record join-run plus group/detail sidecars, and final-event streams; a block is now admitted against its consumer's decode budget as well as its producer's, and the multi-worker estimate charges the job pool zstd actually allocates while clamping workers a block cannot feed; redundancy, preprocessing, and LCP level files remain raw |
| Bounded external path-label runs and leveled multi-pass merge | versioned grouped/prefix-compressed runs, shared-left-context references, parallel in-place run sorting, rolling cache windows, one-run bypass, direct leveled merge of already-sorted worker shards, and forced multi-pass spilling implemented and tested |
| Logical/physical `PathGraph` identity in pruning and joining | implemented; pruning coalesces only by logical ID and never treats a physical shard as a semantic input; durable run-set manifest pending |
| External prefix-doubling join | bounded sort-merge, rolling page-cache windows, compact exact key-group sidecars, RAM-capped persisted 4096-record MSD radix range-pack planning bound to run/sidecar checksums, exact range contracts, recursive two-dimensional heavy-key splitting, blocked nested-loop expansion, and selectable checksum scans implemented and forced-spill tested |
| Direct join-to-label pipeline | implemented and forced-spill tested; the unsorted generated path/rank pair is no longer materialized |
| Process worker scheduler and partition resume | implemented with fork-free spawn, global byte admission, immutable-payload semantic range checkpoints, same-filesystem zero-copy restore, and one-/multi-partition tests; a restored partition is admitted on the workspace's own record and charged its exact stored size instead of a peak re-derived from the resuming run's codec, and a killed sibling's pid-named label-sort runs are swept, including those of the worker whose failure ends the run |
| External keys/start nodes/initial paths | implemented with bounded fixed-record runs, global duplicate reduction, streaming support construction, durable key/start checkpoints, and physical initial shards preserving logical IDs; only key LCP support remains resident during doubling, while mapper, last-character, and start-node supports are delayed to merge/final scan |
| Spillable pruning groups | implemented and forced-spill tested for equal-label priority groups, extended ranges, external same-from sets, and bounded input/output descriptor caches; the external route can partition zero-LCP roots across admitted workers and reconcile the exact logical output in root order, with resource/safe-boundary fallback to serial; legacy construction remains serial because its extend path is not physical-shard aware; framed input pairs are bounded at write time so every shard a merge opens stays resident and each block is decoded once, and a pair committed above the current merge budget is admitted with its overshoot reported rather than refused, because a committed block cannot be renegotiated; the descriptor allowance is divided by what each cache can occupy rather than evenly, since pruning emits one output per logical input while only the input side scales with shard count, and the post-join compactor retains against that same bound (both label merges visit their shards in round-robin label order, where a cache one entry short of the shard count misses on every record rather than on a fraction of them) |
| Final event/component passes | implemented with one-task immutable-payload event checkpoint, same-filesystem zero-copy restore, spillable mapped start-node sets, disk-backed `prev_occ` and suffix-tree stack, external redundancy sort, streaming fast/sparse BWT, packed-sample-ID, sample-boundary, SadaSparse, and SadaCount serialization, direct component-at-a-time packing, bounded concurrent component queues, and predecessor batches preserving ordered state without rereading the source; bounded workers prepare independent node/rank inputs while `prev_occ`, the suffix stack, and serialization remain ordered; all 16 GCSA2 suite targets and all 35 vg integration checks pass, including legacy/parallel 4097-pattern resume, sparse/skewed and batch-boundary fixtures, spill/fallback coverage, and decoded-record equality for the changed-I/O-buffer fallback |
| Streaming LCP levels and direct packing | implemented and resume-tested with one raw level resident at a time, bounded direct final serialization, atomic publication, and byte-identical legacy output including padding edge cases |
| External verification | implemented and forced-spill tested with bounded input blocks, external expected/actual occurrence sorts, callback-based locate, and sequential set comparison; resumable runs and parallel label ranges pending |
| Standalone and `vg index` / `vg autoindex` CLI integration | implemented for memory/disk, process workers, temporary compression, and safe cleanup; forced-spill/resume and raw-versus-framed equivalence integration tested |
| Chromosome-scale benchmark | preexisting chr20 k32-pruned fixture completed and verified at 32.36 GiB RSS under a 128 GiB cgroup with byte-identical legacy outputs; preexisting chr21 k32-pruned fixture completed and verified at 21.34 GiB process RSS under a hard 25 GiB cgroup, including a successful step-4 resume; the 2026-09-30 chr21 serial-pass gate reproduced the accepted public MD5s but took 2,363.3 versus 2,285.5 s, so it is correctness and targeted-mechanism evidence rather than a whole-run speedup |

Current limitations are intentionally explicit. The external route is selected
only when `ConstructionParameters::work_directory` is nonempty. The production
file-building routes no longer retain all completed GCSA components or retain a
completed GCSA while constructing LCP. Fast BWT/edge/sample vectors, packed
sample IDs, sample boundaries, and redundancy pointers stream in bounded
memory, and sparse BWT and occurrence pointers use bounded multi-pass
Elias--Fano packing. The resident-object compatibility constructor still uses
SDSL component builders; production file-building routes do not. Input graph
decoding, mapping support, allocator overhead, and charged page cache still
keep the in-process budget from being a strict whole-process RSS ceiling, so
production runs retain a cgroup limit. The de Bruijn mapper is no longer kept
alive throughout prefix doubling, and the last-character/from-node supports
are not created until the final event scan, but those completed succinct
objects are still resident rather than disk-cached while active.

Join distribution runs and their exact group sidecars are checksummed immutable
files but are not yet registered as resumable workspace tasks, so an incomplete
generation rebuilds them before restoring the persisted MSD plan and completed
join ranges. Initial records and join distribution records still repeat full
rank payloads.
Framed compression currently targets the dominant committed path/rank,
fixed-record join distribution and its group/detail sidecars, and final-event
streams. Key/start preprocessing runs, transient grouped label runs, the externally
sorted redundancy stream, and LCP level files remain in their existing
compact/raw formats. Framed final-event
adoption also performs one physical checksum pass after close; exposing the
codec's stored-byte digest to the workspace would remove that pass. A configured
compressed block establishes the minimum decoder reservation for any later
phase that reads it, so very small resume budgets must still fit one block or
receive a clear refusal.
Label-sort runs group shared left contexts and prefix-compress labels, but a
deeper reference representation spanning the pre-sort join stream remains an
optional compaction rather than a feasibility dependency.
Final component construction, final-path merging, final-event preparation, and
external-route pruning now have bounded parallel paths; the completed full-suite
and exact chr21 gate is recorded above. Ordered previous-occurrence and suffix
state remain serial, and legacy construction retains serial pruning because its
extend path cannot join across physical shards. Process workers accelerate
independent join ranges, and
post-join compaction runs the batches of a pass on up to `-t` threads when each
batch fits one output pair and the disk and memory budgets admit them together.
On the chr21 k32 fixture (2026-09-28, 200 GiB goal, `-t 24`) that took step-4
compaction from 682.8 s to 82.4 s (48 batches, 24 at once) with the reference
`.gcsa` and `.lcp` digests unchanged. When the budget admits two distribution
sorters the two join sides are scanned concurrently, but each side's sort,
record encoding and write happen in `ExternalJoinSorter::finish()`, which still
runs one side after the other; on the same fixture the concurrent scan left
distribution within run-to-run variation. That path also calls
`TempFile::getName()` from both scan threads whenever a side's buffer fills
mid-scan, and `TempFile` is not thread-safe. Even with these limitations, preprocessing and prefix doubling no
longer require a logical chromosome, physical shard, join key, equal-label prune
group, same-from set, or generated label run to fit in RAM, and LCP no longer
retains all raw hierarchy levels. One final path's mapped start-node set is also
no longer required to fit in RAM.

### Whole-pangenome feasibility and throughput priorities

`--memory-limit` in `build_gcsa` and `--gcsa-memory-limit` in `vg index` /
`vg autoindex` are the external working-set goalpost. Sort, join, worker, and
verification reservations are derived from this byte value; autoindex falls
back to its `--target-mem` value when no GCSA-specific value is given. Until
every library allocation is admitted, a deployment cgroup must remain the hard
whole-process ceiling, with the GCSA goal below it for allocator, SDSL, and
charged page-cache headroom.

Join distribution is an earlier, non-overlapping lifetime and can use the full
phase-local join budget while sharing the ceiling only with one source decoder
and its output codec. An explicit join-partition cap still bounds that phase;
automatic mode can use the otherwise idle aggregate workspace. During
join-to-label expansion, label sorting receives 75%
and join blocking 25% of the current process reservation. Process workers
recompute the same ratio after the aggregate goal has been divided among
admitted workers; they do not each inherit the full global target. Explicit
`--sort-run-size` and `--join-partition-size` values remain upper bounds and may
change on resume. Even when both caps equal the global memory limit, they do not
become competing reservations: the concurrent lifetime starts at 75/25. When a
cap or a tiny reservation makes that split infeasible, construction transfers
unused capacity to the other side while preserving both minima and the aggregate
ceiling. `build.json` records the nominal automatic shares; the phase allocator
applies this minimum-aware rebalance to the top-level task and again inside every
admitted process-worker reservation.

The goal is operational, not semantic, and may change on resume. A larger value
admits larger initial/sort runs, join blocks, and more worker reservations,
usually reducing physical shards, merge passes, and HDD traffic. A smaller
value deliberately produces more spill runs and I/O while preserving logical
file identity and final-index semantics. Disk limit and free-space safety—not
the RAM goal—decide whether a valid but very large generation may continue.

Following the completed emission profile and any explicitly authorized real-chr21
triplet-only measurement, a bounded matched diagnosis of the 369.731 versus
151.342 s preprocessing result comes before a chr2+chr18 gate. It must reuse
one identical retained k-mer stream across both binaries to remove
generated-order variability. The current shared-host comparison does not
establish the slowdown's cause, and the targeted sampled-CPU reductions cannot
be used to account for wall time.

Before a whole-pangenome claim, the remaining work is ordered by the project's
feasibility, RSS, then speed policy:

1. Bound the `vg` input producer. `SourceSinkOverlay` now discovers components
   on the fly and retains only visited state, its traversal stack, and component
   tips, including a sparse-ID fallback. `VGset::for_each()` still loads one
   whole physical graph before GCSA2's budget exists. A sharded graph reader
   must preserve one explicit logical GCSA2 input identity across all physical
   chunks.
2. Bound the remaining resident floor: `NodeMapping`, active de Bruijn support,
   input graph decoding, and library/allocator scratch allocate outside the
   shared token budget. Phase-lifetime scheduling now removes inactive key and
   start supports from prefix doubling, but external mapper/mapping lookup and a
   budget-aware `vg` input producer remain the principal greater-than-RAM
   changes. Until then, measured component floors must be reserved below the
   deployment hard cap.
3. Make `DiskBudget` global. Join/path generation is reserved today, but
   checkpoint payloads, final events, verification runs, and LCP levels can
   still discover ENOSPC only when a write fails.
4. Add block-level final-event recovery. The event scan is bounded but a crash
   currently restarts the complete scan; idempotent assignment logs and output
   blocks would retain finer progress.
5. Persist distribution runs and their already-emitted exact group sidecars as
   task artifacts. This avoids rebuilding both products after interruption.
6. Compact the repeated rank payload in join distribution records. Framed
   group/detail sidecars now remove the immediate raw-write bottleneck; delta
   coding remains an optional second-stage reduction before increasing worker
   count. The first-stage block codec reduces physical fixed-record bytes
   without changing the logical format; a later side-specific/reference format
   should reduce the logical bytes themselves. Only after per-path bytes and
   storage headroom are measured should more encoders, merge groups, or verifier
   label ranges run in parallel under combined memory and I/O admission.
7. Extend framing only where measurement justifies it: LCP levels and
   preprocessing ranks are promising, while already grouped/prefix-compressed
   label runs may not repay another codec layer. Path/rank checkpoint adoption
   already consumes its writer's durable checksum; exporting the final-event
   framed writer's physical checksum and fusing compatible occurrence component
   passes remain possible ways to remove other rereads.
8. Add a resource-aware read/compute/write pipeline. Join partitions already
   run as independent `posix_spawn` workers and each stream can use bounded zstd
   threads; final component encoders and verifier label ranges are the next
   safe process tasks. Their combined RAM reservations and device throughput
   must be admitted globally so concurrency stops increasing once the storage
   device is saturated.

Acceptance should report the configured GCSA goal and cgroup cap separately,
along with peak anonymous/cache/PSS memory, live disk, records and bytes per
generated path, merge levels, worker count, filesystem throughput/PSI, and
total `(bytes read + bytes written) / source byte`. Whole-genome runtime should
not be extrapolated until the projected resident floor fits the selected cap.

### The path merge is serial, and why it need not be

`MergedGraph::MergedGraph` (`src/path_graph.cpp:4052`) is the last single-
threaded phase of any size. Measured on the chr21 gate fixture at 96 threads, it
reduced 1,208,172,107 paths in about 26 minutes at 0.50 cores while the step-4
join partitioned the same data across 241 workers in 5m29s. The gap is not I/O:
the framed prefetch reported 6,238 of 6,279 blocks already resident, 41 waits,
zero synchronous fetches and zero errors, so reads are effectively free and the
phase is idle waiting on its own emit loop.

Four values are carried between iterations, and each is written into the output
rather than merely observed. `rank_count` becomes the pointer stored in every
`PathNode` via `setPointer`, so record N's stored value depends on the ranks of
records 0..N-1. `path_count` is the ordinal written into `from_file` records and
into the `next` transform. `from_count` feeds `next_from` the same way.
`path_lcp` is the LCP against the preceding record, which is pairwise-adjacent
across the whole stream. `curr_comp` looks like a fifth but is not: it is
monotone in `firstLabel(0)` and a partition can binary-search its own starting
value out of `next`.

None of these prevent partitioning, because all are prefix computations. Split
the sorted path space into P disjoint key ranges, merge each independently with
counters based at zero, prefix-sum the per-partition totals, and re-base the
stored pointers as the partition outputs are concatenated in key order. The
partition boundaries are exact rather than approximate: the merge is a stable
merge of sorted inputs, so for a split key K the number of paths preceding K is
the sum over input files of the count of records below K in that file, and
`PathGraphInputCache::read` is addressed by record ordinal rather than by byte,
so those per-file counts are a binary search over an existing API. Only the P-1
cross-boundary LCP values need the serial pass, computed from each partition's
first record.

This is the same shape the join already uses -- `sampleJoinKeys`
(`src/external_join.cpp:2498`) samples 8,192 keys into MSD range packs, with
`"MSD range packs do not tile the join-key space"` as an explicit invariant --
and the same shape `PhaseUnfolder` uses on the `vg` side, where workers mint
duplicate ids locally and a serial pass in component order re-bases them.

Before implementation, two costs required explicit accounting. Each
partition needs its own `PathGraphMerger`, whose `ranges` deque and `buffer`
group are both sized from `group_buffer_bytes`; P mergers therefore multiply the
merge reservation or divide it and spill more, and the phase allocator has no
term for that today. And the change moves pointer values in the output, so it
requires its own byte-identity gate against a serial run -- the chr21 fixture
digests `297f5fe840f7ec48337c97eacdc00c0f` and
`6fc46b008b54f34cc30dcbb812807a9a` are the reference, and both the `.gcsa` and
the `.lcp` must match, since a wrong re-base can leave a self-consistent index
whose pointers are uniformly shifted.

Implemented as the bounded group-aware parallel final-path merge described in
the implementation checkpoint above. The prefetch observation remains relevant:
the change targets serial emission rather than adding readers.

GCSA construction runs this merge serially unless
`GCSA_EXPERIMENTAL_PARALLEL_MERGE=1` is set. Copying every partition into place
doubles the merge's scratch writes, and on chr21 the parallel merge took 0:07:04
against 0:07:48 serial while two runs of the same parallel code differed by
0:00:42. `MergedGraph` still takes an explicit worker count, and its tests
exercise both routes.

### Descriptor budget and latency-tolerant I/O (2026-09-30 to 2026-10-01)

**The descriptor-limit abort.** A joint chr2+chr18 build at 64 threads (vg
`f28c31fba`, GCSA2 `4ab21ca`) aborted after 2:33:43 in step-1 compaction with
"compressed block: open failure" from `compactLogicalJoinShards`. Transient
systemd user services on this host start with a soft `RLIMIT_NOFILE` of 1,024,
while concurrent compaction batches were sized against the 128-file merge
ceiling (52 descriptors per chr2 batch, twenty batches wanting 1,040).
`609d7d9` bounds concurrent batches by the ceiling divided by (2w + 4) for
batch width w.

**A separate concurrent open-file budget.** The 128-file ceiling
(`--max-open-files`) still sizes every individual merge, because raising it
shrinks every committed compression block and reshapes the compaction target.
Work that runs concurrently is admitted against a second number,
`ConstructionParameters::concurrent_open_files` (`build_gcsa
--concurrent-open-files N`; default the soft `RLIMIT_NOFILE` minus the ceiling,
never below the ceiling; `getConcurrentOpenFiles()` in `src/support.cpp`).
`raiseOpenFileLimit()` lifts the soft limit to min(hard limit, 1,048,576) and is
called first by `build_gcsa`, `vg index` and vg's autoindex GCSA recipe (vg
branch `gcsa-descriptor-budget`). Raw-shard prune workers (`2fc9af5`) and
compaction batches and join-range planners (`69aa6b6`) are admitted against
it. Key-range distribution is not: `src/external_join.cpp` sets its descriptor
budget to the ceiling minus two and needs nine descriptors per range worker,
so at most 14 range workers run at once.

Measured effect on the first prune of the joint input (four resumes, binary
`69aa6b6`): the old admission ran it serially in 0:22:44.6 and 0:09:56.1; the
budget admitted seven workers, 0:04:27.7 and 0:05:57.2, with identical pruned
counts. The joint rebuild under the budget reproduced the vg-built baseline
index byte for byte (`joint.gcsa` MD5 `39441e40a12472d8e005a0c7926b0e90`, LCP
`f05d53df0d4d7b8aa5b71cd52f963b09`). Evidence:
`hprc_v2_vg_rna/gcsa2_distribution_concurrency_20260928/step1_prune_workers_20260930/`
and `joint_step4_resume_concurrent_20260930/`.

**Latency-tolerant I/O (`1c0e542`, `2fd1694`), opt-in.**
`GCSA_IO_DIRECT_WRITES=1` writes temporary streams with `O_DIRECT` through a
writer pool (`GCSA_IO_DIRECT_WRITE_THREADS`, default 16;
`GCSA_IO_DIRECT_WRITE_DEPTH`, default 4); tmpfs falls back to buffered writes
with a one-time warning. With it set, the final path merge's outputs
(`SequentialRecordWriter`) also bypass the page cache instead of syncing every
64 MiB (`2fd1694`). `GCSA_IO_READAHEAD_BYTES=64M` gives each sequential reader
a read-ahead window served by a shared pool (`GCSA_IO_READAHEAD_THREADS`,
default 16). The motivation is the host's 100 MB `vm.dirty_bytes`, shared by
every user: buffered writes stall in the dirty-page throttle, while reads are
mostly served from the page cache (in one 2:14:47 window, 1,555 GB of
application reads against 185 GB from the device). On a 97 MB synthetic graph
the median whole build fell from 0:01:45.0 to 0:00:57.7 with direct writes and
identical output; read-ahead alone showed no consistent effect. No controlled
chromosome-scale timing exists yet; a chr18 merge-write comparison is recorded
under `chr18_iteration_20261001/` in the evidence root: on chr18 alone the merge took
0:06:31.3 buffered (15 of 195 build-thread samples in the dirty-page throttle) against 0:05:35.7
direct (none of 168), one run each, with byte-identical final-events checkpoints.

### Prefix-depth prune partitions and stitching (2026-10-01, experimental)

The parallel prune splits its input where the first key character changes, so it has at most
seven partitions. `GCSA_EXPERIMENTAL_PRUNE_SPLIT_DEPTH=d` (`d958a7d`) splits wherever adjacent keys
share fewer than `d` characters. Every comparison inside such a partition shares at least `d`
characters, so a merged group can cross a split only through the partition's first range, and
only if that range merged the whole partition. `d958a7d` fell back to the serial prune whenever
that happened, and on chr18 at depths 3 and 4 it always did.

That rule is necessary but not sufficient. `extendRange()` carries the range into the next
partition only if (a) the LCP across the right split exceeds the LCP across the left split, which
is the range's left LCP, and (b) every partition below that right LCP is also one merged range
from the same start node and logical file; otherwise it stops at the first group that fails.
`4bf5201` tests (a) and (b) exactly. On chr18 at depth 3, eight of 143 partitions were one merged
range; seven failed (a) or (b), and one crossed, inside the end-marker root. The keys beginning
`#A`, `#C`, `#G` and `#T` (2,961 on chr18) all start at the end node's last offset, so the serial
merge collapses each of those subtrees into one group and a depth-3 split inside one cuts it; the
`##` keys start at other offsets, so the merge stops at each `#X` boundary.

`0d6443e` replaces that fallback with stitching: the crossed span (the open partition through the
last partition below the right LCP) is pruned again on one thread as one partition, and the check
repeats until nothing crosses. On chr18 at depth 3 that is four stitches, one per `#X` subtree. A node whose label interval
reaches past a split still sends the step to the serial prune; that trigger is possible from step
2 on, when merged nodes carry label intervals, and is untested at depth greater than 1. Equal-label
records from different shards may leave the parallel prune in a different order than the serial
pass; the final index's insensitivity to that order is shown at depth 1 (byte-identical joint
index) and at depth 3: with `2144b0e`, which also stitches spans that a merged node's label interval
crosses (the step 2 to 4 trigger), full chr18 builds at depth 1 and 3 are byte-identical (`.gcsa`
MD5 `02b0b8c0b3f55ab875421a22abdb047d`) with no serial fallback at any step. The depth-3 build was
slower, 0:54:09 against 0:43:26: its step 2 to 4 prunes took 332.8 to 469.3 s on 9 to 11 workers
against 74.0 to 91.4 s on 7. The extend after each prune re-partitions the data (the same join
partition counts at both depths), so the cause is inside the prune; per-partition merger setup over
every input shard, with framed-shard admission held to the 128-file ceiling, is the unmeasured
candidate. Depth 1 remains the default. In the unit fixture a stitched span's records match the serial
pass as a multiset per label but not in order: four unsorted records with one label from two start
nodes leave in a different order, because a stitched span starts a fresh merge heap. The test
compares per-label multisets for that reason. Whether such reordering occurred in the chr18 builds
was not checked, so their identical indexes do not by themselves establish that the index is
insensitive to it.

**Partition setup (2026-10-01, branch `prune-chunk-setup`).** The depth-3 slowdown was merger
setup: each partition binary-searched every input shard for its bounds, and on framed shards each
probe in a new 16 MiB block decoded the whole block. At step 2 on chr18 (143 partitions, 33
shards) that was 2,613 of 2,742 summed worker seconds. `f678b6c` locates all bounds once per shard
with a galloping search and passes them to the workers and stitched spans; `1dde002` gives each
worker runs of adjacent partitions sharing one input cache, since about 18 adjacent depth-3
partitions start in the same block. On the SSD the four chr18 prunes then took 66.9, 37.1, 49.6 and
58.5 s at depth 3 against 78.3, 54.8, 56.7 and 67.4 s at depth 1 (one run each, identical pruned
counts at every step). Framed-shard workers stay limited to 9 to 11 by memory: each keeps one
decoded block pair per shard (about 2.2 GB for 33 shards) within the prune's input-cache budget.

### Repeated reads: single-pass distribution and compressed merge output (2026-10-01, branch `prune-chunk-setup`)

**Distribution read the input once per key range.** `scanJoinRange` had every key-range worker read
and decompress every shard of its logical input and keep only its own range's records, so with 64
ranges the input was read 64 times per doubling step: on chr18 the distribution blocks issued
324-577 GiB of read requests per step for 4.9-8.9 GB of pruned paths (77% of the build's 2,183 GiB
of read requests). Raising the worker cap (`88ee230`: range workers against the concurrent
open-file budget, 64 at once instead of 14) did not help, because the repeated reads, not the cap,
were the limit (running threads 13.1 against 14.0, 29.6 more waiting on page I/O). `00c81a7` reads
each shard once: reader threads take whole shards and route each record to its range's sorter; all
ranges' sorters share the 64 GiB distribution budget, and records keep their sequential-scan
ordinals, on which `joinRecordLess` ends, so runs do not depend on reader scheduling. That alone was
slower, because every reader appends to every range and a sorter sorted and wrote a full run while
holding its lock (1-3 threads running, 70-87 waiting); `c02a9d8` sorts and writes runs outside the
lock. The `88ee230` admission is still what sizes the readers and the merge fan-in.

**Merged-graph output, compressed (opt-in).** `f9d4679` lets the serial merge write paths, ranks and
start nodes framed, each stream compressed on its own thread (`CompressedRecordStream`); the LCP stays
raw; `ReadBuffer` reads framed files. `4842598` fixed `MergedGraphReader::seek`, which tested the
descriptor to decide whether a reader was open, so framed readers never seeked and construction RSS
reached 110.9 GiB on chr18 (`ReadBuffer::isOpen()`), and a double count of framed reads in the
logged scan volume. `e40e008` makes it opt-in, `GCSA_IO_COMPRESS_MERGE=1`.

**Measured on chr18 (depth 3, sampled, nice 0, index MD5 identical in every run).** Against
`1dde002`: distribution plus join planning 0:02:59 against 0:05:12 over four steps (13.5-18.2 mean
cores against 8.1-10.5); read requests for the whole build 482.0 against 2,182.8 GiB; distribution
writes up about 3.5 GiB per step (about 0.5 GiB per sorter) and peak RSS 31.2 against 27.1 GiB. With
compressed merge output, the merge wrote 4.1 against 25.8 GiB in the same time and the final-event
scan's device reads fell from 72.7 to 15.3 GiB, but the scan took 0:08:51 against 0:07:24: its main
thread ran 75% against 71% of a longer scan and waited no more, so it does the decompression and
block checksums itself whenever it outruns the one-buffer read-ahead. Decoding through
`CompressedBlockPrefetchPool` is the untested remedy. Whole-build times (0:36:57 against 0:40:20)
are within this host's spread: preprocessing, which no change touches, differed by 0:02:16
between the two runs.

At the default depth 1 each prune writes seven shards, which capped whole-shard readers at seven
threads (distribution 0:06:09 against 0:05:02 on chr18); `5b4d225` gives readers chunks of shards
(`PathShardReader::restrict`), and chr18 distribution then took 0:02:45. A complete joint
chr2+chr18 build with `5b4d225` (default depth, compressed merge off) matched the 9/30 index MD5 and
its blocks summed to 3:45:44 at 6.08 mean threads: joins 1:16:35 (33.9%, about 9 of 32 workers),
construction 0:40:54, merge 0:29:52, preprocessing 0:22:05, prune 0:19:05 (0:27:38 on 9/30),
distribution 0:17:19 (0:37:12 on 9/30), compaction 0:16:54. The merge's drop from 1:08:05 is host
write pressure (dirty-page throttle in 36% of its samples on 9/30, 0% now), not code. Record:
`hprc_v2_vg_rna/gcsa2_distribution_concurrency_20260928/joint_chr2_chr18_5b4d225_20261002/`.

## Build, test, and usage

From the containing `vg` checkout, use its local toolchain wrapper. The GCSA2
library tests deliberately use tiny byte budgets so the same records create
multiple runs and resume checkpoints:

```bash
JOBS=1 ./build-local.sh -C deps/gcsa2 test
JOBS=8 ./build-local.sh -C deps/gcsa2 all
```

The containing `vg` Makefile must track `deps/gcsa2/src/*.cpp` as archive
prerequisites. Otherwise a source-only GCSA2 change can leave a successfully
relinked `vg` using an older copied `lib/libgcsa2.a`; the local fork carries
that dependency correction.

The `vg` integration is exercised by the ordinary index regressions and a
dedicated two-logical-input equivalence/recovery fixture:

```bash
./build-local.sh bin/vg
(cd test && prove -v t/58_vg_gcsa_external.t)
(cd test && prove -v t/06_vg_index.t t/52_vg_autoindex.t)
```

Standalone GCSA2 accepts binary de Bruijn graph inputs. The external route is
selected by a workspace; operational settings can change on resume, while a
semantic mismatch is rejected:

```bash
deps/gcsa2/bin/build_gcsa \
  --work-dir /large-local-disk/chr6.gcsa-work \
  --resume --keep-work \
  --memory-limit 96G --disk-limit 40T \
  --io-buffer-size 64M --sort-run-size 4G \
  --join-partition-size 16G --process-workers 4 --merge-fan-in 64 \
  --temp-compression auto --compression-block-size 16M \
  --compression-workers 4 --compression-level 1 --clean-obsolete \
  --max-open-files 128 --allow-path-explosion \
  -T 32 -o chr6 chr6
```

`vg index` persists the graph-derived k-mer stream before entering GCSA2. One
named source graph remains one semantic GCSA2 input even when later phases spill
to many physical files:

```bash
vg index -p -V -g chr6.gcsa -k 16 -X 4 -t 32 \
  --gcsa-work-dir /large-local-disk/chr6.gcsa-work \
  --gcsa-resume \
  --gcsa-memory-limit 96G --gcsa-disk-limit 40T \
  --gcsa-process-workers 4 \
  --gcsa-temp-compression auto --gcsa-compression-block-size 16M \
  --gcsa-compression-workers 4 --gcsa-compression-level 1 \
  --gcsa-clean-obsolete \
  --gcsa-sort-run-size 64G --gcsa-join-partition-size 64G \
  -f chr6.mapping chr6.pruned.vg
```

For a hard end-to-end process ceiling, including allocator overhead and charged
filesystem cache, use the benchmark harness. It binds input and executable
hashes before launch, records the exact command and `/usr/bin/time -v`, and
samples phase, anonymous/cache/dirty memory, cgroup CPU, block I/O, tasks, and
live workspace bytes. `phase_summary.tsv` aggregates those counters by phase;
`summary.tsv` records the whole-run maxima and totals:

```bash
TMPDIR=/large-local-disk/chr6-run/vg-tmp \
scripts/benchmark-gcsa-external.sh \
  --vg ./bin/vg --graph chr6.pruned.vg --mapping chr6.mapping \
  --run-dir /large-local-disk/chr6-run \
  --memory-limit 94G --cgroup-limit 96G --disk-limit 40T \
  --sort-run-size 64G --join-partition-size 64G \
  --process-workers 4 --threads 32 \
  --temp-compression auto --compression-block-size 16M \
  --compression-workers 4 --compression-level 1 --clean-obsolete \
  --kmer-length 16 --doubling-steps 4
```
