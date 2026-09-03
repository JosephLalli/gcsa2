# External-memory GCSA2 construction

Status: implementation in progress on top of GCSA2 `2a1d479`. The public
`.gcsa` and `.lcp` formats remain unchanged. This document records both the
measured baseline and the invariants of the disk-first route; the final section
tracks which parts are actually implemented.

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
preprocessing, prefix-doubling joins and label sorting, final raw component
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
and 105,438,580 KiB maximum RSS (100.55 GiB). This is the primary target for
the requested approximately 25 GiB comparison once the host's other heavy
`vg prune`/GCSA work and I/O contention are absent.

Only preexisting pruned fixtures are eligible for acceptance runs. The chr19
production prune exited successfully after 49:45:06 at 154,737,660 KiB maximum
RSS, but its wrapper did not publish the 748,778,080-byte `.partial` graph as a
committed fixture. No competing chr19 index build was launched. The remaining
final-component memory gaps below must be closed or bounded before that race is
considered ready.

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
replay remains bounded. `posix_fadvise()` is best-effort and non-semantic;
`fdatasync()` is the durability boundary. Pruning and `MergedGraph` now use
bounded spill structures and byte-sized buffers, but their reservations, final
SDSL component allocations, and some library-owned compression scratch are not
all admitted through the process-wide `MemoryBudget`. The in-process token
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

The largest remaining I/O opportunity is further record compaction, not extra
writer processes. Current fixed-width join runs repeat complete path labels and
may be larger than the variable-width source path/rank pair. A versioned grouped
codec now writes one left context for consecutive expansions in each label-sort
run, followed by compact right-dependent records that reference it. This removes
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
`--disk-limit` and never compares those bytes with RAM. The `DiskBudget`
primitive can also inspect recursively live workspace bytes and `statvfs()`
free space, but it is not yet shared by every writer. In particular, final-event
copies, LCP raw levels, and some checkpoint copies currently stop on a failed
write/ENOSPC rather than reserving their complete output first. Consequently
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
for provenance but excluded from the fingerprint.

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
newly synced run was written and avoids immediately rereading it. Phase
path/rank checkpoints use immutable raw fixed-record payloads: the checkpoint
hard-links them into the workspace when source and workspace share a filesystem
and otherwise performs one bounded copy. The payload is checksummed once while
it is committed. Normal same-filesystem restore validates task identity, record
count, and byte length and then hard-links the committed inode into the next
`PathGraph`; `--verify-workspace` additionally rereads the complete checksum.
Join-range worker outputs and final-event streams use the same immutable raw
payload protocol, avoiding both the checkpoint rewrite and normal resume copy.
Final events retain a checked header/footer compatibility restore for older
workspaces. Cross-filesystem restore copies and checksum-validates in one pass.
A task with a missing, truncated, or corrupt output is invalidated and rerun
along with its dependants. Cleanup is itself idempotent and never removes the
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
partition outputs. Completed partition task names carry the same run checksums,
so equal-length replacement runs cannot reuse stale output. The already-sorted
left/right join files act as shared pack files:
thousands of logical radix bins are half-open offsets in two files instead of
thousands of open descriptors. An exact streaming group pass remains
authoritative for every range, count, and byte estimate, so a sampling miss can
change task balance but cannot change join semantics.

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
for occurrence, continuation, and sample passes. All streams are
synced and committed together as one workspace task; a crash restarts only this
ordered scan. Mid-scan resume and the assignment-log protocol remain future
work.

Final raw components are built serially, one at a time, from those streams.
The library retains the resident-object constructor for compatibility, while
`GCSA::buildAndStore()` writes each fast/sparse BWT, rank support, edge, sample,
occurrence, and redundancy component immediately in the exact existing
`GCSA::load()` order and releases it before constructing the next component.
Fast BWT, edge-boundary, and sampled-path `bit_vector_il<512>` members are
serialized directly from their event streams, including SDSL's interleaved
cumulative counts and bounded select samples. They never materialize the dense
source bitvector or the full interleaved destination; empty `rank_support_il`
payloads are emitted without rebuilding the corresponding vector.
The partial file is synced and atomically published. Standalone `build_gcsa`,
`vg index`, and `vg autoindex` use this direct packer whenever a workspace is
configured; they reload the index only for explicit verification. This removes
the cumulative completed-index peak, but sparse SDSL builders, packed sample
IDs and boundaries, occurrence/redundancy members, and their conversion scratch
can still allocate outside `MemoryBudget`.

LCP level `i+1` is produced by streaming groups of the configured branching
factor from level `i`; each raw internal level is independently checkpointed
and consumed. The final packed hierarchy necessarily resides in
`LCPArray::data`, in legacy leaf-to-root serialization order.

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

| Slice | State |
| --- | --- |
| Shared byte-token `MemoryBudget` and statvfs-aware `DiskBudget` | memory primitives and major external-phase byte caps implemented; pruning/merged-graph buffers are byte-bounded but not globally admitted; SDSL/library allocations remain outside the shared budget; disk guard covers generated path/join volume but not every final/LCP/checkpoint writer |
| Versioned workspace artifacts, atomic publication, manifest compatibility, recovery tests | phase and join-range checkpoints plus abrupt-exit commit-boundary tests implemented; transient distribution-run task records pending |
| Human-readable operational construction parameters | core budget/run controls implemented and parser-tested; cleanup/cadence controls are not all wired |
| Bounded external path-label runs and leveled multi-pass merge | versioned grouped/prefix-compressed runs, shared-left-context references, parallel in-place run sorting, rolling cache windows, one-run bypass, direct leveled merge of already-sorted worker shards, and forced multi-pass spilling implemented and tested |
| Logical/physical `PathGraph` identity in pruning and joining | implemented; pruning coalesces only by logical ID and never treats a physical shard as a semantic input; durable run-set manifest pending |
| External prefix-doubling join | bounded sort-merge, rolling page-cache windows, RAM-capped persisted 4096-record MSD radix range-pack planning bound to run checksums, exact range contracts, recursive two-dimensional heavy-key splitting, blocked nested-loop expansion, and selectable checksum scans implemented and forced-spill tested |
| Direct join-to-label pipeline | implemented and forced-spill tested; the unsorted generated path/rank pair is no longer materialized |
| Process worker scheduler and partition resume | implemented with fork-free spawn, global byte admission, immutable-payload semantic range checkpoints, same-filesystem zero-copy restore, and one-/multi-partition tests |
| External keys/start nodes/initial paths | implemented with bounded fixed-record runs, global duplicate reduction, streaming support construction, durable key/start checkpoints, and physical initial shards preserving logical IDs |
| Spillable pruning groups | implemented and forced-spill tested for equal-label priority groups, extended ranges, external same-from sets, and bounded input/output descriptor caches |
| Final event/component passes | implemented with one-task immutable-payload event checkpoint, same-filesystem zero-copy restore, spillable mapped start-node sets, disk-backed `prev_occ` and suffix-tree stack, external redundancy sort, streaming fast-vector serialization, and direct component-at-a-time packing in standalone/vg/autoindex; mid-scan assignment logs and remaining SDSL token admission pending |
| Streaming LCP levels | implemented and resume-tested with one raw level resident at a time and byte-identical legacy serialization; final packed hierarchy remains resident |
| External verification | implemented and forced-spill tested with bounded input blocks, external expected/actual occurrence sorts, callback-based locate, and sequential set comparison; resumable runs and parallel label ranges pending |
| Standalone and `vg index` / `vg autoindex` CLI integration | implemented; forced-spill/resume integration tested |
| Chromosome-scale benchmark | preexisting chr20 k32-pruned fixture completed and verified at 32.36 GiB RSS under a 128 GiB cgroup; outputs are byte-identical to legacy; separate 25 GiB forced-spill/recovery path exercised |

Current limitations are intentionally explicit. The external route is selected
only when `ConstructionParameters::work_directory` is nonempty. The production
file-building routes no longer retain all completed GCSA components or retain a
completed GCSA while constructing LCP. Fast BWT/edge/sample vectors stream in
bounded memory, but individual sparse SDSL builders, packed sample data,
occurrence/redundancy structures, and their conversion scratch still allocate
outside the token budget, and the final packed LCP remains resident. These are
bounded by one final component rather than total intermediate path volume, but
they keep the in-process budget from being a strict whole-process RSS ceiling;
production runs retain a cgroup limit.

Join distribution runs are checksummed immutable files but are not yet
registered as resumable workspace tasks, so an incomplete generation rebuilds
them before restoring the persisted MSD plan and completed join ranges. The
exact group-summary pass still scans both shared pack files after sampling; this
extra sequential I/O makes sampled boundaries safe rather than speculative.
Initial records and join distribution records still repeat full rank payloads.
Label-sort runs group shared left contexts and prefix-compress labels, but a
deeper reference representation spanning the pre-sort join stream remains an
optional compaction rather than a feasibility dependency.
Final component construction is serial and the ordered event scan resumes only
at its task boundary. Process workers currently accelerate independent join
ranges only. Even with these limitations, preprocessing and prefix doubling no
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

The goal is operational, not semantic, and may change on resume. A larger value
admits larger initial/sort runs, join blocks, and more worker reservations,
usually reducing physical shards, merge passes, and HDD traffic. A smaller
value deliberately produces more spill runs and I/O while preserving logical
file identity and final-index semantics. Disk limit and free-space safety—not
the RAM goal—decide whether a valid but very large generation may continue.

Before a whole-pangenome claim, the remaining work is ordered by the project's
feasibility, RSS, then speed policy:

1. Bound the `vg` input producer. `SourceSinkOverlay` now discovers components
   on the fly and retains only visited state, its traversal stack, and component
   tips, including a sparse-ID fallback. `VGset::for_each()` still loads one
   whole physical graph before GCSA2's budget exists. A sharded graph reader
   must preserve one explicit logical GCSA2 input identity across all physical
   chunks.
2. Bound the remaining resident floor: `NodeMapping`, de Bruijn support, one
   SDSL final-component builder/conversion, and the final packed LCP allocate
   outside the shared token budget. External mapping lookup, direct streaming
   serialization for the largest SDSL component, and streamed final LCP packing
   are the principal remaining greater-than-RAM changes. Until then, measured
   component floors must be reserved below the deployment hard cap.
3. Make `DiskBudget` global. Join/path generation is reserved today, but
   checkpoint payloads, final events, verification runs, and LCP levels can
   still discover ENOSPC only when a write fails.
4. Add block-level final-event recovery. The event scan is bounded but a crash
   currently restarts the complete scan; idempotent assignment logs and output
   blocks would retain finer progress.
5. Persist distribution and group-summary products as task artifacts, and emit
   exact group summaries while creating join runs. This removes the full
   post-sampling rescan and avoids rebuilding distributions after interruption.
6. Compact the repeated rank payload in join distribution records before
   increasing worker count. Only after per-path bytes and storage headroom are
   measured should more encoders, merge groups, or verifier label ranges run in
   parallel under combined memory and I/O admission.

Acceptance should report the configured GCSA goal and cgroup cap separately,
along with peak anonymous/cache/PSS memory, live disk, records and bytes per
generated path, merge levels, worker count, filesystem throughput/PSI, and
total `(bytes read + bytes written) / source byte`. Whole-genome runtime should
not be extrapolated until the projected resident floor fits the selected cap.

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
  --process-workers 4 --threads 32 --kmer-length 16 --doubling-steps 4
```
