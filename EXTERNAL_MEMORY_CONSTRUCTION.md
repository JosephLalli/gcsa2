# External-memory GCSA2 construction

Status: implementation in progress on top of GCSA2 `2a1d479`. The public
`.gcsa` and `.lcp` formats remain unchanged. This document records both the
measured baseline and the invariants of the disk-first route; the final section
tracks which parts are actually implemented.

## Baseline audit

The existing builder uses temporary files, but a temporary file is generally a
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

`DiskBudget` tracks committed live workspace bytes, bytes reserved by active
writers, configured maximum bytes, `statvfs()` free bytes, and a free-space
safety margin. A generation estimate is diagnostic. Construction stops only
when an actual reservation cannot be honored, a write fails, an explicit growth
cutoff is reached, or input is invalid. RAM size is never a disk-limit proxy.

## Durable workspace

The workspace is persistent and semantic:

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

`build.json` contains the workspace format, GCSA2 version, final format
versions, input identities/sizes/checksums, mapping identity, k-mer length,
doubling steps, sample period, LCP branching, alphabet, and logical input IDs.
Those fields form the semantic fingerprint. Memory, threads, run size,
partition target, merge fan-in, open-file limit, and I/O buffer size are stored
for provenance but excluded from the fingerprint.

Every binary artifact has an explicit little-endian header and footer with
magic, format version, artifact kind, logical ID, physical shard ID, record and
payload byte counts, checksum, sort order, and key range. On-disk phase records
are explicit formats; raw C++ layouts are used only for the unchanged legacy
final index format.

### Commit and recovery protocol

1. Write and checksum an artifact incrementally at a unique `.partial` path.
2. Rewrite/finalize its header and append its footer.
3. Flush, `fdatasync()`, close, rename to the semantic final name, and sync its
   parent directory.
4. After every output artifact for a task exists, validate their metadata and
   atomically commit a task record under `tasks/` using the same sync protocol.
5. Only a task record makes its outputs part of the restart frontier.

On resume, partial files and artifacts without a committed task are removed or
quarantined; immutable inputs and committed predecessors remain. Task records
are replayed in dependency order. Headers, lengths, footers, and task metadata
are always checked. Payload checksums are selected by validation level. A task
with a missing, truncated, or corrupt output is invalidated and rerun along with
its dependants. Cleanup is itself idempotent and never removes the newest
committed predecessor of an incomplete task.

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

A `PathGraph` is a manifest of immutable shards, each carrying logical ID,
physical ID, sort order, counts, byte size, key range, and checksum. The next
phase consumes a run set directly. Compaction limits fan-in but does not force a
single logical input into one giant file.

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

A sampling pass records key/range histograms, heavy keys, leading label ranks,
and fanout estimates. It commits a deterministic range partition plan. MSD
radix distribution writes left and right pack files with per-bin offsets, so the
number of logical bins is independent of open file count.

For each partition, the scheduler reserves a declared working set. If one side
fits, it is indexed and the other is streamed. Otherwise the partition is
recursively split on more key bits. A single unsplittable heavy key is processed
as bounded left blocks crossed with bounded right blocks; generated records are
immediately sent to label-run construction. This uses
`O(|L_key| * |R_key|)` generation work, as required by the output, but only
`O(block_L + block_R + output_buffer)` RAM.

Sorted paths bypass the join and enter label-run generation unchanged. A later
optimization may encode one left reference followed by many right references,
delaying repeated-label materialization.

The first production slice uses an external sort-merge implementation of the
same plan. For every logical input it writes explicit, versioned left-by-`to`
and right-by-`from` run records. Run creation uses leveled compaction, so file
name metadata is `O(F log_F R)` for merge fan-in `F` and `R` initial runs. A
matching key is processed without materializing either key group: the right
range is replayed for each left record. This is deliberately slower than the
planned sampled range partitions, but it already handles an individual key
larger than RAM and preserves joins across physical shards. Generated records
flow directly to the bounded external label sorter.

### Pruning

`SpillableGroup` keeps only a bounded prefix plus a summary: first/last label,
common LCP, first `from`, all-from-equal, first logical ID,
all-logical-equal, predecessor union, count, and required minima/maxima. Records
beyond the cap spill sequentially. Summary-only decisions avoid rereads; output
branches replay the group once. Physical shard IDs never participate in
`SameFromFile` semantics.

### Final components and LCP

The final merged-graph scan emits immutable BWT-set, edge-degree, sample,
occurrence, redundancy, and LCP-leaf event streams. Unordered events are
externally sorted and aggregated. Dense previous-occurrence state is a
disk-backed block array: each scan block commits output events and an idempotent
sorted assignment log before replaying the log into the synced base array.

Each final component is built and released separately. A packer writes the
unchanged `GCSA::serialize()` component order to a partial `.gcsa`, validates
it, and atomically publishes it. LCP level `i+1` is produced by streaming groups
of the configured branching factor from level `i`; levels are checkpointed and
packed in legacy `LCPArray::serialize()` order without simultaneous residency.

## Implementation ledger

The ledger is intentionally conservative: a feature moves to **complete** only
when its production call path and forced-spill/recovery tests pass.

| Slice | State |
| --- | --- |
| Shared byte-token `MemoryBudget` and statvfs-aware `DiskBudget` | in review |
| Versioned workspace artifacts, atomic publication, manifest compatibility, recovery tests | in review |
| Human-readable operational construction parameters | in review |
| Bounded external path-label runs and leveled multi-pass merge | implemented; forced-spill tested |
| Logical/physical `PathGraph` identity in pruning and joining | implemented; durable run-set manifest pending |
| External prefix-doubling join | bounded sort-merge implemented; sampled range partitioning pending |
| External keys/start nodes/initial paths | not implemented |
| Spillable pruning groups | not implemented |
| Final event/component passes and transactional `prev_occ` | not implemented |
| Streaming LCP levels and legacy packer | not implemented |
| Standalone and vg CLI integration | not implemented |
| Chromosome-scale benchmark | not run |

Current limitations are intentionally explicit. The external route is selected
only when `ConstructionParameters::work_directory` is nonempty. Key extraction,
start-node extraction, initial path construction, pruning groups, final GCSA
components, and LCP levels still use the legacy resident algorithms. Join runs
are checksummed immutable files but are not yet registered as resumable
workspace tasks. The current heavy-key fallback uses one-record left blocks and
repeated right-range reads; it is bounded but can perform much more I/O than a
larger blocked implementation. These limitations mean the full construction
does not yet satisfy the end-to-end RAM invariant, even though prefix-doubling
extension itself no longer requires a chromosome or join key to fit in memory.
