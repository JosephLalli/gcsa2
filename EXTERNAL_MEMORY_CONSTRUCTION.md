# External-memory construction foundation

The external route is a fresh, serial alternative to the resident GCSA2
builder. `ConstructionParameters::setWorkDirectory()` selects it, and
`GCSA::buildAndStore()` writes an ordinary GCSA followed by an ordinary LCP.
The resident constructor and both public file formats remain unchanged.

This foundation has no workspace, resume, checkpoint, worker command,
compression, telemetry, or pair-atomic publication protocol. The memory value
is a construction-buffer target. It excludes the input graph, final succinct
supports, allocator and library overhead, thread runtime, and page cache, so it
is not a process RSS guarantee.

## Stage and ownership map

| Stage | Owning modules and consumer | Format and ownership invariant | Failure and focused check |
| --- | --- | --- | --- |
| Configuration and entry | `support.*`, `utils.*`, `build_gcsa.cpp`, `gcsa.h`, `gcsa.cpp`; consumed by `build_gcsa` and vg | `setWorkDirectory()` stores the selected path privately. `GCSA::buildAndStore()` applies it to the process-global `TempFile` owner when construction starts. | Invalid byte strings, setter-order clamping, and the configuration setter's lack of an early global side effect are covered by `test_parameters`. A nonexternal call to `buildAndStore()` throws. |
| Bounded input decoding | `files.*`, `external_preprocessing.*`; consumed by the external GCSA constructor | `InputGraph::scanKMerBlocks()` preserves one logical input identity across callback blocks. Preprocessing owns sorted key/start streams and deletes them on destruction. | Record counts, text splitting, key merging, start nodes, and physical-to-logical shard mapping are covered by `test_external_preprocessing`. |
| Initial supports | `dbg.*`, the `LCP` additions in `path_graph.*`; consumed by `ExternalInputPreprocessor` | File constructors scan the sorted unique-key stream through bounded buffers and construct the same resident query supports as the vector constructors. | The preprocessing test compares counts, labels, and starts; whole-construction gates compare final bytes. |
| Fixed-record sorting | `external_sort.*`; consumed by preprocessing, spillable start-node sets, and redundancy sorting | Raw fixed-width runs are private scratch. The sorter owns registered staging files, uses bounded run/merge buffers, and installs one output only after success. | `test_external_sort` forces multiple runs and merge passes, checks budget statistics, reduction, malformed input, and descriptor limits. |
| Prefix doubling | `external_join.cpp`, `path_graph_external.h`, and external portions of `path_graph.cpp`; consumed by the GCSA doubling loop | A `logical_file_id_t` carries graph semantics; a `physical_shard_id_t` only identifies a spill shard. Paired path/rank streams remain aligned through sort-merge joins. | `test_external_join` compares the resident result, exercises same-logical and separate-logical shards, forced runs, blocked key groups, and the byte budget. `test_external_path_sort` covers paired-stream corruption and descriptor admission. |
| Prune and final merge | Spillable range/deque state and bounded input descriptors in `path_graph.cpp`; consumed by `PathGraph::prune()` and `MergedGraph` | Arbitrarily long equal-label/from sets may spill. The current generation remains owned until its successor is complete. Stored-byte accounting saturates rather than wrapping. | `test_path_graph_prune` compares resident and spilled results, including a late spill with nonzero absolute offsets and many physical shards. |
| Final ordered scan | `resources.*`, `disk_array.*`, `final_events.*`, and the external branch in `gcsa.cpp` | One ordered scan writes ephemeral little-endian event streams. `FinalEventFiles` owns and removes them. `MemoryBudget` reservations cover the disk arrays, event buffers, and spillable node sets used by this stage. | `test_disk_array`, `test_internal_buffers`, and `test_final_events` cover eviction, reservation release, forced node-set spill, malformed streams, and exact SDSL component bytes. |
| GCSA publication | `storeFinalComponents()` in `final_events.cpp`; consumed only by `GCSA::buildAndStore()` and its exact-format test | Components are encoded directly in the existing `GCSA::load()` order. One partial output is synced and renamed; private event streams are validated by size and semantic metadata held by the same construction. | `test_final_events` reloads the published index and compares its complete serialization with the resident equivalent. Output-failure gates require no published replacement and empty scratch. |
| LCP publication | `lcp.cpp`, with the private friend entry in `lcp.h`; consumed by `GCSA::buildAndStore()` | Raw leaf generations are temporary; direct serialization preserves the existing LCP format. GCSA publishes first, then LCP. A failed LCP can leave a complete GCSA and the input leaf available for retry. | The serial fixture gate compares resident/external LCP bytes and reloads them. Publication-failure tests must check the ordered partial-success contract explicitly. |
| Build and installed surface | `Makefile` and `install.sh` | The archive contains the external modules. Only the established nine GCSA headers are installed; `disk_array.h`, `external_preprocessing.h`, `external_sort.h`, `final_events.h`, `path_graph_external.h`, and `resources.h` remain source-private. | Compile and run an installed-header consumer against the built archive, then build the actual vg caller. |

## Installed interface additions

The installed additions each have a production consumer:

- `GCSA::buildAndStore()` is the paired external entry used by the native
  builder and vg.
- Byte-valued disk/memory setters, `setWorkDirectory()`, `externalMemory()`,
  and `parseBytes()` configure that entry without changing resident defaults.
- `InputGraph::scanKMerBlocks()` is the bounded decode seam used by external
  preprocessing and by the optional bounded verifier.
- The file-backed `DeBruijnGraph` and `LCP` constructors build existing support
  types from preprocessing output.
- Logical/physical shard IDs, stored-byte accounting, and the bounded
  prune/merge arguments preserve semantics when one logical input spills into
  several files.
- Byte-sized `ReadBuffer`/`WriteBuffer` configuration and cache retirement keep
  the existing path readers usable by bounded merge stages.

Formatting helpers, construction counters, an unread event-metadata file, and
test-only forwarding parameters are intentionally absent. The external helper
headers are compiled with the library but are not copied by `install.sh`.

## Deliberate specializations

The generic fixed-record sorter is reused for keys, nodes, redundancy, and
spillable node sets. The path sorter remains separate because one logical
record spans a `PathNode` stream and a variable-length rank sidecar; its atomic
install, descriptor admission, and corruption checks apply to that pair.

The external final encoders reproduce SDSL byte layouts directly because the
ordinary SDSL builders materialize path-proportional vectors. Their exact-byte
tests invoke the production encoders. Consolidating them into a general codec
would add an interface without another consumer.

The legacy asynchronous `ReadBuffer` remains in use for path streams. Small raw
event readers use checked exact I/O and exceptions so component publication can
fail before rename. Replacing either with the other would change failure and
seek behavior beyond this foundation.

## Publication and failure contract

Scratch names belong to `TempFile` or to a local RAII owner. Registered staging
files are removed on exceptions and on the inherited direct-exit paths. Output
files use same-directory partial files, sync, rename, and directory sync.

GCSA and LCP are individually atomic and ordered, not pair-atomic. An early
failure must preserve any prior destination. If GCSA publication succeeds and
LCP publication fails, the new GCSA remains. Retry must use the retained LCP
leaf; no resume manifest is implied.

## Required validation

Use the normal vg/owner build instructions in a fresh checkout and put
`TMPDIR` on SSD storage. After vg has built its dependencies, the owner checks
are available as `make -C deps/gcsa2 test` and `make -C deps/gcsa2 all`; rebuild
vg after installing matching GCSA headers and its archive. The local performance
fork uses `build-local.sh` to invoke these targets safely around frozen binaries;
that machine-specific wrapper is not required by this library contribution.

The bounded acceptance fixtures below do not require chromosome datasets.

At one thread, build empty, zero-doubling, and cyclic fixtures through both
routes; require byte-identical GCSA/LCP files, successful reload/query checks,
and empty scratch. Also inject too-small memory, unwritable output, and paired
path/rank install failures; require the specified failure, unchanged prior
outputs, unchanged source path/rank files, and empty scratch. The installed
consumer must compile using only the headers copied by `install.sh`.

Optional commits must be replayed onto this foundation and rerun their own
focused tests. A clean replay is compatibility evidence; it does not approve
the optional behavior.

## Owner decisions still required

- Is ordered GCSA-then-LCP publication acceptable to callers, including the
  observable complete-GCSA/failed-LCP state?
- Is the process-global `TempFile` directory acceptable for the library entry,
  given that simultaneous constructions in one process cannot select different
  scratch directories?
- Should the file-backed `DeBruijnGraph` and `LCP` constructors eventually be
  hidden behind a friend or factory, despite the extra header churn?
- Does the serial foundation pass the supported owner platforms and the actual
  vg integration gate? Historical fixture receipts do not answer that question
  for a new tree.
