/*
  Copyright (c) 2026 GCSA2 contributors

  Bounded fixed-record external sorting.
*/

#ifndef GCSA_EXTERNAL_SORT_H
#define GCSA_EXTERNAL_SORT_H

#include <gcsa/utils.h>

#include <functional>
#include <iosfwd>
#include <string>

namespace gcsa
{

/*
  This is deliberately a byte-record API instead of a template. Construction
  stages write compact binary records already, and a non-template interface
  keeps one tested run/merge implementation available to key, node, and future
  event streams.

  Comparator returns a negative, zero, or positive value. If ordering among
  records that compare equal matters, callers must provide a total order. The
  sorter otherwise uses deterministic run and input positions only as an
  implementation tie-breaker; reducers must treat equal-key groups as sets.
*/
struct ExternalFixedRecordSortStats
{
  size_type runs, merge_operations, merge_passes;
  size_type max_records_resident, max_bytes_resident;

  ExternalFixedRecordSortStats() :
    runs(0), merge_operations(0), merge_passes(0),
    max_records_resident(0), max_bytes_resident(0)
  {
  }
};

class ExternalFixedRecordSorter
{
public:
  typedef std::function<int(const void*, const void*)> Comparator;

  /*
    The reducer is called once for each sorted input record. first_in_group and
    last_in_group delimit comparator-equal records. It may keep bounded state
    for the current group and write an arbitrary result to output. This makes
    unique/merged streams possible without a second global vector.
  */
  typedef std::function<void(const void* record, bool first_in_group,
    bool last_in_group, std::ostream& output)> Reducer;

  static size_type minimumBudget(size_type record_bytes);

  /*
    total_order asserts that compare returns zero only for byte-identical
    records. Run formation may then sort the records in place rather than a
    permutation of their offsets. It is safe exactly when the comparator
    inspects the whole record: no equal-but-distinct records exist, so the
    input-order tie-break it drops cannot be observed in the output. Leave it
    false for a comparator that orders on a key only.

    Two things improve. Every comparison loses an indirect call to two random
    offsets. And the offset array, which is as large as the data for an
    eight-byte record, is not allocated: an input that fits in one run holds
    only its own bytes, while an input that does not gets more records per run
    out of the same budget, so fewer runs and fewer merge passes.
  */

  // Sort fixed-width binary records from input_name into output_name.
  static void sort(const std::string& input_name, const std::string& output_name,
    size_type record_bytes, size_type byte_budget, size_type requested_fan_in,
    const Comparator& compare, ExternalFixedRecordSortStats* stats = nullptr,
    bool total_order = false);

  // Sort and consume equal-key groups while writing output_name. The reducer,
  // not this class, defines the resulting output record format.
  static void sortAndReduce(const std::string& input_name,
    const std::string& output_name, size_type record_bytes,
    size_type byte_budget, size_type requested_fan_in,
    const Comparator& compare, const Reducer& reducer,
    ExternalFixedRecordSortStats* stats = nullptr, bool total_order = false);
};

} // namespace gcsa

#endif // GCSA_EXTERNAL_SORT_H
