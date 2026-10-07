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
  // ASCENDING_U64 is an eight-byte little-endian unsigned order.
  enum class RecordOrder { COMPARATOR, ASCENDING_U64 };

  /*
    The reducer is called once for each sorted input record. first_in_group and
    last_in_group delimit comparator-equal records. It may keep bounded state
    for the current group and write an arbitrary result to output. This makes
    unique/merged streams possible without a second global vector.
  */
  typedef std::function<void(const void* record, bool first_in_group,
    bool last_in_group, std::ostream& output)> Reducer;

  static size_type minimumBudget(size_type record_bytes);

  // Callers may set total_order only when comparator equality implies
  // byte-identical records. Run formation can then reorder the records
  // directly instead of preserving equal-record input order through offsets.
  // Sort fixed-width binary records from input_name into output_name.
  static void sort(const std::string& input_name, const std::string& output_name,
    size_type record_bytes, size_type byte_budget, size_type requested_fan_in,
    const Comparator& compare, ExternalFixedRecordSortStats* stats = nullptr,
    bool total_order = false, RecordOrder order = RecordOrder::COMPARATOR);

  // Sort and consume equal-key groups while writing output_name. The reducer,
  // not this class, defines the resulting output record format.
  static void sortAndReduce(const std::string& input_name,
    const std::string& output_name, size_type record_bytes,
    size_type byte_budget, size_type requested_fan_in,
    const Comparator& compare, const Reducer& reducer,
    ExternalFixedRecordSortStats* stats = nullptr, bool total_order = false,
    RecordOrder order = RecordOrder::COMPARATOR);
};

} // namespace gcsa

#endif // GCSA_EXTERNAL_SORT_H
