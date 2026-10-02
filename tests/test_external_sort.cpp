#include <gcsa/external_sort.h>

#include <algorithm>
#include <array>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

using namespace gcsa;

namespace
{

struct Record
{
  std::uint64_t key, value;
};

struct Summary
{
  std::uint64_t key, sum, count;
};

void require(bool value) { if(!value) { std::abort(); } }

int compareKey(const void* left, const void* right)
{
  const Record& a = *static_cast<const Record*>(left);
  const Record& b = *static_cast<const Record*>(right);
  return (a.key < b.key ? -1 : (a.key > b.key ? 1 : 0));
}

int compareRecord(const void* left, const void* right)
{
  int result = compareKey(left, right);
  if(result != 0) { return result; }
  const Record& a = *static_cast<const Record*>(left);
  const Record& b = *static_cast<const Record*>(right);
  return (a.value < b.value ? -1 : (a.value > b.value ? 1 : 0));
}

std::string readBytes(const std::string& name)
{
  std::ifstream input(name.c_str(), std::ios_base::binary);
  return std::string(std::istreambuf_iterator<char>(input),
    std::istreambuf_iterator<char>());
}

// A record whose comparator is a total order but whose width is not a whole
// number of 64-bit words. The sorter must ignore the declaration rather than
// mis-sort, so the flag stays an optimisation and never a correctness input.
struct OddRecord
{
  std::uint32_t key;
  std::uint8_t  tag;
};

int compareOddRecord(const void* left, const void* right)
{
  const OddRecord& a = *static_cast<const OddRecord*>(left);
  const OddRecord& b = *static_cast<const OddRecord*>(right);
  if(a.key != b.key) { return (a.key < b.key ? -1 : 1); }
  return (a.tag < b.tag ? -1 : (a.tag > b.tag ? 1 : 0));
}

std::vector<Record> readRecords(const std::string& name)
{
  std::ifstream input(name.c_str(), std::ios_base::binary);
  std::vector<Record> result;
  Record record;
  while(input.read(reinterpret_cast<char*>(&record), sizeof(record))) { result.push_back(record); }
  require(input.eof());
  return result;
}

} // namespace

int main()
{
  char root[] = "/tmp/gcsa-fixed-sort-XXXXXX";
  require(mkdtemp(root) != nullptr);
  TempFile::setDirectory(root);
  const std::string input = std::string(root) + "/input.bin";
  const std::string sorted = std::string(root) + "/sorted.bin";
  const std::string reduced = std::string(root) + "/reduced.bin";

  std::vector<Record> records;
  std::map<std::uint64_t, Summary> expected;
  for(std::uint64_t i = 0; i < 500; i++)
  {
    Record record = { (97 * i + 11) % 37, (271 * i + 5) % 1009 };
    records.push_back(record);
    Summary& summary = expected[record.key];
    summary.key = record.key; summary.sum += record.value; summary.count++;
  }
  std::ofstream output(input.c_str(), std::ios_base::binary);
  require(static_cast<bool>(output));
  output.write(reinterpret_cast<const char*>(records.data()), records.size() * sizeof(Record));
  output.close();

  const size_type budget = ExternalFixedRecordSorter::minimumBudget(sizeof(Record));
  ExternalFixedRecordSortStats sort_stats;
  ExternalFixedRecordSorter::sort(input, sorted, sizeof(Record), budget, 2,
    compareRecord, &sort_stats);
  std::vector<Record> observed = readRecords(sorted);
  require(observed.size() == records.size());
  require(std::is_sorted(observed.begin(), observed.end(), [](const Record& a, const Record& b)
  {
    return (a.key < b.key || (a.key == b.key && a.value < b.value));
  }));
  require(sort_stats.runs > 2 && sort_stats.merge_operations > 0);
  require(sort_stats.merge_passes >= 2);
  require(sort_stats.max_records_resident < records.size());
  require(sort_stats.max_bytes_resident <= budget);

  Summary current = { 0, 0, 0 };
  ExternalFixedRecordSortStats reduce_stats;
  ExternalFixedRecordSorter::sortAndReduce(input, reduced, sizeof(Record), budget, 2,
    compareKey,
    [&current](const void* value, bool first, bool last, std::ostream& stream)
    {
      const Record& record = *static_cast<const Record*>(value);
      if(first) { current = { record.key, 0, 0 }; }
      current.sum += record.value; current.count++;
      if(last) { stream.write(reinterpret_cast<const char*>(&current), sizeof(current)); }
    }, &reduce_stats);
  std::vector<Summary> summaries;
  std::ifstream summary_input(reduced.c_str(), std::ios_base::binary);
  Summary summary;
  while(summary_input.read(reinterpret_cast<char*>(&summary), sizeof(summary))) { summaries.push_back(summary); }
  require(summary_input.eof());
  require(summaries.size() == expected.size());
  for(size_type i = 0; i < summaries.size(); i++)
  {
    require(summaries[i].key == expected[summaries[i].key].key);
    require(summaries[i].sum == expected[summaries[i].key].sum);
    require(summaries[i].count == expected[summaries[i].key].count);
  }
  require(reduce_stats.runs > 2 && reduce_stats.max_bytes_resident <= budget);

  // Declaring a total order is an optimisation, never a change of result: the
  // same comparator must produce the same bytes with and without it. Run
  // formation is the only thing that differs, so a budget small enough to force
  // several runs is what actually exercises the difference.
  {
    const std::string permuted = std::string(root) + "/permuted.bin";
    const std::string in_place = std::string(root) + "/in-place.bin";
    ExternalFixedRecordSortStats permuted_stats, in_place_stats;
    ExternalFixedRecordSorter::sort(input, permuted, sizeof(Record), budget, 2,
      compareRecord, &permuted_stats, false);
    ExternalFixedRecordSorter::sort(input, in_place, sizeof(Record), budget, 2,
      compareRecord, &in_place_stats, true);
    require(readBytes(permuted) == readBytes(in_place));
    require(permuted_stats.runs > 2);
    // Sixteen-byte records, so dropping the eight-byte offset per record buys
    // half again as many records per run.
    require(in_place_stats.runs < permuted_stats.runs);
    require(in_place_stats.max_bytes_resident <= budget);

    // sortAndReduce sees the same equal-key groups either way. Reducing on the
    // key alone while sorting on the whole record is exactly the case the flag
    // is meant to survive: groups are formed by the comparator that was passed.
    const std::string permuted_reduced = std::string(root) + "/permuted-reduced.bin";
    const std::string in_place_reduced = std::string(root) + "/in-place-reduced.bin";
    auto count_groups = [](const void* value, bool first, bool, std::ostream& stream)
    {
      if(first) { stream.write(static_cast<const char*>(value), sizeof(Record)); }
    };
    ExternalFixedRecordSorter::sortAndReduce(input, permuted_reduced, sizeof(Record),
      budget, 2, compareRecord, count_groups, nullptr, false);
    ExternalFixedRecordSorter::sortAndReduce(input, in_place_reduced, sizeof(Record),
      budget, 2, compareRecord, count_groups, nullptr, true);
    require(readBytes(permuted_reduced) == readBytes(in_place_reduced));

    std::remove(permuted.c_str()); std::remove(in_place.c_str());
    std::remove(permuted_reduced.c_str()); std::remove(in_place_reduced.c_str());
  }

  // Declaring ASCENDING_U64 must not change a byte: the fast path sorts the
  // same eight-byte records as scalars instead of calling the comparator, so
  // its output has to match the comparator path exactly, and must also match a
  // caller that declares the order but whose width the fast path rejects.
  {
    const std::string u64_input = std::string(root) + "/u64.bin";
    const std::string via_cmp = std::string(root) + "/u64-cmp.bin";
    const std::string via_fast = std::string(root) + "/u64-fast.bin";
    std::vector<std::uint64_t> values;
    for(std::uint64_t i = 0; i < 5000; i++)
    {
      values.push_back((1181783497276652981ULL * (i + 1)) ^ (i << 17));
    }
    {
      std::ofstream out(u64_input.c_str(), std::ios_base::binary);
      out.write(reinterpret_cast<const char*>(values.data()),
        values.size() * sizeof(std::uint64_t));
    }
    auto compare_u64 = [](const void* left, const void* right)
    {
      std::uint64_t a, b;
      std::memcpy(&a, left, sizeof(a)); std::memcpy(&b, right, sizeof(b));
      return (a < b ? -1 : (a > b ? 1 : 0));
    };
    const size_type u64_budget =
      ExternalFixedRecordSorter::minimumBudget(sizeof(std::uint64_t)) * 4;
    ExternalFixedRecordSorter::sort(u64_input, via_cmp, sizeof(std::uint64_t),
      u64_budget, 2, compare_u64, nullptr, true,
      ExternalFixedRecordSorter::RecordOrder::COMPARATOR);
    ExternalFixedRecordSorter::sort(u64_input, via_fast, sizeof(std::uint64_t),
      u64_budget, 2, compare_u64, nullptr, true,
      ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64);
    require(readBytes(via_cmp) == readBytes(via_fast));
    std::sort(values.begin(), values.end());
    std::string sorted_bytes(reinterpret_cast<const char*>(values.data()),
      values.size() * sizeof(std::uint64_t));
    require(readBytes(via_fast) == sorted_bytes);
    std::remove(u64_input.c_str()); std::remove(via_cmp.c_str());
    std::remove(via_fast.c_str());
  }

  // A run of at least 2^20 records formed on the main thread is sorted in
  // parallel. Under a declared total order equal records are byte-identical,
  // so the parallel result must equal a serial std::sort byte for byte, with
  // many duplicates, for both the scalar path and the comparator path.
  {
    const std::string big_input = std::string(root) + "/big-u64.bin";
    const std::string big_sorted = std::string(root) + "/big-u64-sorted.bin";
    std::vector<std::uint64_t> values;
    for(std::uint64_t i = 0; i < (std::uint64_t(3) << 20); i++)
    {
      values.push_back(((1181783497276652981ULL * (i + 1)) ^ (i << 17)) % 1000003);
    }
    {
      std::ofstream out(big_input.c_str(), std::ios_base::binary);
      out.write(reinterpret_cast<const char*>(values.data()),
        values.size() * sizeof(std::uint64_t));
    }
    auto compare_u64 = [](const void* left, const void* right)
    {
      std::uint64_t a, b;
      std::memcpy(&a, left, sizeof(a)); std::memcpy(&b, right, sizeof(b));
      return (a < b ? -1 : (a > b ? 1 : 0));
    };
    ExternalFixedRecordSortStats big_stats;
    ExternalFixedRecordSorter::sort(big_input, big_sorted, sizeof(std::uint64_t),
      256 * MEGABYTE, 2, compare_u64, &big_stats, true,
      ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64);
    require(big_stats.runs == 1);
    std::sort(values.begin(), values.end());
    require(readBytes(big_sorted) == std::string(reinterpret_cast<const char*>(values.data()),
      values.size() * sizeof(std::uint64_t)));
    std::remove(big_input.c_str()); std::remove(big_sorted.c_str());

    typedef std::array<std::uint64_t, 4> Wide;
    const std::string wide_input = std::string(root) + "/big-wide.bin";
    const std::string wide_sorted = std::string(root) + "/big-wide-sorted.bin";
    std::vector<Wide> wide;
    for(std::uint64_t i = 0; i < (std::uint64_t(3) << 19); i++)
    {
      std::uint64_t h = (1181783497276652981ULL * (i % 700001 + 1)) ^ (i % 700001);
      wide.push_back(Wide{ h % 97, h % 1009, h, h >> 7 });
    }
    {
      std::ofstream out(wide_input.c_str(), std::ios_base::binary);
      out.write(reinterpret_cast<const char*>(wide.data()), wide.size() * sizeof(Wide));
    }
    auto compare_wide = [](const void* left, const void* right)
    {
      Wide a, b;
      std::memcpy(&a, left, sizeof(a)); std::memcpy(&b, right, sizeof(b));
      return (a < b ? -1 : (b < a ? 1 : 0));
    };
    ExternalFixedRecordSortStats wide_stats;
    ExternalFixedRecordSorter::sort(wide_input, wide_sorted, sizeof(Wide),
      256 * MEGABYTE, 2, compare_wide, &wide_stats, true);
    require(wide_stats.runs == 1);
    std::sort(wide.begin(), wide.end());
    require(readBytes(wide_sorted) == std::string(reinterpret_cast<const char*>(wide.data()),
      wide.size() * sizeof(Wide)));
    std::remove(wide_input.c_str()); std::remove(wide_sorted.c_str());
  }

  // A width the in-place sorter cannot handle falls back to the permutation.
  {
    const std::string odd_input = std::string(root) + "/odd.bin";
    const std::string odd_sorted = std::string(root) + "/odd-sorted.bin";
    std::vector<OddRecord> odd;
    for(std::uint32_t i = 0; i < 400; i++)
    {
      OddRecord record = {};
      record.key = (53 * i + 7) % 29;
      record.tag = static_cast<std::uint8_t>(i % 251);
      odd.push_back(record);
    }
    {
      std::ofstream odd_output(odd_input.c_str(), std::ios_base::binary);
      for(const OddRecord& record : odd)
      {
        odd_output.write(reinterpret_cast<const char*>(&record), sizeof(record));
      }
    }
    const size_type odd_budget =
      ExternalFixedRecordSorter::minimumBudget(sizeof(OddRecord)) * 3;
    ExternalFixedRecordSorter::sort(odd_input, odd_sorted, sizeof(OddRecord),
      odd_budget, 2, compareOddRecord, nullptr, true);
    std::sort(odd.begin(), odd.end(), [](const OddRecord& a, const OddRecord& b)
      { return compareOddRecord(&a, &b) < 0; });
    std::ifstream odd_input_stream(odd_sorted.c_str(), std::ios_base::binary);
    OddRecord observed;
    for(size_type i = 0; i < odd.size(); i++)
    {
      require(static_cast<bool>(odd_input_stream.read(
        reinterpret_cast<char*>(&observed), sizeof(observed))));
      require(observed.key == odd[i].key && observed.tag == odd[i].tag);
    }
    require(!odd_input_stream.read(reinterpret_cast<char*>(&observed), sizeof(observed)));
    std::remove(odd_input.c_str()); std::remove(odd_sorted.c_str());
  }

  bool rejected = false;
  try
  {
    ExternalFixedRecordSorter::sort(input, sorted, sizeof(Record), budget - 1, 2,
      compareRecord);
  }
  catch(const std::runtime_error&) { rejected = true; }
  require(rejected);

  std::remove(input.c_str()); std::remove(sorted.c_str()); std::remove(reduced.c_str());
  rmdir(root);
  return 0;
}
