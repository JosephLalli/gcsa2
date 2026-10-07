#include <gcsa/external_sort.h>

#include <algorithm>
#include <array>
#include <cstdio>
#include <cstdlib>
#include <cstring>
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

typedef std::array<std::uint8_t, 5> OddRecord;

static_assert(sizeof(OddRecord) == 5, "fallback record must not be word-aligned");

int compareOddRecord(const void* left, const void* right)
{
  const auto* a = static_cast<const std::uint8_t*>(left);
  const auto* b = static_cast<const std::uint8_t*>(right);
  for(size_type i = 0; i < sizeof(OddRecord); i++)
  {
    if(a[i] != b[i]) { return (a[i] < b[i] ? -1 : 1); }
  }
  return 0;
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
  const char* configured = std::getenv("TMPDIR");
  std::string pattern = std::string(
    configured != nullptr && configured[0] != '\0' ? configured : "/tmp") +
    "/gcsa-fixed-sort-XXXXXX";
  std::vector<char> root(pattern.begin(), pattern.end()); root.push_back('\0');
  require(mkdtemp(root.data()) != nullptr);
  TempFile::setDirectory(root.data());
  const std::string input = std::string(root.data()) + "/input.bin";
  const std::string sorted = std::string(root.data()) + "/sorted.bin";
  const std::string reduced = std::string(root.data()) + "/reduced.bin";

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

  {
    const std::string permuted = std::string(root.data()) + "/permuted.bin";
    const std::string in_place = std::string(root.data()) + "/in-place.bin";
    ExternalFixedRecordSortStats permuted_stats, in_place_stats;
    ExternalFixedRecordSorter::sort(input, permuted, sizeof(Record), budget, 2,
      compareRecord, &permuted_stats, false);
    ExternalFixedRecordSorter::sort(input, in_place, sizeof(Record), budget, 2,
      compareRecord, &in_place_stats, true);
    require(readBytes(permuted) == readBytes(in_place));
    require(permuted_stats.runs > 2);
    require(in_place_stats.runs < permuted_stats.runs);
    require(in_place_stats.max_bytes_resident <= budget);

    const std::string permuted_reduced = std::string(root.data()) + "/permuted-reduced.bin";
    const std::string in_place_reduced = std::string(root.data()) + "/in-place-reduced.bin";
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

  {
    const std::string u64_input = std::string(root.data()) + "/u64.bin";
    const std::string comparator_sorted = std::string(root.data()) + "/u64-comparator.bin";
    const std::string scalar_sorted = std::string(root.data()) + "/u64-scalar.bin";
    std::vector<std::uint64_t> values;
    for(std::uint64_t i = 0; i < 5000; i++)
    {
      values.push_back((1181783497276652981ULL * (i + 1)) ^ (i << 17));
    }
    {
      std::ofstream output(u64_input.c_str(), std::ios_base::binary);
      output.write(reinterpret_cast<const char*>(values.data()),
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
    ExternalFixedRecordSorter::sort(u64_input, comparator_sorted,
      sizeof(std::uint64_t), u64_budget, 2, compare_u64, nullptr, true);
    ExternalFixedRecordSorter::sort(u64_input, scalar_sorted,
      sizeof(std::uint64_t), u64_budget, 2, compare_u64, nullptr, true,
      ExternalFixedRecordSorter::RecordOrder::ASCENDING_U64);
    require(readBytes(comparator_sorted) == readBytes(scalar_sorted));
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    std::sort(values.begin(), values.end());
    require(readBytes(scalar_sorted) == std::string(
      reinterpret_cast<const char*>(values.data()), values.size() * sizeof(std::uint64_t)));
#endif
    std::remove(u64_input.c_str()); std::remove(comparator_sorted.c_str());
    std::remove(scalar_sorted.c_str());
  }

  {
    const std::string odd_input = std::string(root.data()) + "/odd.bin";
    const std::string odd_sorted = std::string(root.data()) + "/odd-sorted.bin";
    std::vector<OddRecord> odd;
    for(std::uint32_t i = 0; i < 400; i++)
    {
      std::uint32_t key = (53 * i + 7) % 29;
      OddRecord record = {{
        static_cast<std::uint8_t>(key >> 24),
        static_cast<std::uint8_t>(key >> 16),
        static_cast<std::uint8_t>(key >> 8),
        static_cast<std::uint8_t>(key),
        static_cast<std::uint8_t>(i % 251)
      }};
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
      { return compareOddRecord(a.data(), b.data()) < 0; });
    std::ifstream odd_input_stream(odd_sorted.c_str(), std::ios_base::binary);
    OddRecord observed;
    for(size_type i = 0; i < odd.size(); i++)
    {
      require(static_cast<bool>(odd_input_stream.read(
        reinterpret_cast<char*>(observed.data()), sizeof(observed))));
      require(observed == odd[i]);
    }
    require(!odd_input_stream.read(reinterpret_cast<char*>(observed.data()), sizeof(observed)));
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
  rmdir(root.data());
  return 0;
}
