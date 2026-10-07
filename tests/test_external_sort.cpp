#include <gcsa/external_sort.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
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
