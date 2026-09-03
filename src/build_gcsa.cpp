/*
  Copyright (c) 2018, 2019 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#include <algorithm>
#include <cerrno>
#include <filesystem>
#include <fcntl.h>
#include <getopt.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

#include <gcsa/algorithms.h>
#include <gcsa/path_graph.h>

using namespace gcsa;

//------------------------------------------------------------------------------

const size_type INDENT = 20;

template<class IndexType>
bool
atomicStore(const IndexType& index, const std::string& final_name)
{
  std::string partial = final_name + "." +
    std::to_string(static_cast<unsigned long long>(getpid())) + ".partial";
  if(!sdsl::store_to_file(index, partial)) { return false; }
  int descriptor = ::open(partial.c_str(), O_RDONLY);
  if(descriptor < 0 || ::fdatasync(descriptor) != 0)
  {
    if(descriptor >= 0) { ::close(descriptor); }
    ::unlink(partial.c_str()); return false;
  }
  if(::close(descriptor) != 0 || ::rename(partial.c_str(), final_name.c_str()) != 0)
  {
    ::unlink(partial.c_str()); return false;
  }
  std::filesystem::path parent = std::filesystem::path(final_name).parent_path();
  if(parent.empty()) { parent = "."; }
  descriptor = ::open(parent.c_str(), O_RDONLY | O_DIRECTORY);
  if(descriptor < 0 || ::fsync(descriptor) != 0)
  {
    if(descriptor >= 0) { ::close(descriptor); }
    return false;
  }
  return (::close(descriptor) == 0);
}

int
main(int argc, char** argv)
{
  if(argc == 3 && std::string(argv[1]) == "gcsa-worker-task")
  {
    return externalPathJoinWorker(argv[2]);
  }
  if(argc < 2)
  {
    Version::print(std::cerr, "GCSA2 builder");
    std::cerr << "Usage: build_gcsa [options] base_name [base_name2 ..]" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Input/output options:" << std::endl;
    std::cerr << "  -b    Read the input in binary format (default)" << std::endl;
    std::cerr << "  -t    Read the input in text format" << std::endl;
    std::cerr << "  -o X  Use X as the base name for output (default: the first input)" << std::endl;
    std::cerr << "Index construction options:" << std::endl;
    std::cerr << "  -d N  Doubling steps (default " << ConstructionParameters::DOUBLING_STEPS << ", max " << ConstructionParameters::MAX_STEPS << ")" << std::endl;
    std::cerr << "  -m X  Use node mapping from file X" << std::endl;
    std::cerr << "  -s N  Use sample period N (default " << ConstructionParameters::SAMPLE_PERIOD << ")" << std::endl;
    std::cerr << "  -B N  Set LCP branching factor to N (default " << ConstructionParameters::LCP_BRANCHING << ")" << std::endl;
    std::cerr << "  -L    Load the index instead of building it" << std::endl;
    std::cerr << "  -v    Verify the index by querying it with the kmers" << std::endl;
    std::cerr << "Other options:" << std::endl;
    std::cerr << "  -D X  Use X as the directory for temporary files (default: " << TempFile::DEFAULT_TEMP_DIR << ")" << std::endl;
    std::cerr << "  -l N  Limit disk space usage to N gigabytes (default " << ConstructionParameters::SIZE_LIMIT << ")" << std::endl;
    std::cerr << "  -T N  Set the number of threads to N (default and max " << omp_get_max_threads() << " on this system)" << std::endl;
    std::cerr << "  -V N  Set verbosity level to N (default 3)" << std::endl;
    std::cerr << "External-memory options (sizes accept K/M/G/T and KiB/GiB forms):" << std::endl;
    std::cerr << "      --work-dir PATH            durable construction workspace" << std::endl;
    std::cerr << "      --resume                   resume compatible committed phases" << std::endl;
    std::cerr << "      --keep-work                retain workspace after success" << std::endl;
    std::cerr << "      --memory-limit SIZE        external working-set goal (RAM/HDD tradeoff)" << std::endl;
    std::cerr << "      --disk-limit SIZE          spill-generation disk budget" << std::endl;
    std::cerr << "      --io-buffer-size SIZE      byte size of sequential I/O buffers" << std::endl;
    std::cerr << "      --sort-run-size SIZE       maximum label-sort working set" << std::endl;
    std::cerr << "      --join-partition-size SIZE maximum join working set" << std::endl;
    std::cerr << "      --merge-fan-in N           maximum merge inputs" << std::endl;
    std::cerr << "      --max-open-files N         construction descriptor ceiling" << std::endl;
    std::cerr << "      --process-workers N        parallel external join worker processes" << std::endl;
    std::cerr << "      --checkpoint-records N     record checkpoint interval" << std::endl;
    std::cerr << "      --checkpoint-bytes SIZE    byte checkpoint interval" << std::endl;
    std::cerr << "      --verify-workspace         verify payload checksums on reuse" << std::endl;
    std::cerr << "      --clean-obsolete           remove eligible predecessor artifacts" << std::endl;
    std::cerr << "      --stop-after PHASE         stop after a committed phase" << std::endl;
    std::cerr << "      --allow-path-explosion     permit growth up to the disk ceiling" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }
  Verbosity::set(Verbosity::FULL);

  int c = 0;
  bool binary = true, load_index = false, verify = false;
  std::string index_file, lcp_file, mapping_file;
  ConstructionParameters parameters;
  enum LongOption
  {
    OPT_WORK_DIR = 1000, OPT_RESUME, OPT_KEEP_WORK, OPT_MEMORY_LIMIT,
    OPT_DISK_LIMIT, OPT_IO_BUFFER, OPT_SORT_RUN, OPT_JOIN_PARTITION,
    OPT_MERGE_FAN_IN, OPT_MAX_OPEN_FILES, OPT_CHECKPOINT_RECORDS,
    OPT_CHECKPOINT_BYTES, OPT_VERIFY_WORKSPACE, OPT_CLEAN_OBSOLETE,
    OPT_STOP_AFTER, OPT_ALLOW_PATH_EXPLOSION, OPT_PROCESS_WORKERS
  };
  static struct option long_options[] =
  {
    { "work-dir", required_argument, nullptr, OPT_WORK_DIR },
    { "resume", no_argument, nullptr, OPT_RESUME },
    { "keep-work", no_argument, nullptr, OPT_KEEP_WORK },
    { "memory-limit", required_argument, nullptr, OPT_MEMORY_LIMIT },
    { "disk-limit", required_argument, nullptr, OPT_DISK_LIMIT },
    { "io-buffer-size", required_argument, nullptr, OPT_IO_BUFFER },
    { "sort-run-size", required_argument, nullptr, OPT_SORT_RUN },
    { "join-partition-size", required_argument, nullptr, OPT_JOIN_PARTITION },
    { "merge-fan-in", required_argument, nullptr, OPT_MERGE_FAN_IN },
    { "max-open-files", required_argument, nullptr, OPT_MAX_OPEN_FILES },
    { "process-workers", required_argument, nullptr, OPT_PROCESS_WORKERS },
    { "checkpoint-records", required_argument, nullptr, OPT_CHECKPOINT_RECORDS },
    { "checkpoint-bytes", required_argument, nullptr, OPT_CHECKPOINT_BYTES },
    { "verify-workspace", no_argument, nullptr, OPT_VERIFY_WORKSPACE },
    { "clean-obsolete", no_argument, nullptr, OPT_CLEAN_OBSOLETE },
    { "stop-after", required_argument, nullptr, OPT_STOP_AFTER },
    { "allow-path-explosion", no_argument, nullptr, OPT_ALLOW_PATH_EXPLOSION },
    { "output", required_argument, nullptr, 'o' },
    { "threads", required_argument, nullptr, 'T' },
    { nullptr, 0, nullptr, 0 }
  };
  while((c = getopt_long(argc, argv, "bto:d:m:s:B:LvD:l:T:V:",
    long_options, nullptr)) != -1)
  {
    switch(c)
    {
    case 'b':
      binary = true; break;
    case 't':
      binary = false; break;
    case 'o':
      index_file = std::string(optarg) + GCSA::EXTENSION;
      lcp_file = std::string(optarg) + LCPArray::EXTENSION;
      break;
    case 'd':
      parameters.setSteps(std::stoul(optarg)); break;
    case 'm':
      mapping_file = optarg; break;
    case 's':
      parameters.setSamplePeriod(std::stoul(optarg)); break;
    case 'B':
      parameters.setLCPBranching(std::stoul(optarg)); break;
    case 'L':
      load_index = true; break;
    case 'v':
      verify = true; break;
    case 'D':
      TempFile::setDirectory(optarg); break;
    case 'l':
      parameters.setLimit(std::stoul(optarg)); break;
    case 'T':
      omp_set_num_threads(Range::bound(std::stoul(optarg), 1, omp_get_max_threads())); break;
    case 'V':
      Verbosity::set(std::stoul(optarg)); break;
    case OPT_WORK_DIR:
      parameters.setWorkDirectory(optarg); break;
    case OPT_RESUME:
      parameters.setResume(); break;
    case OPT_KEEP_WORK:
      parameters.setKeepWork(); break;
    case OPT_MEMORY_LIMIT:
      parameters.setMemoryLimitBytes(parseBytes(optarg)); break;
    case OPT_DISK_LIMIT:
      parameters.setLimitBytes(parseBytes(optarg)); break;
    case OPT_IO_BUFFER:
      parameters.setIOBufferSize(parseBytes(optarg)); break;
    case OPT_SORT_RUN:
      parameters.setSortRunSize(parseBytes(optarg)); break;
    case OPT_JOIN_PARTITION:
      parameters.setJoinPartitionSize(parseBytes(optarg)); break;
    case OPT_MERGE_FAN_IN:
      parameters.setMergeFanIn(std::stoull(optarg)); break;
    case OPT_MAX_OPEN_FILES:
      parameters.setMaxOpenFiles(std::stoull(optarg)); break;
    case OPT_PROCESS_WORKERS:
      parameters.setProcessWorkers(std::stoull(optarg)); break;
    case OPT_CHECKPOINT_RECORDS:
      parameters.setCheckpointRecords(std::stoull(optarg)); break;
    case OPT_CHECKPOINT_BYTES:
      parameters.setCheckpointBytes(parseBytes(optarg)); break;
    case OPT_VERIFY_WORKSPACE:
      parameters.setVerifyWorkspace(); break;
    case OPT_CLEAN_OBSOLETE:
      parameters.setCleanObsolete(); break;
    case OPT_STOP_AFTER:
      parameters.setStopAfter(optarg); break;
    case OPT_ALLOW_PATH_EXPLOSION:
      parameters.setAllowPathExplosion(); break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(optind >= argc)
  {
    std::cerr << "build_gcsa: No input files specified" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(index_file.empty())
  {
    index_file = std::string(argv[optind]) + GCSA::EXTENSION;
    lcp_file = std::string(argv[optind]) + LCPArray::EXTENSION;
  }
  if(parameters.getResume() && parameters.getWorkDirectory().empty())
  {
    std::cerr << "build_gcsa: --resume requires --work-dir" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  parameters.setWorkerExecutable(argv[0]);
  if(parameters.externalMemory())
  {
    std::error_code error;
    std::filesystem::create_directories(parameters.getWorkDirectory(), error);
    if(error)
    {
      std::cerr << "build_gcsa: Cannot create workspace "
                << parameters.getWorkDirectory() << ": " << error.message() << std::endl;
      std::exit(EXIT_FAILURE);
    }
    TempFile::setDirectory(parameters.getWorkDirectory());
  }

  Version::print(std::cout, "GCSA2 builder");
  for(int i = optind; i < argc; i++)
  {
    printHeader("Input", INDENT);
    std::cout << argv[i];
    if(binary) { std::cout << InputGraph::BINARY_EXTENSION << " (binary format)" << std::endl; }
    else { std::cout << InputGraph::TEXT_EXTENSION << " (text format)" << std::endl; }
  }
  if(!(mapping_file.empty()))
  {
    printHeader("Node mapping", INDENT); std::cout << mapping_file << std::endl;
  }
  printHeader("Output", INDENT); std::cout << index_file << ", " << lcp_file << std::endl;
  if(!load_index)
  {
    printHeader("Doubling steps", INDENT); std::cout << parameters.doubling_steps << std::endl;
    printHeader("Sample period", INDENT); std::cout << parameters.sample_period << std::endl;
    printHeader("Branching factor", INDENT); std::cout << parameters.lcp_branching << std::endl;
    printHeader("Temp directory", INDENT); std::cout << TempFile::temp_dir << std::endl;
    printHeader("Size limit", INDENT); std::cout << inGigabytes(parameters.size_limit) << " GB" << std::endl;
    printHeader("Memory goal", INDENT); std::cout << formatBytes(parameters.memory_limit) << std::endl;
    if(parameters.externalMemory())
    {
      printHeader("Work directory", INDENT); std::cout << parameters.getWorkDirectory() << std::endl;
      printHeader("Resume", INDENT); std::cout << (parameters.getResume() ? "yes" : "no") << std::endl;
    }
    printHeader("Threads", INDENT); std::cout << omp_get_max_threads() << std::endl;
    printHeader("Verbosity", INDENT); std::cout << Verbosity::levelName() << std::endl;
  }
  std::cout << std::endl;

  InputGraph graph(argc - optind, argv + optind, binary, parameters, Alphabet(), mapping_file);

  GCSA index;
  LCPArray lcp;
  bool stored_directly = false;
  if(load_index)
  {
    if(!sdsl::load_from_file(index, index_file))
    {
      std::cerr << "build_gcsa: Cannot load the index from " << index_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(!sdsl::load_from_file(lcp, lcp_file))
    {
      std::cerr << "build_gcsa: Cannot load the LCP array from " << lcp_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else
  {
    double start = readTimer();
    try
    {
      if(parameters.externalMemory())
      {
        GCSA::buildAndStore(graph, parameters, index_file);
        stored_directly = true;
      }
      else { index = GCSA(graph, parameters); }
      lcp = LCPArray(graph, parameters);
    }
    catch(const ConstructionStopped& stopped)
    {
      std::cout << stopped.what() << std::endl;
      return 0;
    }
    double seconds = readTimer() - start;
    std::cout << "Index built in " << seconds << " seconds" << std::endl;
    std::cout << "Memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
    std::cout << "I/O volume: " << inGigabytes(readVolume()) << " GB read, "
              << inGigabytes(writeVolume()) << " GB write" << std::endl;
    std::cout << std::endl;
    if(!stored_directly && !atomicStore(index, index_file))
    {
      std::cerr << "build_gcsa: Cannot write the index to " << index_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if(!atomicStore(lcp, lcp_file))
    {
      std::cerr << "build_gcsa: Cannot write the LCP array to " << lcp_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // Loading the just-written index solely for statistics would defeat the
  // staged packer's resident-memory reduction. Verification necessarily needs
  // the query index; otherwise report the durable output without reloading it.
  if(stored_directly && verify)
  {
    if(!sdsl::load_from_file(index, index_file))
    {
      std::cerr << "build_gcsa: Cannot reload the index for verification from "
                << index_file << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if(!stored_directly || verify) { printStatistics(index, lcp); }
  else
  {
    std::cout << "Index components stored directly in " << index_file
              << " (not reloaded for statistics)" << std::endl << std::endl;
  }

  if(verify)
  {
    if(parameters.externalMemory())
    {
      const size_type verification_budget = std::max(
        verifyIndexMinimumBudget(),
        std::min(parameters.getMemoryLimitBytes(), static_cast<size_type>(64 * MEGABYTE)));
      verifyIndex(index, &lcp, graph, verification_budget,
        parameters.getMergeFanIn());
    }
    else { verifyIndex(index, &lcp, graph); }
  }

  std::cout << "Final memory usage: " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
