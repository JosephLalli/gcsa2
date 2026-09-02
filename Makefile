SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src

# This enables various debugging options in build_gcsa.
#VERIFY_FLAGS=-DVERIFY_CONSTRUCTION

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info The compiler is Apple Clang that needs libomp for OpenMP support.)

        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # Find libomp installed by Homebrew or MacPorts.
        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found a keg-only libomp installed by Homebrew at $(HOMEBREW_PREFIX)/opt/libomp.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by MacPorts at /opt/local.)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts.)
        endif

        # We also need to link it.
        LIBS += -lomp
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(VERIFY_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -Iinclude -I$(INC_DIR)

HEADERS=$(wildcard include/gcsa/*.h)
LIBOBJS=$(addprefix $(BUILD_OBJ)/,algorithms.o checkpoint.o dbg.o disk_array.o external_join.o external_preprocessing.o external_sort.o files.o final_events.o gcsa.o internal.o lcp.o path_graph.o support.o utils.o resources.o workspace.o)
LIBRARY=$(BUILD_LIB)/libgcsa2.a

PROGRAMS=$(addprefix $(BUILD_BIN)/,build_gcsa convert_graph gcsa_format try_extend)
OBSOLETE=build_gcsa convert_graph gcsa_format try_extend

.PHONY: all clean directories test workspace-test external-path-sort-test external-join-test external-sort-test external-preprocessing-test final-events-test lcp-streaming-test disk-array-test internal-buffer-test parameter-test construction-resume-test
all: directories $(LIBRARY) $(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)

$(BUILD_BIN):
	mkdir -p $@

$(BUILD_LIB):
	mkdir -p $@

$(BUILD_OBJ):
	mkdir -p $@

$(BUILD_OBJ)/%.o:$(SOURCE_DIR)/%.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

$(BUILD_BIN)/%:$(BUILD_OBJ)/%.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

$(BUILD_OBJ)/test_workspace.o:tests/test_workspace.cpp include/gcsa/resources.h include/gcsa/workspace.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

workspace-test: directories $(BUILD_OBJ)/test_workspace.o $(BUILD_OBJ)/resources.o $(BUILD_OBJ)/workspace.o
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_workspace $(BUILD_OBJ)/test_workspace.o $(BUILD_OBJ)/resources.o $(BUILD_OBJ)/workspace.o -pthread
	$(BUILD_BIN)/test_workspace

$(BUILD_OBJ)/test_external_path_sort.o:tests/test_external_path_sort.cpp include/gcsa/path_graph.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

external-path-sort-test: directories $(BUILD_OBJ)/test_external_path_sort.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_external_path_sort $(BUILD_OBJ)/test_external_path_sort.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_external_path_sort

$(BUILD_OBJ)/test_external_sort.o:tests/test_external_sort.cpp include/gcsa/external_sort.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

external-sort-test: directories $(BUILD_OBJ)/test_external_sort.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_external_sort $(BUILD_OBJ)/test_external_sort.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_external_sort

$(BUILD_OBJ)/test_external_preprocessing.o:tests/test_external_preprocessing.cpp include/gcsa/external_preprocessing.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

external-preprocessing-test: directories $(BUILD_OBJ)/test_external_preprocessing.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_external_preprocessing $(BUILD_OBJ)/test_external_preprocessing.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_external_preprocessing

$(BUILD_OBJ)/test_final_events.o:tests/test_final_events.cpp include/gcsa/final_events.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

final-events-test: directories $(BUILD_OBJ)/test_final_events.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_final_events $(BUILD_OBJ)/test_final_events.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_final_events

$(BUILD_OBJ)/test_lcp_streaming.o:tests/test_lcp_streaming.cpp include/gcsa/lcp.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

lcp-streaming-test: directories $(BUILD_OBJ)/test_lcp_streaming.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_lcp_streaming $(BUILD_OBJ)/test_lcp_streaming.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_lcp_streaming

$(BUILD_OBJ)/test_disk_array.o:tests/test_disk_array.cpp include/gcsa/disk_array.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

disk-array-test: directories $(BUILD_OBJ)/test_disk_array.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_disk_array $(BUILD_OBJ)/test_disk_array.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_disk_array

$(BUILD_OBJ)/test_external_join.o:tests/test_external_join.cpp include/gcsa/path_graph.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

external-join-test: directories $(BUILD_OBJ)/test_external_join.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_external_join $(BUILD_OBJ)/test_external_join.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_external_join

$(BUILD_OBJ)/test_internal_buffers.o:tests/test_internal_buffers.cpp include/gcsa/internal.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

internal-buffer-test: directories $(BUILD_OBJ)/test_internal_buffers.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_internal_buffers $(BUILD_OBJ)/test_internal_buffers.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_internal_buffers

$(BUILD_OBJ)/test_parameters.o:tests/test_parameters.cpp include/gcsa/support.h include/gcsa/utils.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

parameter-test: directories $(BUILD_OBJ)/test_parameters.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_parameters $(BUILD_OBJ)/test_parameters.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_parameters

$(BUILD_OBJ)/test_construction_resume.o:tests/test_construction_resume.cpp include/gcsa/checkpoint.h
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

construction-resume-test: directories $(BUILD_OBJ)/test_construction_resume.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $(BUILD_BIN)/test_construction_resume $(BUILD_OBJ)/test_construction_resume.o $(LIBRARY) $(LIBS)
	$(BUILD_BIN)/test_construction_resume

test: workspace-test external-path-sort-test external-join-test external-sort-test external-preprocessing-test final-events-test lcp-streaming-test disk-array-test internal-buffer-test parameter-test construction-resume-test

clean:
	rm -rf $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)
	rm -f *.o *.a $(OBSOLETE)
	cd benchmark && $(MAKE) clean
