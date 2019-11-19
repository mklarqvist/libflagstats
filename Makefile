###################################################################
# Copyright (c) 2019
# Author(s): Marcus D. R. Klarqvist
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.
###################################################################

OPTFLAGS  := -O3 -march=native
CFLAGS     = -std=c99 $(OPTFLAGS) $(DEBUG_FLAGS)
CPPFLAGS   = -std=c++0x $(OPTFLAGS) $(DEBUG_FLAGS)
CPP_SOURCE = benchmark/flagstats.cpp benchmark/utility.cpp benchmark/generate.cpp linux/instrumented_benchmark.cpp
C_SOURCE   = 
OBJECTS    = $(CPP_SOURCE:.cpp=.o) $(C_SOURCE:.c=.o)

POSPOPCNT_PATH  := ./libalgebra
LZ4_PATH :=
ZSTD_PATH :=
INCLUDE_PATHS := -I$(PWD)
LIBRARY_PATHS :=
ifneq ($(LZ4_PATH),)
	INCLUDE_PATHS += -I$(LZ4_PATH)/include
	LIBRARY_PATHS += -L$(LZ4_PATH)/lib
endif
ifneq ($(ZSTD_PATH),)
	INCLUDE_PATHS += -I$(ZSTD_PATH)/include
	LIBRARY_PATHS += -L$(ZSTD_PATH)/lib
endif

# dedup
INCLUDE_PATHS := $(sort $(INCLUDE_PATHS))
LIBRARY_PATHS := $(sort $(LIBRARY_PATHS))

# Default target
all: benchmark utility generate

# Generic rules
utility: benchmark/utility.cpp
	$(CXX) $(CPPFLAGS) -o $@ $<

generate: benchmark/generate.cpp
	$(CXX) $(CPPFLAGS) -o $@ $<

bench.o: benchmark/flagstats.cpp
	$(CXX) $(CPPFLAGS) -I$(POSPOPCNT_PATH) $(INCLUDE_PATHS) -c -o $@ $<

benchmark: bench.o
	$(CXX) $(CPPFLAGS) bench.o -I$(POSPOPCNT_PATH) $(INCLUDE_PATHS) $(LIBRARY_PATHS) -o bench -llz4 -lzstd

instrumented_benchmark : linux/instrumented_benchmark.cpp
	$(CXX) $(CPPFLAGS) -I. -I$(POSPOPCNT_PATH) -I./linux -o $@ $<

instrumented_benchmark_align64 : linux/instrumented_benchmark.cpp
	$(CXX) $(CPPFLAGS) -DALIGN -I. -I$(POSPOPCNT_PATH) -I./linux -c -o $@ $<

clean:
	rm -f $(OBJECTS)
	rm -f bench bench.o utility generate instrumented_benchmark instrumented_benchmark_align64

.PHONY: all clean
