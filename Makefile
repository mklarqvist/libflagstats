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
CPP_SOURCE = lz4_import.cpp
C_SOURCE   = 
OBJECTS    = $(CPP_SOURCE:.cpp=.o) $(C_SOURCE:.c=.o)
POSPOPCNT  := /home/mdrk/repos/FastFlagStats

# Default target
all: flagstats

# Generic rules
%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CPPFLAGS) -c -o $@ $<

pospopcnt.o: $(POSPOPCNT)/pospopcnt.c
	$(CC) $(CFLAGS) -c -o $@ $<

lz4_import.o: lz4_import.cpp $(POSPOPCNT)/pospopcnt.c
	$(CXX) $(CPPFLAGS) -I$(POSPOPCNT) -c -o $@ $<

flagstats: lz4_import.o pospopcnt.o
	$(CXX) $(CPPFLAGS) lz4_import.o pospopcnt.o -I$(POSPOPCNT) -o flagstats -llz4 -lzstd

clean:
	rm -f $(OBJECTS)
	rm -f flagstats

.PHONY: all clean
