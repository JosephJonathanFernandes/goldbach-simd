CXX      ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wno-format
LDFLAGS  ?= -pthread

ifeq ($(OS),Windows_NT)
  EXT := .exe
  LDFLAGS :=
else
  EXT :=
endif

.PHONY: all test bench clean

all: goldbach$(EXT) baseline$(EXT) research_suite$(EXT) carry_exact$(EXT)

goldbach$(EXT): goldbach.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

research_suite$(EXT): research_suite.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

carry_exact$(EXT): carry_exact.cpp
	$(CXX) $(CXXFLAGS) $< -o $@ $(LDFLAGS)

baseline$(EXT): baseline.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

test: goldbach$(EXT)
	./goldbach$(EXT) 4 2
	./goldbach$(EXT) 100 2
	./goldbach$(EXT) 1048576 2
	./goldbach$(EXT) 1e8

bench: goldbach$(EXT)
	./goldbach$(EXT) 1e9 1
	./goldbach$(EXT) 1e10

clean:
	rm -f goldbach goldbach.exe baseline baseline.exe *.ckpt
