# File       : Makefile
# Description: Makefile utility to compile ising model sequential program
# Copyright 2022 Harvard University. All Rights Reserved.
CXX ?= g++
CXXFLAGS = -g -O0 -DNDEBUG -Wall -Wextra -Wpedantic -Isrc
SRC = ising_model.cpp
DEP = ising_model.h
.PHONY: clean

all: sequential

sequential: $(SRC) $(DEP)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f sequential