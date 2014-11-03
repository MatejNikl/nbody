CXX      := g++
BINS     := nbody
CXXFLAGS := -O3 -g3 -ansi -pedantic -Wall -Wextra

all: $(BINS)

visual: CXXFLAGS += -DVISUAL
visual: $(BINS)

.SECONDEXPANSION:
$(BINS): $$@.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm $(BINS) output*.bmp

.PHONY: clean
