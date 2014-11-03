
CXX      := g++
BINS     := nbody
CXXFLAGS := -O0 -g3

all: $(BINS)

visual: CXXFLAGS += -DVISUAL
visual: $(BINS)

.SECONDEXPANSION:
$(BINS): $$@.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm $(BINS)

.PHONY: clean
