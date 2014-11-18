CXX      := g++
BINS     := nbody
CXXFLAGS := -O3 -g3 -std=c++11 -pedantic -Wall -Wextra

all: $(BINS)

visual: CXXFLAGS += -DVISUAL
visual: $(BINS)

.SECONDEXPANSION:
$(BINS): $$@.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	-rm $(BINS) output*.bmp

video:
	avconv -r 25 -i output_%05d.bmp -vcodec qtrle -pix_fmt rgb24 -r 25 $(filter-out $@,$(MAKECMDGOALS))

%:
	@:

.PHONY: clean video
