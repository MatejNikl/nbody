CXX      := g++
BIN      := nbody
CXXFLAGS := -Ofast -std=c++11 -pedantic -Wall -Wextra -fopenmp -march=native -mtune=native

IMAGE    := output_00000.bmp

all: $(BIN)

visual: CXXFLAGS += -DVISUAL
visual: $(BIN)

.SECONDEXPANSION:
$(BIN): $$@.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	rm -f $(BIN)

run:
	./$(BIN) 1000 1000

visualrun: visual
	./$(BIN) 50 25000 800 600 50 0.002 0 1000 1000 0 0.001 -1000

video: $(IMAGE)
	avconv -r 25 -i output_%05d.bmp -vcodec qtrle -pix_fmt rgb24 -r 25 $(filter-out $@,$(MAKECMDGOALS))
	rm output_*.bmp

$(IMAGE): visualrun

.PHONY: clean run visualrun video
