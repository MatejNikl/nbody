CXX      := g++
BIN      := nbody
CXXFLAGS := -O0 -g3 -std=c++11 -pedantic -Wall -Wextra -fopenmp -march=native -mtune=native

BUILD    := build
SRC      := src

SRCS     := $(sort $(wildcard $(SRC)/*.cpp))
OBJS     := $(addprefix $(BUILD)/,$(notdir $(SRCS:.cpp=.o)))
DEPS     := $(OBJS:.o=.d)

all: $(BIN)

# include compiler-generated dependencies, so obj files get recompiled when
# included headers change
-include $(DEPS)

visual: CXXFLAGS += -DVISUAL
visual: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BUILD)/%.o: $(SRC)/%.cpp Makefile | $(BUILD)
	$(CXX) -c $(CXXFLAGS) -MMD $< -o $@

$(BUILD):
	mkdir $(BUILD)

clean:
	rm -rf $(BIN) $(BUILD)

run: $(BIN)
	./$(BIN) 1000 1000

.PHONY: clean run
