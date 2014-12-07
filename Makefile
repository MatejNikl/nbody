CXX      := g++
BIN      := nbody
VBIN     := v_nbody
CXXFLAGS := -O3 -g3 -std=c++11 -pedantic -Wall -Wextra -march=native -mtune=native

BUILD    := build
SRC      := src

SRCS     := $(sort $(wildcard $(SRC)/*.cpp))
OBJS     := $(addprefix $(BUILD)/,$(notdir $(SRCS:.cpp=.o)))
VOBJS    := $(addprefix $(BUILD)/v_,$(notdir $(SRCS:.cpp=.o)))
DEPS     := $(OBJS:.o=.d) $(VOBJS:.o=.d)

TEST     := test
REF_CONF := $(TEST)/rconfig.txt
REF_IN   := $(TEST)/rin.txt
REF_OUT  := $(TEST)/rout.txt
TEST_OUT := $(TEST)/out.txt

perf performance: $(BIN)
vis visual: $(VBIN)
all: perf vis

# include compiler-generated dependencies, so obj files get recompiled when
# included headers change
-include $(DEPS)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(VBIN): $(VOBJS)
	$(CXX) $(CXXFLAGS) -DVISUAL $^ -o $@

$(BUILD)/%.o: $(SRC)/%.cpp Makefile | $(BUILD)
	$(CXX) -c $(CXXFLAGS) -MMD $< -o $@

$(BUILD)/v_%.o: $(SRC)/%.cpp Makefile | $(BUILD)
	$(CXX) -c $(CXXFLAGS) -DVISUAL -MMD $< -o $@

$(BUILD):
	mkdir $(BUILD)

clean:
	rm -rf $(BIN) $(VBIN) $(BUILD)

run: $(BIN)
	./$(BIN) 1000 1000

test: $(BIN)
	./$(BIN) $(REF_CONF) $(REF_IN)
	diff $(REF_OUT) $(TEST_OUT)


.PHONY: perf performance vis visual all clean run test
