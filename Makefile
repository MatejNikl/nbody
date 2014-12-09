CXX      := g++
BIN      := nbody

CXXFLAGS := -O3 -g3 -std=c++11 -pedantic -Wall -Wextra -march=native -mtune=native
CXXFLAGS += -rdynamic -ldl#-DACCURATE_VEC

BUILD    := build
SRC      := src

SRCS     := $(sort $(wildcard $(SRC)/*.cpp))
OBJS     := $(addprefix $(BUILD)/,$(notdir $(SRCS:.cpp=.o)))
DEPS     := $(OBJS:.o=.d)

TEST     := test
REF_CONF := $(TEST)/rconfig.txt
REF_IN   := $(TEST)/rin.txt
REF_OUT  := $(TEST)/rout.txt
TEST_OUT := $(TEST)/out.txt

SRCS     := $(sort $(wildcard $(SRC)/*.cpp))
OBJS     := $(addprefix $(BUILD)/,$(notdir $(SRCS:.cpp=.o)))
VOBJS    := $(addprefix $(BUILD)/v_,$(notdir $(SRCS:.cpp=.o)))
DEPS     := $(OBJS:.o=.d) $(VOBJS:.o=.d)

# first target is the default one and it should be something reasonable,
# not the one from included dependencies
all: $(BIN)

# include compiler-generated dependencies, so obj files get recompiled when
# included headers change
-include $(DEPS)

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

test: $(BIN)
	./$(BIN) $(REF_CONF) $(REF_IN)
	diff -u $(REF_OUT) $(TEST_OUT)

showsims: $(BIN)
	nm -C $(BIN) | sed -n 's/.*\ssimulator_\(.*\)/\1/p'

.PHONY: all clean run test showsims
