CXX := g++
SRC_DIR := src
INCLUDE_DIR := include

SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp, $(SRC_DIR)/%.o, $(SRC_FILES))

.PHONY: all clean

all: cpp_server

cpp_server: $(OBJ_FILES) $(CXX) -o $@ $^ -I$(INCLUDE_DIR)

$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp $(CXX) -c -o $@ $< -I$(INCLUDE_DIR)

clean: rm -f cpp_server $(OBJ_FILES)
