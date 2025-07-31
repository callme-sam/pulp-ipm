
BUILD_DIR = build
OBJ_DIR = $(BUILD_DIR)/obj
SRC_DIR = src
INC_DIR = inc

SRC = main.c $(wildcard $(SRC_DIR)/*.c)
OBJ = $(patsubst %.c, $(OBJ_DIR)/%.o, $(SRC))
DEP = $(patsubst %.c, $(OBJ_DIR)/%.d, $(SRC))
TARGET = $(BUILD_DIR)/ipm_solver

CC = gcc
CFLAGS = -Wall -Wextra -I$(INC_DIR)
LDFLAGS = -lgsl -lgslcblas -lm
MAKEFLAGS += --no-print-directory

# Makefile Verbosity
VERBOSE ?= 0
ifeq ($(VERBOSE),1)
  V =
else
  V = @
endif

# Log Verbosity
LOG_LEVEL ?= 3	# default is INFO
CFLAGS += -DCURRENT_LOG_LEVEL=$(LOG_LEVEL)

# Rules

.PHONY: all clean run test

all: $(TARGET)

-include $(DEP)

# Link
$(TARGET): $(OBJ) | $(BUILD_DIR)
	$(V)echo "Linking $@"
	$(V)$(CC) -o $@ $^ $(LDFLAGS)

# Compile
$(OBJ_DIR)/%.o: %.c
	$(V)mkdir -p $(dir $@)
	$(V)echo "Compiling $<"
	$(V)$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

# Build directory
$(BUILD_DIR):
	$(V)echo "Creating directory $@"
	$(V)mkdir -p $@

# Clean
clean:
	$(V)echo "Cleaning build directory"
	$(V)rm -rf build/

# Run
run: $(TARGET)
	$(V)echo "Running $(TARGET)"
	$(V)./$(TARGET)

# Test
test:
	$(V)$(MAKE)  -C test test

test_clean:
	$(V)$(MAKE) -C test clean