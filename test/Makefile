BIN_DIR = build
DATA_DIR = $(BIN_DIR)/test_data
SRC_DIR = ../src
INC_DIR = ../inc

CC = gcc
CFLAGS = -I../inc -Wall -O3
LDFLAGS = -lm

SRC = $(wildcard $(SRC_DIR)/*.c)

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

.PHONY: all test clean

all: $(BIN_DIR)/generate_test_data $(BIN_DIR)/test_lp_solver

$(BIN_DIR)/generate_test_data: generate_test_data.c test_utils.c $(SRC) | $(BIN_DIR)
	$(V)$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/test_lp_solver: test_lp_solvers.c test_utils.c $(SRC) | $(BIN_DIR)
	$(V)$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR):
	$(V)mkdir -p $(BIN_DIR)
	$(V)mkdir -p $(DATA_DIR)

test: all
	$(V)echo "Running Python test suite..."
	$(V)python3 run_tests.py

clean:
	$(V)rm -rf $(BIN_DIR)
