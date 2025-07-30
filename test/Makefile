CC = gcc
CFLAGS = -I../inc -Wall -O3
LDFLAGS = -lgsl -lgslcblas -lm

BIN_DIR = build
DATA_DIR = $(BIN_DIR)/test_data

.PHONY: all test clean

all: $(BIN_DIR)/generate_test_data $(BIN_DIR)/test_lp_solver

$(BIN_DIR)/generate_test_data: generate_test_data.c test_utils.c ../src/utils.c | $(BIN_DIR)
	@$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR)/test_lp_solver: test_lp_solvers.c test_utils.c ../src/ipm.c ../src/utils.c | $(BIN_DIR)
	@$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)
	@mkdir -p $(DATA_DIR)

test: all
	@echo "Running Python test suite..."
	@python3 run_tests.py

clean:
	@rm -rf $(BIN_DIR)
