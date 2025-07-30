CC = gcc
CFLAGS = -I../inc -Wall -O3
LDFLAGS = -lgsl -lgslcblas -lm

all: generate_test_data test_lp_solver

generate_test_data: generate_test_data.c test_utils.c ../src/utils.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

test_lp_solver: test_lp_solvers.c test_utils.c ../src/ipm.c ../src/utils.c
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f generate_test_data test_lp_solver *.o