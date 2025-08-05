#!/usr/bin/env python3
import numpy as np
import subprocess
import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt

# Add parent directory to path to import lp_solver
sys.path.append(str(Path(__file__).parent.parent / "python" / "lp_solver"))
from lp_solver import LPSolver

# Configuration
NUM_TESTS = 100
MIN_M = 2       # Minimum number of constraints
MAX_M = 5       # Maximum number of constraints
MIN_N = 6       # Minimum number of variables (n > m)
MAX_N = 10      # Maximum number of variables
TOLERANCE = 1e-3  # Tolerance for comparing results

def run_c_solver(problem_file, test_num):
    """Run the C solver and return output files"""
    c_sts_file = f"build/test_data/c_status_{test_num}.txt"
    c_val_file = f"build/test_data/c_opt_val_{test_num}.txt"
    c_vec_file = f"build/test_data/c_x_opt_{test_num}.txt"

    # Build path to C executable
    c_executable = Path(__file__).parent / "build" / "test_lp_solver"

    # Run C program
    cmd = f"{c_executable} {problem_file} {c_sts_file} {c_val_file} {c_vec_file}"
    subprocess.run(cmd, shell=True, check=True)

    return c_sts_file, c_val_file, c_vec_file

def run_python_solver(problem_file):
    """Run the Python solver and return solution"""
    # Load problem data
    A, b, c = [], [], []
    with open(problem_file, "rb") as f:
        # Read matrix A
        rows = np.fromfile(f, dtype=np.uint64, count=1)[0]
        cols = np.fromfile(f, dtype=np.uint64, count=1)[0]
        A = np.fromfile(f, dtype=np.double, count=rows*cols).reshape((rows, cols))

        # Read vector b
        size_b = np.fromfile(f, dtype=np.uint64, count=1)[0]
        b = np.fromfile(f, dtype=np.double, count=size_b)

        # Read vector c
        size_c = np.fromfile(f, dtype=np.uint64, count=1)[0]
        c = np.fromfile(f, dtype=np.double, count=size_c)

    # Solve with Python
    solver = LPSolver()
    solver.solve(A, b, c)

    if solver.status == 'optimal':
        return solver.status, solver.value, solver.x_opt
    else:
        return solver.status, None, None

def compare_results(c_status, py_status, c_val, py_val, c_x, py_x, test_num):
    """Compare results from C and Python implementations"""
    val_diff = None
    x_diff = None

    print(f"\nTest {test_num} Results:")

    if (c_status == py_status == 'optimal'):
        print(f"Python optimal value: {py_val:.6f}")
        print(f"C optimal value:     {c_val:.6f}")
        val_diff = abs(c_val - py_val)
        print(f"Absolute difference: {val_diff:.6f}")
        print(f"Within tolerance ({TOLERANCE}): {'Yes' if val_diff < TOLERANCE else 'No'}")

        x_diff = np.linalg.norm(c_x - py_x)
        print(f"Solution vector difference (L2 norm): {x_diff:.6f}")
        print(f"Within tolerance ({TOLERANCE}): {'Yes' if x_diff < TOLERANCE else 'No'}")


    elif (c_status == py_status == 'infeasible'):
        print(f"Both python and c solver returned infeasible!")

    else:
        print("One or both solvers failed to find a solution")
        print(f"C solver status: {c_status}")
        print(f"Py solver status: {py_status}")

    result = {
            "test_num": test_num,
            "c_status": c_status,
            "py_status": py_status,
            "val_diff": val_diff,
            "x_diff": x_diff
        }

    return result

def print_summary(results):
    """Print test summary with colored output"""
    # ANSI color codes
    GREEN = '\033[0;32m'
    RED = '\033[0;31m'
    YELLOW = '\033[1;33m'
    BOLD = '\033[1m'
    ENDC = '\033[0m'  # Resets the color

    print(f"\n{BOLD}=== Test Summary ==={ENDC}")
    total_tests = len(results)
    successful_tests = 0
    failed_tests = 0

    for res in results:
        print(f"\nTest {res['test_num']}:")

        if res['c_status'] == res['py_status'] == 'optimal':
            if res['val_diff'] < TOLERANCE and res['x_diff'] < TOLERANCE:
                print(f"{GREEN}  ✓ SUCCESS: Both solvers returned optimal solutions within tolerance{ENDC}")
                print(f"  Value difference: {res['val_diff']:.6f}")
                print(f"  Solution difference: {res['x_diff']:.6f}")
                successful_tests += 1
            else:
                print(f"{RED}  ✗ FAILURE: Solutions differ beyond tolerance{ENDC}")
                print(f"  Value difference: {res['val_diff']:.6f} (Tolerance: {TOLERANCE})")
                print(f"  Solution difference: {res['x_diff']:.6f} (Tolerance: {TOLERANCE})")
                failed_tests += 1

        elif res['c_status'] == res['py_status'] == 'infeasible':
            print(f"{GREEN}  ✓ SUCCESS: Both solvers agree the problem is infeasible{ENDC}")
            successful_tests += 1

        else:
            print(f"{RED}  ✗ FAILURE: Solvers disagree on problem status{ENDC}")
            print(f"  C solver status: {res['c_status']}")
            print(f"  Python solver status: {res['py_status']}")
            failed_tests += 1

    success_rate = successful_tests/total_tests*100
    color = GREEN if success_rate == 100 else YELLOW if success_rate >= 1 else RED
    print(f"\n{BOLD}{color}=== Final Statistics ==={ENDC}")
    print(f"{color}Total tests: {total_tests}{ENDC}")
    print(f"{color}Successful tests: {successful_tests}{ENDC}")
    print(f"{color}Failed tests: {failed_tests}{ENDC}")
    print(f"{color}Success rate: {success_rate:.1f}%{ENDC}")


def main():
    results = []

    # Run tests
    for i in range(1, NUM_TESTS + 1):
        # Generate random problem dimensions (n > m)
        m = np.random.randint(MIN_M, MAX_M + 1)
        n = np.random.randint(max(m + 1, MIN_N), MAX_N + 1)

        print(f"\nRunning Test {i} with m={m}, n={n}")

        # Generate problem file
        problem_file = f"build/test_data/test_{i}.bin"

        # Generate and save problem data
        generate_cmd = f"build/generate_test_data {m} {n} {problem_file}"
        subprocess.run(generate_cmd, shell=True, check=True)

        # Run Python solver
        print("Running Python solver...")
        py_status, py_val, py_x = run_python_solver(problem_file)

        # Run C solver
        print("Running C solver...")
        c_sts_file, c_val_file, c_vec_file = run_c_solver(problem_file, i)

        # Load C solution
        c_status = open(c_sts_file).read().strip() if os.path.exists(c_sts_file) else None
        c_val = np.loadtxt(c_val_file) if os.path.exists(c_val_file) else None
        c_x = np.loadtxt(c_vec_file) if os.path.exists(c_vec_file) else None

        # Compare results
        result = compare_results(c_status, py_status, c_val, py_val, c_x, py_x, i)

        # Save results for later
        results.append(result)

    print_summary(results)

if __name__ == "__main__":
    main()