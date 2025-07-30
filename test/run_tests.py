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
NUM_TESTS = 5
MIN_M = 2       # Minimum number of constraints
MAX_M = 5       # Maximum number of constraints
MIN_N = 6       # Minimum number of variables (n > m)
MAX_N = 10      # Maximum number of variables
TOLERANCE = 1e-3  # Tolerance for comparing results

def run_c_solver(problem_file, test_num):
    """Run the C solver and return output files"""
    c_val_file = f"test/build/test_data/c_opt_val_{test_num}.txt"
    c_vec_file = f"test/build/test_data/c_x_opt_{test_num}.txt"
    
    # Build path to C executable
    c_executable = Path(__file__).parent / "build" / "test_lp_solver"
    
    # Run C program
    cmd = f"{c_executable} {problem_file} {c_val_file} {c_vec_file}"
    subprocess.run(cmd, shell=True, check=True)
    
    return c_val_file, c_vec_file

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
        return solver.value, solver.x_opt
    else:
        return None, None

def compare_results(c_val, py_val, c_x, py_x, test_num):
    """Compare results from C and Python implementations"""
    print(f"\nTest {test_num} Results:")
    
    if c_val is not None and py_val is not None:
        print(f"Python optimal value: {py_val:.6f}")
        print(f"C optimal value:     {c_val:.6f}")
        val_diff = abs(c_val - py_val)
        print(f"Absolute difference: {val_diff:.6f}")
        print(f"Within tolerance ({TOLERANCE}): {'Yes' if val_diff < TOLERANCE else 'No'}")
    else:
        print("One or both solvers failed to find a solution")
    
    if c_x is not None and py_x is not None:
        x_diff = np.linalg.norm(c_x - py_x)
        print(f"Solution vector difference (L2 norm): {x_diff:.6f}")
        print(f"Within tolerance ({TOLERANCE}): {'Yes' if x_diff < TOLERANCE else 'No'}")
    else:
        print("Could not compare solution vectors")

def main():
    # Create test directory if it doesn't exist
    os.makedirs('test/build/test_data', exist_ok=True)
    
    # Compile C test solver
    print("Compiling C test solver...")
    compile_cmd = "make -C test"
    subprocess.run(compile_cmd, shell=True, check=True)
    
    # Run tests
    for i in range(1, NUM_TESTS + 1):
        # Generate random problem dimensions (n > m)
        m = np.random.randint(MIN_M, MAX_M + 1)
        n = np.random.randint(max(m + 1, MIN_N), MAX_N + 1)
        
        print(f"\nRunning Test {i} with m={m}, n={n}")
        
        # Generate problem file
        problem_file = f"test/build/test_data/test_{i}.bin"
        
        # Generate and save problem data
        generate_cmd = f"test/build/generate_test_data {m} {n} {problem_file}"
        subprocess.run(generate_cmd, shell=True, check=True)
        
        # Run Python solver
        print("Running Python solver...")
        py_val, py_x = run_python_solver(problem_file)
        
        # Run C solver
        print("Running C solver...")
        c_val_file, c_vec_file = run_c_solver(problem_file, i)
        
        # Load C solution
        c_val = np.loadtxt(c_val_file) if os.path.exists(c_val_file) else None
        c_x = np.loadtxt(c_vec_file) if os.path.exists(c_vec_file) else None
        
        # Compare results
        compare_results(c_val, py_val, c_x, py_x, i)

if __name__ == "__main__":
    main()