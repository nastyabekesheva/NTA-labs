from Factorizer import Factorizer
from indexcalculus import IC
import sys
import time

if __name__ == "__main__":
    factors = []
    try:
        a, b, n, type = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
        print("Solving...")
        start_time = time.time()
        ic = IC(a, b, n, type)
        result = ic.solve()
        end_time = time.time()

        print("Solution found.")
        print("x = ", result)
        elapsed_time = end_time - start_time
        print("Elapsed time:", "{:.6f}".format(elapsed_time), "seconds")

    except ValueError:
        print("Specified parameter is not an int.")
