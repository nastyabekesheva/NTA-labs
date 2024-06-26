from Factorizer import Factorizer
import sys
import time

if __name__ == "__main__":
    factors = []
    flag = True
    try:
        a, b, p = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
        print("Solving...")
        start_time = time.time()
        F = Factorizer(a, b, p, method="sph")
        result = F.factorize()
        end_time = time.time()

        print("Solution found.")
        print("x = ", result)
        elapsed_time = end_time - start_time
        print("Elapsed time:", "{:.6f}".format(elapsed_time), "seconds")

    except ValueError:
        print("Specified parameter is not an int.")

    