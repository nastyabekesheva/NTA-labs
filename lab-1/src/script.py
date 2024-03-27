from Factorizer import Factorizer
from Test import Test
import sys
import time

def print_factors_multiplication(factors, N):
    product = 1
    for num in factors:
        product *= num
    if int(N / product)  > 0:
        factors_str = ' * '.join(map(str, factors))
        print(f"{N} = " + factors_str + f" * {int(N / product)}")
    else:
        factors_str = ' * '.join(map(str, factors))
        print(f"{N} = " + factors_str)

if __name__ == "__main__":
    factors = []
    flag = True
    original_N = sys.argv[1]
    N = sys.argv[1]
    if N and float(N).is_integer():
        try:
            N = int(N)
            original_N = int(original_N)
        except ValueError:
            print("Specified parameter is not an int.")
            flag = False

        print(f"Factorizing: {N}")

        start_time = time.time()
        
        T = Test(N, 100)

        if T.test() == True and N == 1:
            #print(f'{N} is prime.')
            print_factors_multiplication([N], original_N)
            flag = False
        
        
        while flag:
            td = Factorizer(N, method="trial-division")
            x = td.factorize()
            if x:
                #print(f"{N} = {x} * {y}.")
                factors.append(x)
                N //= x
            else:
                flag = False

        print_factors_multiplication(factors, original_N)
        flag = True
        rp = Factorizer(N, method="rho-pollard")
        x = rp.factorize()
        if x:
            factors.append(x)
            N //= x
            T = Test(N, 100)
            if T.test() == True and N == 1:
                print(f'{N} is prime.')
                print_factors_multiplication(factors, original_N)
                flag = False

        if  N == 1:
                print(f'{N} is prime.')
                print_factors_multiplication(factors, original_N)
                flag = False
        
        while flag:

            print_factors_multiplication(factors, original_N)

            bm = Factorizer(N, method="rho-pollard")
            x = bm.factorize()
            if x:
                #print(f"{N} = {x} * {y}.")
                factors.append(x)
                N //= x

                T1 = Test(x, 1000)
                T2 = Test(N, 1000)

                if T1.test() == True and T2.test() == True and x == 1 and N == 1:
                    #factors.append(N)
                    end_time = time.time()
                    elapsed_time = end_time - start_time

                    print("Elapsed time:", elapsed_time, "seconds")

                    flag = False
                    


            else:
                end_time = time.time()
                elapsed_time = end_time - start_time

                print("Elapsed time:", elapsed_time, "seconds")
                flag = False
    else:
        print("Specified parameter is not an int.")
