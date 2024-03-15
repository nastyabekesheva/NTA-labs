from Factorizer import Factorizer
from Test import Test
import sys

def print_factors_multiplication(factors, N):
    factors_str = ' * '.join(map(str, factors))
    print(f"{N} = " + factors_str)

if __name__ == "__main__":
    factors = []
    flag = True
    original_N = sys.argv[1]
    N = sys.argv[1]
    if N and float(N).is_integer():
        N = int(float(N))
        original_N = int(float(original_N))

        T = Test(N, 100)

        if T.test() == True and N == 1:
            #print(f'{N} is prime.')
            print_factors_multiplication([N], original_N)
            flag = False
        

        while flag:
            td = Factorizer(N, method="trial-division")
            x, y = td.factorize()
            if x:
                #print(f"{N} = {x} * {y}.")
                factors.append(x)
                N //= x
            else:
                flag = False

        flag = True
        rp = Factorizer(N, method="rho-pollard")
        x, y = rp.factorize()
        if x:
            #print(f"{N} = {x} * {y}.")
            factors.append(x)
            N //= x
            T = Test(N, 100)
            if T.test() == True and N == 1:
                factors.append(N)
                print(f'{N} is prime.')
                print_factors_multiplication(factors, original_N)
                flag = False

        if  N == 1:
                factors.append(N)
                print(f'{N} is prime.')
                print_factors_multiplication(factors, original_N)
                flag = False

        while flag:

            print(f"{N} - BM.")
            bm = Factorizer(N, method="brillhart-morrison")
            x, y = bm.factorize()
            if x:
                #print(f"{N} = {x} * {y}.")
                factors.append(x)
                N //= x

                T1 = Test(x, 1000)
                T2 = Test(y, 1000)

                if T1.test() == True and T2.test() == True and x == 1 and y == 1:
                    factors.append(y)
                    print_factors_multiplication(factors, original_N)
                    flag = False
    else:
        print("Specified parameter is not an int.")
