from Factorizer import Factorizer
from Test import Test
import sys

if __name__ == "__main__":
    factors = []
    N = sys.argv[1]
    if N and N.is_integer():
        N = int(N)
        T = Test(N, 100)

        if T.test() == False:
            print(f'{N} is prime.')
            return [N]
        

        while True:
            td = Factorizer(N, method="trial-division")
            x, y = td.factorize()
            if x:
                print(f"{N} = {x} * {y}.")
                factors.append(x)
                N //= x
            else:
                break

        rp = Factorizer(N, mathod="rho-pollard")
        x, y = rp.factorize()
            if x:
                print(f"{N} = {x} * {y}.")
                factors.append(x)
                N //= x
                T = Test(N, 100)

                if T.test() == False:
                    print(f'{N} is prime.')
                    return factors
        while True:
            bm = Factorizer(N, mathod="brillhart-morrison")
            x, y = bm.factorize()
            if x:
                print(f"{N} = {x} * {y}.")
                factors.append(x)
                N //= x

                T1 = Test(x, 1000)
                T2 = Test(y, 1000)

                if T1.test() == False and T2.test() = False:
                    factors.append(y)
                    return factors
                    break
    else:
        print("Specified parameter is not an int.")

