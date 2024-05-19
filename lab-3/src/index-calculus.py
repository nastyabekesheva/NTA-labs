import Factorizer
import numpy as np

def is_subset(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    return set1 <= set2

class IC:
    def __init__(self, a, b, n):
        self.n = int(n)
        self.a = int(a)
        self.b = int(b)

    def solve(self):
        self.base = self.__generate_factor_base__()


    def __generate_factor_base__(self):
        c = 3.38
        B = c * np.exp(0.5 * np.sqrt(np.log2(self.n) * np.log2(np.log2(self.n))))

        prime = [True for _ in range(B + 1)]
        p = 2

        while (p * p <= B):
            if prime[p]:
                for i in range(p * p, B + 1, p):
                    prime[i] = False
            p += 1

        prime_numbers = [p for p in range(2, B + 1) if prime[p]]
        return prime_numbers
    
    def __smooth_check__(self, k):
        ak = pow(self.alpha, k, self.n)
        factors = Factorizer.main(ak)

        if is_subset(factors, self.base):
            #return factors

            eq = np.zeros(len(self.base) + 1)

            for i in range(len(self.base)):
                eq[i] = factors.count(self.base[i])

            eq[-1] = k

            return eq







