import random
import math
import sympy

class Test:
    def __init__(self, p, k):
        self.p = int(p)
        self.k = k
    
    def test(self):

        if self.p == 2:
            return True
        if self.p % 2 == 0 or self.p == 1:
            return False
    
        for _ in range(self.k):
            x = random.randint(2, self.p - 1)
            if sympy.gcd(x, self.p) != 1:
                return False

            jacobi = sympy.jacobi_symbol(x, self.p)
            if pow(x, (self.p - 1) // 2, self.p) != jacobi % self.p:
                return False
        
        return True

T = Test(7, 100)
print(T.test())