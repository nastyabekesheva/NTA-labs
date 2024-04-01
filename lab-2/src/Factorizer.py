import math 
import sympy
import numpy as np
import time

class Factorizer:
    def __init__(self, p, a, b, method="brute-force"):
        self.p = int(p)
        self.a = int(a)
        self.b = int(b)
        self.method = method

    def factorize(self):
        if self.method == "brute-force":
            return self.__brute_force__()
        
    def __brute_force__(self):
        x = 1
        start_time = time.time()

        while x < self.p and  time.time() - start_time < 300:
            tmp = pow(self.a, x, self.p)

            if tmp == self.b:
                print(x)
                return x
            
            x += 1
                
F = Factorizer(863, 405, 162)
print(F.factorize())