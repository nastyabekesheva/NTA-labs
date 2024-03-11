import random
import math

class Test:
    def __init__(self, p):
        self.p = p
    
    def test(self):
        x = random.randint(1, self.p)
        if math.gcd(x, self.p) > 1:
            return False
        