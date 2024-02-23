import math 

class Factorizer:
    def __init__(self, n, f=None, method="rho-pollard", patience=50):
        self.n = n
        self.f = f
        self.mathod = method
        self.patience = patience

    def factorize(self):
        if self.method == "rho-pollard":
            if self.f:
                orbit = self.__create_orbit__(2)
                d = self.__rho_pollard__(orbit)
                i = 1
                while d == None:
                    orbit = self.__create_orbit__(2+i)
                    d = self.__rho_pollard__(orbit)
                    i += 1
                    if i == self.patience:
                        return None, None
                    
                return d, int(self.n/d)
            
            return None, None
        elif self.method == "trial-division":
            primes = self.__generate_primes__()
            
            for p in primes:
                d = self.__trial_division__(p)
                if d != None:
                    return d, int(self.n/d)
                
            return None, None


    def __create_orbit__(self, initial_x):
        orbit = [initial_x]
        new_x = self.f(initial_x, self.n)
        i = 1
        while new_x not in orbit:
            orbit.append(new_x)
            new_x = self.f(orbit[i], self.n)
            i += 1

        return orbit
    
    def __rho_pollard__(self, orbit):
        for k in range(len(orbit)):
            if 2*k < len(orbit):
                d = math.gcd(abs(orbit[2*k]-orbit[k]), self.n)
                #print(f"x_{2*k}, x_{k}: gcd({abs(orbit[2*k]-orbit[k])}, {n}) = {d}")
                if d != 1 and 2*k != k:
                    return d
            
        return None
    
    def __generate_primes__(self):
        # sieve of eratosthenes
        prime = [True] * (self.n + 1)
        i = 2
        while i**2 <= self.n:
            if prime[i] == True: # if prime[p] is not changed, then it is a prime         
                for j in range(i**2, self.n + 1, i): # update all multiples of p
                    prime[j] = False
            i += 1

        primes = [p for p in range(2, self.n + 1) if prime[p]]
        return primes
    
    def __trial_division__(self, m):   
        digits = [int(digit) for digit in str(self.n)]
        n = 0

        r = [0]
        for i in range(0, len(digits)):
            r.append(r[0] * 10 % m)

        for i in range(0, len(digits)):
            n += digits[i] * r[i] % m

        if n == self.n:
            return m
