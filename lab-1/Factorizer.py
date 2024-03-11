import math 
import sympy
import numpy as np

def xor(a, b):
    return np.logical_xor(a, b).astype(int)

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
        
    def __brillhart_morrison__(self):
        base = self.__generate_factor_base__()
        cf = self.__generate_continued_fraction__(len(base))
        sn = self.__generate_smooth_numbers__(cf)
        s = [self.__find_prime_degree__(base, sn[i]**2) for i in range(len(sn))]
        x, y = 1, 1
        for i in range(len(base)):
            x *= sn[i]**s[i]
        for i in base:
            y *= sn[i]**s[i]
        
        
    def __generate_factor_base__(self):
        prime = [True] * (self.n + 1)
        i = 2
        while i**2 <= self.n:
            if prime[i] == True:        
                for j in range(i**2, self.n + 1, i): 
                    prime[j] = False
            i += 1

        L = sympy.exp((sympy.log(self.n) / sympy.log(sympy.log(self.n)))**0.5)**(1 / math.sqrt(2))
        primes = []
        p = 2
        while p < L**(1/np.sqrt(2)):
            if prime[p] and sympy.legendre_symbol(self.n, p) == 1:
                primes.append(p)
            p += 1
        return primes
    
    def __generate_continued_fraction__(self, k):
        v = [1]
        alpha = [math.sqrt(self.n)]
        a = [int(math.sqrt(self.n))]
        u = [int(math.sqrt(self.n))]

        for i in range(k):
            v.append((self.n-u[i]**2) / v[i])
            alpha.append((math.sqrt(self.n)+u[i]) / v[i+1])
            a.append(int(alpha[i]))
            u.append(a[i+1]*v[i+1]-u[i])
        
        return a
    
    def __generate_smooth_numbers__(self, a):
        b = [1, 0]
        for i in range(len(a)):
            b.append(a[i]*b[i-1]+b[i-2])

        return b
    
    def __solve_slr__(self, s):
        m, n = s.shape
        full_solution = []

        undetermined_rows = np.ones((m, 1), dtype=bool)
        while np.any(undetermined_rows):
            vector = np.zeros((n, 1), dtype=int)
            solution = []

            for i in range(n):
                for x in range(m):
                    if undetermined_rows[x] and vector[i] == s[x, i] == 1:
                        vector = xor(vector, s[x])
                        solution.append(s[x])

            if np.all(vector == 0):
                full_solution.append(solution)

            undetermined_rows = np.any(s[undetermined_rows] == 1, axis=1)

        return full_solution
    
    def __find_prime_degree__(self, factor_base, target_number):
        # Construct the matrix A
        A = np.array(factor_base).reshape(-1, 1)

        # Construct the vector b
        b = np.array([target_number])

        # Solve the linear equation Ax = b
        x, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

        # Round the degrees to integers
        x = np.round(x).astype(int)

        # Ensure non-negative values for x
        x[x < 0] = 0

        return x
