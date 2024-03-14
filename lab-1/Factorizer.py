import math 
import sympy
import numpy as np

def xor(a, b):
    return np.logical_xor(a, b).astype(int)

def is_array_nonzero(arr):
    return np.any(arr != 0)

def replace_none_with_zero(arr, k):
    return [x if x is not None else np.zeros(k) for x in arr]

def any_non_zero(arr):
    return any(any(x) for x in arr)

def remove_none_rows(matrix):
    return [row for row in matrix if row is not None]

def check_matching_positions(arrays):
    result = np.all(arrays, axis=0)
    
    return result.any()

def find_index_of_nonzero_array(matrix):
    return np.nonzero(np.any(matrix, axis=1))[0] # Return None if no row has non-zero elements

def fx(x, n):
    return (x**2+1) % n

class Factorizer:
    def __init__(self, n, f=fx, method="rho-pollard", patience=50):
        self.n = n
        self.f = f
        self.method = method
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
        
        elif self.method == "brillhart-morrison":
            a = 1 / math.sqrt(2)
            res = None
            i = 0

            while res == None or self.patience >= i:
                x, y = self.__brillhart_morrison__(a)
                if x and y.is_integer():
                    d = np.gcd((x+int(y)), self.n)

                    if d > 1 and d < self.n:
                        return d, int(self.n/d)
                i += 1
                a += 0.5
                
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
        while i**2 <= self.n and i <= 47:
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
        
    def __brillhart_morrison__(self, a):
        base = self.__generate_factor_base__(a)
        cf = self.__generate_continued_fraction__(len(base)-1)
        sn = self.__generate_smooth_numbers__(cf)
        s = [self.__find_mod_prime_degrees__(base, pow(sn[i], 2, self.n)) for i in range(len(sn))]
        s = np.array(replace_none_with_zero(s, len(base)))
        if s.any():
            if len(find_index_of_nonzero_array(s)) > 1:
                s2 = s % 2
                xs = self.__solve_slr__(np.array(s2))
                if len(xs) > 0:
                    x, y = 1, 1
                    for i in range(len(base)):
                        x = (x * sn[i]**xs[i]) % self.n
                    pows = []
                    for i in range(len(base)-1):
                        tmp = 0
                        for j in range(len(s)):
                            #f = sum([xs[0][j]*s[i][j] for j in range(len(base))]) 
                            tmp += xs[j]*s[j][i]
                        pows.append(tmp)
                        #y = (y * base[i]**f) % self.n
                    if pows[0] % 2 != 0:
                        return None, None
                    else:
                        pows[0] = 0
                        pows = [i/2 for i in pows]
                    for i in range(len(base)-1):
                        y = (y * base[i]**pows[i]) % self.n
                    return x, y
                else:
                    return None, None
            else:
                return None, None
        else:
            return None, None
        
        
    def __generate_factor_base__(self, a):
        prime = [True] * (self.n + 1)
        i = 2
        while i**2 <= self.n:
            if prime[i] == True:        
                for j in range(i**2, self.n + 1, i): 
                    prime[j] = False
            i += 1

        L = sympy.exp((sympy.log(self.n) / sympy.log(sympy.log(self.n)))**0.5)**(a)
        primes = [-1, 2]
        p = 3  # Start from 3 to ensure an odd number
        while p < L**(1/np.sqrt(2)):
            if prime[p] and sympy.legendre_symbol(self.n, p) == 1:
                primes.append(p)
            p += 2  # Increment by 2 to consider only odd numbers
        return primes
    
    def __generate_continued_fraction__(self, k):
        v = [1]
        alpha = [math.sqrt(self.n)]
        a = [int(math.sqrt(self.n))]
        u = [int(math.sqrt(self.n))]

        for i in range(k):
            v.append((self.n-u[i]**2) / v[i])
            alpha.append((math.sqrt(self.n)+u[i]) / v[i+1])
            a.append(int(alpha[i+1]))
            u.append(a[i+1]*v[i+1]-u[i])
        
        return a[0:]
    
    def __generate_smooth_numbers__(self, a):
        b = []
        for i in range(len(a)):
            if i == 0:
                b.append((a[i]*1) % self.n)
            elif i == 1:
                b.append((a[i]*b[0]+1) % self.n)
            else:
                b.append((a[i]*b[i-1]+b[i-2]) % self.n)

        return b
    
    def __matrix_reduction__(self, matrix):
        det = []
        nz = find_index_of_nonzero_array(matrix)
        m, n = matrix.shape
        for j in range(n):
            if len(det) == len(nz)-1:
                break
            i = np.where(matrix[:, j] == 1)[0]
            if i.size != 0: 
                i = i[0]
                det.append(i)
                rang = list(range(n))
                rang.remove(j)
                for k in rang:
                    if matrix[i, k] == 1:
                        matrix[:, k] = xor(matrix[:, j], matrix[:, k])
        ud = set(range(m)).difference(set(det))
        undet = []
        for i in ud:
            undet.append(i)

        return (matrix, sorted(det), sorted(undet))

    def __solve_slr__(self, matrix):
        matrix, det, undet = self.__matrix_reduction__(matrix)
        m, n = matrix.shape
        '''full_solution = []
        for x in undet:
            solution = np.zeros(n)
            for y in det:
                if np.where((matrix[x] == 1) & (matrix[y] == 1))[0].size != 0:
                    matrix[x] = xor(matrix[x], matrix[y])
                    solution = xor(solution, matrix[x])
            if (np.nonzero(solution)[0].size != 0):
                full_solution.append(solution)'''
        undet_nz = find_index_of_nonzero_array(matrix[undet])[0]
        solution = []
        for i in range(m):
            if i == undet[undet_nz]:
                solution.append(1)
            else:
                if check_matching_positions([matrix[i], matrix[undet[undet_nz]]]):
                    solution.append(1)
                else:
                    solution.append(0)

        for i in range(n-len(solution)):
            solution.append(0)

        return solution
    
    def __find_mod_prime_degrees__(self, base, num):
        x1 =  self.__find_prime_degrees__(base, num)
        x2 =  self.__find_prime_degrees__(base, num-self.n)

        if is_array_nonzero(np.array(x1)) and x1 != None:
            return np.array(x1)
        elif is_array_nonzero(np.array(x2)) and x2 != None:
            return np.array(x2)
        else:
            return None
    
    def __find_prime_degrees__(self, base, num):

        d = [0] * len(base)
        if -1 in base:
            if num < 0:
                d[0] = 1
                num = abs(num)
            
        for i in range(len(base)):
            while num % base[i] == 0 and base[i] != -1:
                num /= base[i]
                d[i] += 1
        if num == 1:
            return d
        else:
            return None

F = Factorizer(17873, method="brillhart-morrison", patience=50)
print(F.factorize())
