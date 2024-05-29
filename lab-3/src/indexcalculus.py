import Factorizer
import numpy as np
import threading
import math

def is_subset(list1, list2):
    set1 = set(list1)
    set2 = set(list2)
    return set1 <= set2

def chech_linear_dependency(A, b):
    if A == []:
        return True
    matrix = np.row_stack((A, b))
    rank = np.linalg.matrix_rank(np.array(matrix))
    
    if rank == len(matrix):
        return True
    else:
        return False
    
def columns_with_only_zeros(matrix):
    cols_with_only_zeros = np.all(matrix == 0, axis=0)
    indices = np.where(cols_with_only_zeros)[0]
    
    return indices.tolist()

class IC:
    def __init__(self, a, b, n):
        self.n = int(n)
        self.a = int(a)
        self.b = int(b)
        self.lock = threading.Lock()
        self.stop_event = threading.Event()

    def solve(self):
        self.base = self.__generate_factor_base__()
        sle = self.__generate_sle__()   
        solution = self.__solve_sle__()  
        l, factors = self.__generate_smooth_bal__()
        x = 0

        for i in range(len(self.base)):
            x += factors[i] * int(solution[i]) % (self.n - 1)

        x = (x - l) % (self.n - 1)

        return int(x)

    def __generate_factor_base__(self):
        c = 3.38
        B = int(c * np.exp(0.5 * np.sqrt(np.log2(self.n) * np.log2(np.log2(self.n)))))
        #B = 6

        prime = [True for _ in range(B)]
        p = 2

        while (p * p < B):
            if prime[p]:
                for i in range(p * p, B, p):
                    prime[i] = False
            p += 1

        prime_numbers = [p for p in range(2, B) if prime[p]]
        return prime_numbers
    
    def __smooth_check__(self, k):
        ak = pow(self.a, k, self.n)
        '''factors = Factorizer.main(ak)

        if factors != None:
            if is_subset(factors, self.base):
                eq = np.zeros(len(self.base) + 1)

                for i in range(len(self.base)):
                    eq[i] = factors.count(self.base[i])

                eq[-1] = k

                return eq, True
            else:
                return [], False
        else:
            return [], False'''
        exponents = [0] * (len(self.base) + 1)
        
        for i, prime in enumerate(self.base):
            while ak % prime == 0:
                exponents[i] += 1
                ak //= prime
        
        if ak == 1:
            exponents[-1] = k
            return exponents, True
        else:
            return [], False
        
    def __smooth_check_bal__(self, l):
        ball = self.b * pow(self.a, l, self.n) % self.n
        factors = Factorizer.main(ball)

        if factors != None:
            if is_subset(factors, self.base):
                eq = np.zeros(len(self.base) + 1)

                for i in range(len(self.base)):
                    eq[i] = factors.count(self.base[i])

                eq[-1] = l

                return eq, True
            else:
                return [], False
        else:
            return [], False
        
    def worker(self, number, sle, stop_event):
        if stop_event.is_set():
            return
        number, is_smooth_number = self.__smooth_check__(number)
        if is_smooth_number:
            sle.append(number)
            if len(sle) >= len(self.base) + 15:
                stop_event.set()

    def __worker__(self, number_range):
        for number in number_range:
            if self.stop_event.is_set():
                break
            number, is_smooth_number  = self.__smooth_check__(number)
            if is_smooth_number:
                with self.lock:
                    self.sle.append(number)
                    if len(self.sle) >= len(self.base) + 15:
                        self.stop_event.set()
                        break
        
        
    def __generate_sle__(self):
        self.sle = []

        # Calculate chunk size
        chunk_size = max(1, self.n // (4 * (len(self.base) + 15)))
        print(f"Chunk size: {chunk_size}")

        # Create tasks as ranges of numbers
        tasks = [range(i, min(i + chunk_size, self.n)) for i in range(0, self.n, chunk_size)]

        # Create and start threads
        threads = []
        for task in tasks:
            thread = threading.Thread(target=self.__worker__, args=(task,))
            thread.start()
            threads.append(thread)

        # Wait for all threads to complete
        for thread in threads:
            thread.join()

        self.sle = np.array(self.sle)

        # Return the first len(self.base) + 15 smooth numbers found
        return np.array(self.sle[:len(self.base) + 15])




        '''sle = []
        i = 0

        while len(sle) <= len(self.base)+35:
            number, is_smooth_number = self.__smooth_check__(i)
            if is_smooth_number:
                sle.append(number)

            i += 1'''
        return np.array(sle, dtype=int)
        
    
    def __solve_sle__(self):
        '''
        GF = galois.GF(self.n)

        A = GF(sle[:, :-1])
        b = GF(sle[:, -1])

        x = np.linalg.solve(A, b) 

        return x'''
        
        n, m = self.sle.shape

        processed = []

        # Gauss elimination
        for j in range(m - 1):
            # Calculate GCD of each element in the column with mod
            gcd_col_mod = np.fromiter((math.gcd(elem, (self.n - 1)) for elem in self.sle[:, j]), dtype=int)
            
            for i, gcd_val in enumerate(gcd_col_mod):
                if i in processed:
                    continue
                
                if gcd_val == 1:
                    inv_elem = pow(int(self.sle[i, j]), -1, (self.n - 1))
                    processed.append(i)

                    self.sle[i] = (self.sle[i] * inv_elem) % (self.n - 1)
                    
                    # Vectorized row elimination
                    mask = np.arange(n) != i
                    self.sle[mask] = (self.sle[mask] - np.outer(self.sle[mask, j], self.sle[i])) % (self.n - 1)
                    break

        solution = []
        for j in range(m - 1):
            non_zero_indices = np.nonzero(self.sle[:, j])[0]  # Find indices of non-zero elements in the j-th column
            if non_zero_indices.size > 0:
                solution.append(self.sle[non_zero_indices[0], -1])  # Append the last non-zero element in the j-th column
            else:
                solution.append(0)  # If all elements are zero in the j-th column, append 0
    
        return solution
    
    def __generate_smooth_bal__(self):
        i = 0
        while True:
            number, is_smooth_number = self.__smooth_check_bal__(i)
            if is_smooth_number:
                return number[-1], number[:-1]
            else:
                i += 1

if __name__ == '__main__':
    ic = IC(14837213830, 38662976351, 44392159481)
    print(ic.solve())



