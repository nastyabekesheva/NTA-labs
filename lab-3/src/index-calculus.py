import Factorizer
import numpy as np
from multiprocessing import Pool, Manager
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
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

    def solve(self):
        self.base = self.__generate_factor_base__()
        sle = self.__generate_sle__()   
        print(sle)
        solution = self.__solve_sle__(sle)  
        l, factors = self.__generate_smooth_ball__()
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
        factors = Factorizer.main(ak)

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
            return [], False
        
    def __smooth_check_ball__(self, l):
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
            if len(sle) >= len(self.base) +10:
                stop_event.set()
        
        
    def __generate_sle__(self):
        '''with Pool(processes=4) as pool:
            sle = []
            # Use imap_unordered to process the numbers in parallel
            for number, is_smooth_number in pool.map(self.__smooth_check__, range(self.n)):
                if is_smooth_number:
                    #if chech_linear_dependency(sle, number):
                    sle.append(number)
                if len(sle) >= len(self.base) + 15:
                    pool.terminate()  # Terminate the pool if we have enough results
                    break '''  
        i = 0
        sle = []

        num_to_find = len(self.base) + 15

        with Manager() as manager:
            sle = manager.list()  # Shared list among processes
            stop_event = manager.Event()  # Event to signal when to stop

            with Pool() as pool:
                tasks = [(i, sle, stop_event) for i in range(self.n)]  # Arbitrarily large range

                pool.starmap(self.worker, tasks)
                # Convert the manager list to a regular list
                sle = list(sle)[:len(self.base) + 15]




        '''sle = []
        i = 0

        while len(sle) <= len(self.base)+15:
            number, is_smooth_number = self.__smooth_check__(i)
            if is_smooth_number:
                sle.append(number)

            i += 1

        print(len(sle))'''
        return np.array(sle, dtype=int)
        
    
    def __solve_sle__(self, sle):
        '''
        GF = galois.GF(self.n)

        A = GF(sle[:, :-1])
        b = GF(sle[:, -1])

        x = np.linalg.solve(A, b) 

        return x'''
        
        n, m = sle.shape

        processed = []

        # Gauss elimination
        for j in range(m - 1):
            # Calculate GCD of each element in the column with mod
            gcd_col_mod = np.fromiter((math.gcd(elem, (self.n - 1)) for elem in sle[:, j]), dtype=int)
            
            for i, gcd_val in enumerate(gcd_col_mod):
                if i in processed:
                    continue
                
                if gcd_val == 1:
                    inv_elem = pow(int(sle[i, j]), -1, (self.n - 1))
                    processed.append(i)

                    sle[i] = (sle[i] * inv_elem) % (self.n - 1)
                    
                    # Vectorized row elimination
                    mask = np.arange(n) != i
                    sle[mask] = (sle[mask] - np.outer(sle[mask, j], sle[i])) % (self.n - 1)
                    break

        solution = []
        for j in range(m - 1):
            non_zero_indices = np.nonzero(sle[:, j])[0]  # Find indices of non-zero elements in the j-th column
            if non_zero_indices.size > 0:
                solution.append(sle[non_zero_indices[0], -1])  # Append the last non-zero element in the j-th column
            else:
                solution.append(0)  # If all elements are zero in the j-th column, append 0
    
        return solution
    
    def __generate_smooth_ball__(self):
        i = 0
        while True:
            number, is_smooth_number = self.__smooth_check_ball__(i)
            if is_smooth_number:
                return number[-1], number[:-1]
            else:
                i += 1
        '''with Pool(processes=4) as pool:
            for i, (number, is_smooth_number) in enumerate(pool.imap_unordered(self.__smooth_check_ball__, range(self.n))):
                if is_smooth_number:
                    pool.terminate()
                    return number[-1], number[:-1]'''
                








if __name__ == '__main__':
    start_time = time.time()

    ic = IC(409634, 294022, 294022)
    #ic = IC(10, 17, 47)
    print(ic.solve())

    end_time = time.time()
    execution_time = end_time - start_time
    print("Execution time:", execution_time, "seconds")


