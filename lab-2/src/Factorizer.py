import time

def find_index(arr, element):
    for i in range(len(arr)):
        if arr[i] == element:
            return i
    return -1  # Element not found

def array_product(arr):
    product = 1
    for num in arr:
        product *= num
    return product

class Factorizer:
    def __init__(self, a, b, p, method="brute-force"):
        self.p = int(p)
        self.a = int(a)
        self.b = int(b)
        self.method = method

    def factorize(self):
        if self.method == "brute-force":
            return self.__brute_force__()
        elif self.method == "sph":
            return self.__silver_pohlig_hellman__()

        
    def __brute_force__(self):
        x = 1
        start_time = time.time()

        while x < self.p and  time.time() - start_time < 300:
            tmp = pow(self.a, x, self.p)

            if tmp == self.b:
                return x
            
            x += 1

    def __silver_pohlig_hellman__(self):
        self.n = self.p - 1
        canonical_n = self.__canonical_form__(self.n)
        r_table = self.__form_table__(canonical_n)
        y = []

        for pi, li in canonical_n.items():
            tmp = pow(self.b, self.n//pi, self.p)
            x = find_index(r_table[pi], tmp)
            for k in range(0, li):
                tmp = pow(self.b * pow(pow(self.a, -1, self.p), x, self.p), self.n // (pi**(k+1)), self.p)
                x += find_index(r_table[pi], tmp) * pow(pi, k) % self.p
            y.append(x)

        P = [pow(pi, li) for pi, li in canonical_n.items()]

        return self.__crt__(y, P)


    def __canonical_form__(self, n):
        result = {}
        divisor = 2
        while n > 1:
            if n % divisor == 0:
                if divisor not in result:
                    result[divisor] = 1
                else:
                    result[divisor] += 1
                n //= divisor
            else:
                divisor += 1
        return result
    
    def __form_table__(self, canonical_n):
        table = {}
        for pi in canonical_n.keys():
            table[pi] = []
            for j in range(0, pi):
                table[pi].append(pow(self.a, self.n*j//pi, self.p))

        return table
    
    def __crt__(self, a, n):
        N = array_product(n)
        y = [N// i for i in n]
        z = [pow(y[i], -1, n[i]) for i in range(len(n))]

        x = sum([a[i]*y[i]*z[i] for i in range(len(n))]) % N
        return x
    


a = 2181
b = 4539
p = 4733
                
F = Factorizer(a, b, p, method="sph")
print(F.factorize())