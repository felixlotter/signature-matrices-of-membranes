import numpy as np
import math
import timeit


def row_sum(matrix):
    return np.sum(matrix, axis=1)

def column_sum(matrix):
    return np.sum(matrix, axis=0)

def weighted_shift(matrix): # this is (the truncation of) integration of polynomials given in matrix form
    m, n = matrix.shape
    out = np.zeros((m, n))
    for i in range(1,m):
        for j in range(1,n):
            coef = 1/(math.factorial(i) * math.factorial(j))
            out[i, j] = coef * (matrix[i-1, j-1])
    return out

def sig(matrixSeq, word):
    d = len(matrixSeq)
    m, n = matrixSeq[0].shape
    k = len(word)
    #initialize sheet in O(m n k^2)
    sheet = np.zeros((m,n,k+1,k+1))
    for i,j in np.ndindex(m,n):
            sheet[i][j][0,0]+=1
    for p in range(k):
        letter = word[p] - 1
        #multiply with new coefs and integrate in O(m n k^2)
        for i, j in np.ndindex(m,n):
            sheet[i,j] *= matrixSeq[letter][i,j]
            sheet[i,j] = weighted_shift(sheet[i,j])

        #"marginal" accumulated sums in O(m n k)
        for i in range(m-1):
                rsum = [row_sum(sheet[i,j]) for j in range(n)]
                for j in range(n):
                    sheet[i+1,j,:,0] += rsum[j]
        for j in range(n-1):
                csum = [column_sum(sheet[i,j]) for i in range(m)]
                for i in range(m):
                    sheet[i,j+1,0,:] += csum[i]
    #finally, the result is obtained as the total sum of the bottom right matrix on the sheet in O(m n)
    out = row_sum([column_sum(sheet[m-1,n-1])])[0]
    return(out)

membrane = [np.array([[4,2],[1,3]]), np.array([[2,1],[0,-4]])] #creates a membrane in 2-dimensional space, given as a sequence of matrices (each matrix corresponding to increments in one coordinate)
word = [2,2]
sig(membrane,word)

# Random test input

m = 100
n = 100
d = 2
k = 5

membrane = [np.random.rand(m,n)*2-1 for _ in range(d)] #creates a membrane in k-dimensional space, given as a sequence of matrices (each matrix corresponding to increments in one coordinate)

word = list(np.random.randint(1, d+1, k)) #creates a random word of length k, with letters in {1,...,d}
wordstr = ''.join(str(x) for x in word)

start_time = timeit.default_timer()
out = sig(membrane, word)
elapsed_time = timeit.default_timer() - start_time

# Print the output and elapsed time
print(f"Evaluated signature of bilinear membrane of order {(m,n)} at word {wordstr} in {elapsed_time:.6f} seconds. Result: {out}")

