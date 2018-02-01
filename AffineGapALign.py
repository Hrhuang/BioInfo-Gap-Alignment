from pprint import pprint
from itertools import permutations
import sys
from pprint import pprint


global match
global mismatch
global Gap
global StartingGap
global MIN


match = 6.
mismatch = -2.
StartingGap = -10.
Gap = -1
MIN = -float("inf")

def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


#return match or mismatch score
def _match(s, t, i, j):
    if t[i-1] == s[j-1]:
        return match
    else:
        return mismatch

#initializers for matrices
def _init_x(i, j):
    if i > 0 and j == 0:
        return MIN
    else:
        if j > 0:
            return -10 + (-0.5 * j)
        else:
            return 0

def _init_y(i, j):
    if j > 0 and i == 0:
        return MIN
    else:
        if i > 0:
            return -10 + (-0.5 * i)
        else:
            return 0

def _init_m(i, j):
    if j == 0 and i == 0:
        return 0
    else:
        if j == 0 or i == 0:
            return MIN
        else:
            return 0

def _format_tuple(inlist, i, j):
    return 0

def distance_matrix(s, t):
    dim_i = len(t) + 1
    dim_j = len(s) + 1
    #score of best alignment of x[1..i] and y[1..j] ending with a space in X.
    X = [[_init_x(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
    #score of best alignment of x[1..i] and y[1..j] ending with a space in Y.
    Y = [[_init_y(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]
    #score of best alignment of x[1..i] and y[1..j] ending with a charactercharacter match or mismatch.
    M = [[_init_m(i, j) for j in range(0, dim_j)] for i in range(0, dim_i)]

    for j in range(1, dim_j):
        for i in range(1, dim_i):
            X[i][j] = max((StartingGap + Gap + M[i][j-1]), (Gap + X[i][j-1]), (StartingGap + Gap + Y[i][j-1]))
            Y[i][j] = max((StartingGap + Gap + M[i-1][j]), (StartingGap + Gap + X[i-1][j]), (Gap + Y[i-1][j]))
            M[i][j] = max(_match(s, t, i, j) + M[i-1][j-1], X[i][j], Y[i][j])

    return [X, Y, M]

def printMat(A):    
    
    print '['
    for row in A:
       print row
    print ']'

def backtrace(s, t, X, Y, M):
    sequ1 = ''
    sequ2 = ''
    i = len(t)
    j = len(s)

    print("\nThe optimal alignment score: ")        
    print(max(M[i][j], X[i][j], Y[i][j])) 
    print ("\n") 

    while (i>0 or j>0):
        if (i>0 and j>0 and M[i][j] == M[i-1][j-1] + _match(s, t, i, j)):
            sequ1 += s[j-1]
            sequ2 += t[i-1]
            i -= 1; j -= 1
        elif (i>0 and M[i][j] == Y[i][j]):
            sequ1 += '_'
            sequ2 += t[i-1]
            i -= 1
        elif (j>0 and M[i][j] == X[i][j]):
            sequ1 += s[j-1]
            sequ2 += '_'
            j -= 1

      

    sequ1r = ' '.join([sequ1[j] for j in range(-1, -(len(sequ1)+1), -1)])
    sequ2r = ' '.join([sequ2[j] for j in range(-1, -(len(sequ2)+1), -1)])

    return [sequ1r, sequ2r]




if len(sys.argv) == 2:
    filename = sys.argv[1]
else:
    sys.exit("Usage: " + sys.argv[0] + " [fasta_filename]")

reads   = []
names   = []

with open(filename) as fp:
    for (name, seq) in read_fasta(fp):
        names.append(name)
        reads.append(seq)

perms = permutations(reads, 2)

for perm in perms:
  
    print("Sequence1 Name: " + names[reads.index(perm[0])])
    print("Sequence: " + perm[0])
  
    print("Sequence2 Name: " + names[reads.index(perm[1])])
    print("Sequence: " + perm[1])
    print("\n")
    [X, Y, M] = distance_matrix(perm[0], perm[1])

    #Printing the matrixes
    print("The matrix for Ix(I,j):\n")
    printMat(X)

    print("\n")

    print("The matrix for Iy(I,j):\n")
    printMat(Y)
    print("\n")

    print("The matrix for F(I,j):\n")
    printMat(M)
    print("\n")

    [str1, str2] = backtrace(perm[0], perm[1], X, Y, M)
    print("The optimal alignment:")
    print(str1)
    print(str2)
    print("\n")
