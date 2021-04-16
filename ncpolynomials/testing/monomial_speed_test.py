import sys
sys.path.append('/home/pbrown/Dropbox/python/')
from ncpolynomials.simplification_utils import get_all_unique_monomials, get_monomials, generate_operators
import ncpol2sdpa as ncp
import timeit

test1 = """
from ncpolynomials.simplification_utils import generate_operators, get_monomials
n = 2
d = 6
A = generate_operators('A', n, 1)
subs1 = {}
for i in range(n):
    subs1.update({A[i] * A[i] : A[i]})
X = get_monomials(A, d)
"""
test2 = """
import ncpol2sdpa as ncp
n = 2
d = 6
B = ncp.generate_operators('B', n, 1)
subs2 = {}
for i in range(n):
    subs2.update({B[i] * B[i] : B[i]})
Y = ncp.nc_utils.get_monomials(B, d)
"""

time1 = timeit.timeit(test1, number=50)/50
time2 = timeit.timeit(test2, number=50)/50
print(time1, time2)
