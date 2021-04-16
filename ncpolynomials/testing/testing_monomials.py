import sys
sys.path.append('/home/pbrown/Dropbox/python/')
from ncpolynomials.simplification_utils import *
from ncpolynomials.polynomials import *

from copy import copy, deepcopy
A = generate_operators('A',3,0)[0]
print(type(A.terms[0]))
B = Operator(A.terms[0])
print(A*A == A*A)

y = A*A + A
subs = {A*A : A}
print(y.simplify(subs))
