import sys
sys.path.append('/home/pbrown/Dropbox/python/')
from ncpolynomials.simplification_utils import generate_operators
from ncpolynomials.polynomials import Monomial
id = Monomial([])
X = generate_operators('X', 2, 1)
Y = generate_operators('Y', 2, 1)
Z = generate_operators('Z', 2, 1)

Aops = [X[0], Y[0], Z[0]]
Bops = [X[1], Y[1], Z[1]]

subs = {}
for x in range(3):
    for y in range(3):
        subs.update({Bops[y] * Aops[x] : Aops[x] * Bops[y]})

    subs.update({Aops[x]*Aops[x] : id})
    subs.update({Bops[x]*Bops[x] : id})

poly1 = (id - Z[0] + Z[1] - Z[1]*Z[0])
poly1 = poly1 * X[0]
poly1 = poly1 * (id + X[0]*Y[0] + X[1]*Y[1] + X[1]*Y[1]*X[0]*Y[0])
poly1 = poly1 * (id + Y[0]*X[0] + Y[1]*X[1] + Y[0]*X[0]*Y[1]*X[1])
poly1 = poly1 * X[1]
poly1 = poly1 * (id + Z[0] - Z[1] - Z[0]*Z[1])
poly1 = poly1.simplify(subs)

poly2 = id + Z[0] - Z[1] - Z[0]*Z[1]
poly2 = poly2 * X[1]
poly2 = poly2 * (id + X[0]*Y[0] + X[1]*Y[1] + X[1]*Y[1]*X[0]*Y[0])
poly2 = poly2 * (id + Y[0]*X[0] + Y[1]*X[1] + Y[0]*X[0]*Y[1]*X[1])
poly2 = poly2 * X[0]
poly2 = poly2 * (id - Z[0] + Z[1] - Z[0]*Z[1])
poly2 = poly2.simplify(subs)

P = poly1 - poly2
Q = P.simplify(subs)
print(len(Q))
print(Q)
