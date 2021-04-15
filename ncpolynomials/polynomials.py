"""
Defining the classes of Operators and Monomials

Operators can be multiplied to form Monomials

Monomials can be summed to form polynomials...

Author: Peter J. Brown (02/12/20)
"""
from itertools import product
from copy import deepcopy
from . import simplification_utils as tools


class Operator(object):
    """
    Operator Class

    Description:
                    Object representing an abstract operator in some *-algebra.

    Attributes:
                _name        String labelling the operator
                _hermitian   Bool indicating whether the operator is Hermitian
                _adjoint     Bool indicating whether the operator is the adjoint
    """

    # Class constructor
    def __init__(self, name = 'Id', hermitian = True, adjoint = False):
        if isinstance(name, Operator):
            self._name = str(name.name)
            self._hermitian = name.hermitian
            self._adjoint = name.adjoint
        else:
            self._name = name
            self._hermitian = hermitian
            self._adjoint = adjoint

    def __hash__(self):
        return hash((self.name, self.adjoint))

    def __eq__(self, other):
        if (self.name == other.name) and (self.adjoint == other.adjoint):
            return True
        else:
            return False

    def __mul__(self, other):
        return Monomial(self) * Monomial(other)
    def __rmul__(self, other):
        if isinstance(other, (int,float,complex)):
            return Monomial([self], 1) * other
    def __div__(self, other):
        if isinstance(other, (int,float,complex)):
            return (1/other) * self
    def __add__(self, other):
        return Polynomial(self) + Polynomial(other)
    def __sub__(self, other):
        return Polynomial(self) + (-Polynomial(other))
    def __radd__(self, other):
        return self + other
    def __neg__(self):
        return Monomial([self], -1)


    # Getters/Setters
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self, val):
        if isinstance(val, str):
            self._name = val
        else:
            raise TypeError('Attribute \'name\' should be a string')

    @property
    def hermitian(self):
        return self._hermitian
    @property
    def adjoint(self):
        return self._adjoint
    @adjoint.setter
    def adjoint(self, val):
        if isinstance(val, bool):
            self._adjoint = val
        else:
            raise TypeError('Attribute \'adjoint\' should be a boolean value')


    def adj(self):
        X = Operator(self)
        if X.hermitian == False:
            X.adjoint = not X.adjoint
        return X

    def simplify(self):
        # If we have taken the adjoint of a hermitian operator then remove it
        if (self.hermitian == 1) and (self.adjoint == 1):
            self.adjoint = False
        return self

    def __repr__(self):
        string = self.name
        if self.adjoint:
            string += '\''
        return string


class Monomial(object):
    """
    Monomial Class

    Description:
                    Object representing a product of operators

    Attributes:
                _terms      list of operators in the product
                _coef       coefficient of the monomial
    """

    # Class constructor
    def __init__(self, terms = [], coef = 1):
        if isinstance(terms, Operator):
            self.terms = [deepcopy(terms)]
            self.coef = 1
        elif isinstance(terms, Monomial):
            self.terms = deepcopy(terms.terms)
            self.coef = terms.coef
        elif isinstance(terms, (int, float, complex)):
            self.terms = []
            self.coef = terms
        else:
            self.terms = deepcopy(terms)
            self.coef = coef

    def __len__(self):
        return len(self.terms)

    def __hash__(self):
        x = tuple(self.terms + [self.coef])
        return hash(x)

    def __eq__(self, other):
        equal = True
        if (len(self) != len(other)) or (self.coef != other.coef):
            equal = False
        else:
            for i in range(len(self.terms)):
                if (self.terms[i].name != other.terms[i].name) or (self.terms[i].adjoint != other.terms[i].adjoint):
                    equal = False
                    break
        return equal

    def __mul__(self, other):
        if isinstance(other, (int,float, complex, Operator)):
            return self * Monomial(other)
        elif isinstance(other, Monomial):
            terms = self.terms + other.terms
            coef = self.coef * other.coef
            return Monomial(terms, coef)
        elif isinstance(other, Polynomial):
            return Polynomial(self) * other
        else:
            raise TypeError('Bad type for multiplication with Monomial')
    def __rmul__(self, other):
        if isinstance(other, (int,float,complex)):
            return self * other
    def __div__(self, other):
        if isinstance(other, (int,float,complex)):
            return (1/other) * self
    def __add__(self, other):
        return Polynomial(self) + Polynomial(other)
    def __radd__(self, other):
        return self + other
    def __sub__(self, other):
        return Polynomial(self) + (-Polynomial(other))
    def __rsub__(self, other):
        return -(self - other)

    def __neg__(self):
        return Monomial(self.terms, -self.coef)

    def __mod__(self, other):
        """
        Overloads the % comparison operator.
        We will use this to define equality ignoring the coefficients
        """
        equal = True
        other = Monomial(other)
        if len(self) != len(other):
            equal = False
        else:
            for i in range(len(self.terms)):
                if (self.terms[i].name != other.terms[i].name) or (self.terms[i].adjoint != other.terms[i].adjoint):
                    equal = False
                    break
        return equal

    # Getters/Setters
    @property
    def terms(self):
        return self._terms
    @terms.setter
    def terms(self, val):
        if isinstance(val, list):
            self._terms = val
        else:
            raise TypeError('Attribute \'terms\' should be a list of operators')

    @property
    def coef(self):
        return self._coef
    @coef.setter
    def coef(self, val):
        if isinstance(val, (int, float, complex)):
            self._coef = val
        else:
            raise TypeError('Attribute \'coef\' should be a number value')

    @property
    def degree(self):
        return len(self)

    @property
    def hermitian(self):
        X = self.adj()
        return self % X

    def adj(self):
        X = Monomial(self)
        # Reverse the order
        X.terms.reverse()
        # Apply adjoints to the operators in list in not hermitian
        X.terms = [op.adj() if op.hermitian == False else op for op in X.terms]
        # If the coefficient is complex then we should conjugate
        if isinstance(X.coef, complex):
            X.coef = X.coef.conjugate()
        return X

    def apply_substitution(self, old, new):
        """
        Tries to find Monomial old in Monomial self.
        If found replaces with Monomial new
        """
        old_term = Monomial(old)
        new_term = Monomial(new)
        success = False
        if len(old_term) <= len(self):
            for m in range(len(self) - len(old_term) + 1):
                if Monomial(self.terms[m : m + len(old_term)], 1) % Monomial(old_term.terms, 1):
                    success = True
                    left_terms = self.terms[0 : m]
                    right_terms = self.terms[m + len(old_term) : len(self)]
                    self.terms = left_terms + new_term.terms + right_terms

                    # if the coefficient of old_term was not 1 then we should divide through
                    self.coef = self.coef * new_term.coef / old_term.coef
                    break

        return success

    def simplify(self, subs):
        """
        Takes a dictionary of substitution rules and applies the substitutions
        until no more simplifications are found.
        """
        success = False

        for old, new in subs.items():
            success = self.apply_substitution(old, new)
            # If we do make a sub then we should try the previous subs again
            # until we don't find one anymore
            if success == True:
                self.simplify(subs)
                break

        # If the coefficient is 0 then remove the terms
        if self.coef == 0:
            self.terms = []
        return Monomial(self)



    def __repr__(self):
        string = ''
        if len(self.terms) > 0:
            for term in self.terms:
                string += term.__repr__()
        else:
            string += 'Id'
        return str(self.coef) + '*' + string


class Polynomial(object):
    """
    Polynomial Class

    Description:
                    Object representing sums of monomials

    Attributes:
                _terms      list of monomials in the sum
    """

    # Class constructor
    def __init__(self, terms = []):
        if isinstance(terms, (int, float, complex, Operator, Monomial)):
            self.terms = [Monomial(terms)]
        elif isinstance(terms, Polynomial):
            self.terms = deepcopy(terms.terms)
        else:
            self.terms = deepcopy(terms)

    def __len__(self):
        return len(self.terms)

    def __hash__(self):
        x = tuple(self.terms)
        return hash(x)

    def __eq__(self, other):
        equal = True
        other = Polynomial(other)
        if (len(self) != len(other)):
            equal = False
        else:
            for i in range(len(self.terms)):
                if Monomial(self.terms[i]) != Monomial(other.terms[i]):
                    equal = False
                    break
        return equal

    def __add__(self, other):
        if isinstance(other, (int, float, complex, Operator, Monomial)):
            return self + Polynomial(other)
        elif isinstance(other, Polynomial):
            return Polynomial(self.terms + other.terms)
        else:
            raise TypeError('Bad type for addition with Polynomial')
    def __radd__(self, other):
        # Addition is commutative
        return self + other
    def __sub__(self, other):
        return Polynomial(self) + (-Polynomial(other))
    def __rsub__(self, other):
        return -(self - other)

    def __mul__(self, other):
        if isinstance(other, (int, float, complex, Operator, Monomial)):
            return self * Polynomial(other)
        elif isinstance(other, Polynomial):
            newterms = [term1 * term2 for term1, term2 in product(self.terms, other.terms)]
            return Polynomial(newterms)
        else:
            raise TypeError('Bad type for multiplication with Polynomial')
    def __rmul__(self, other):
        # If multiplication is with number then it is reflexive
        if isinstance(other, (int,float,complex)):
            return self * other

    def __div__(self, other):
        if isinstance(other, (int,float,complex)):
            return (1/other) * self

    def __neg__(self):
        return Polynomial([-term for term in self.terms])

    # Getters/Setters
    @property
    def terms(self):
        return self._terms
    @terms.setter
    def terms(self, val):
        if isinstance(val, list):
            self._terms = val
        else:
            raise TypeError('Attribute \'terms\' should be a list of Monomials')

    @property
    def degree(self):
        return max([len(term) for term in self.terms])

    def adj(self):
        return Polynomial([term.adj() for term in self.terms])

    def simplify(self, subs):
        """
        Simplifies a polynomial using the substitutions subs
        """
        success = False
        # First simplify each term in the poly individually
        self.terms = [term.simplify(subs) for term in self.terms[:]]

        # Now we collect like terms
        umonos, idx = tools.unique_monomials(self.terms)
        self.terms = [Monomial(umonos[i].terms, sum([self.terms[j].coef for j in idx[i]])) for i in range(len(idx))]
        # After simplifying we may have introduced some zero terms
        self.terms = [term for term in self.terms[:] if term.coef != 0]

        return Polynomial(self)

    def iszero(self):
        # Checks if a polynomial is zero
        # Warning it is advisable to simplify the polynomial first
        coefs = [term.coef for term in self.terms if abs(term.coef) > 1e-10]
        if len(coefs) > 0:
            return False
        else:
            return True



    def __repr__(self):
        string = ''
        first = True
        if len(self.terms) > 0:
            for term in self.terms:
                if not first:
                    string += ' + ' + term.__repr__()
                else:
                    string += term.__repr__()
                    first = False
        else:
            string += '0'
        return string
