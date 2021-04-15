"""
Some functions related to the simplification of nc polys
"""
from . import polynomials as poly


def flatten(lol):
    """Flatten a list of lists to a list.
    :param lol: A list of lists in arbitrary depth.
    :type lol: list of list.
    :returns: flat list of elements.
    """
    new_list = []
    for element in lol:
        if element is None:
            continue
        elif not isinstance(element, list) and not isinstance(element, tuple):
            new_list.append(element)
        elif len(element) > 0:
            new_list.extend(flatten(element))
    return new_list

def num2base(n, b, d):
    if n == 0:
        digits = [0 for _ in range(d)]
    else:
        digits = []
        while n:
            digits.append(int(n % b))
            n //= b
        digits = digits[::-1]
        if len(digits) < d:
            digits = [0 for _ in range(d-len(digits))] + digits
    return digits

def gen_indexing(b,d):
    # Generates all base b strings of length d
    indexes = [num2base(n,b,d) for n in range(b**d)]
    return indexes


def generate_operators(name, n_vars=1, hermitian=1):
    """
    Generates a list of operators with names name0, name1, name2,...
    """
    ops = [poly.Monomial([poly.Operator(name + str(k), hermitian)], 1) for k in range(n_vars)]
    return ops

def get_monomials(base_monomials, degree):
    """
    Grabs all monomials of degree 1 from the base set and then generates all
    products up to degree d.
    """
    id = poly.Monomial([])
    monoset, _ = unique_monomials(pick_monomials_of_degree(base_monomials, 1))
    monos = [id]

    # Generate all base-len(monoset) strings for length up to d
    for d in range(1, degree+1):
        indexes = gen_indexing(len(monoset), d)
        for index in indexes:
            temp = id
            for k in index:
                temp = temp * monoset[k]
            monos.append(temp)
    return monos

def get_all_unique_monomials(base_monomials, degree = 1, subs = {}, extra_monomials = []):
    """
    Generates all monomials up to some degree using the base_monomials set.
    Then adding the extra_monomials it simplifies all monomials and picks the
    remaining unique ones out.
    """
    if degree == 0:
        return [poly.Monomial([])]
    monos = get_monomials([poly.Monomial([])] + base_monomials, degree) + extra_monomials
    monos = [mon.simplify(subs) for mon in monos]
    monos, _ = unique_monomials(monos)
    return monos


def pick_monomials_of_degree(mono_list, degree):
    """
    Returns all monomials in mono_list that have degree=degree
    """
    monolist = [mono for mono in mono_list if mono.degree == degree]

    return monolist



def unique_monomials(mono_list):
    """
    Given a list of monomials it returns a list of unique monomials (umonos),
    and then a list of their indexes in the original list

    NOTE: this function does not care about the coefficients of the monomials.
    """
    umonos = []
    idx = []

    for mono_ind, mono in enumerate(mono_list):
        # Check if mono is same as any of the unique monos
        found = False
        for umono_ind, umono in enumerate(umonos):
            if umono % mono:
                # If mono already in umono
                found = True
                idx[umono_ind].append(mono_ind)

        if not found:
            umonos.append(mono)
            idx.append([mono_ind])

    return umonos, idx
