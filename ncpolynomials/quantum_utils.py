"""
Some functions related to the simplification of nc polys
"""
from . import polynomials as poly
from . import simplification_utils as su

def generate_measurements(label, io_config):
    """
    Returns a list of list of monomials that are the party's measurement ops
    label - name for the operators
    io_config - list of int indicating # of outputs for each input.
    """
    measurements = []
    for i, nout in enumerate(io_config):
        measurements.append(su.generate_operators(label + str(i), io_config[i] - 1, hermitian=True))

    return measurements

def projective_measurement_constraints(*parties):
    """
    Given a collection of parties measurement operators it returns the relevant
    constraints induced by projective measurements. I.e.
    - orthogonality
    - idempotency
    - commutativity
    """
    substitutions = {}

    # idempotency and orthogonality
    if isinstance(parties[0][0][0], list):
        parties = parties[0]
    for party in parties:
        for measurement in party:
            for projector1 in measurement:
                for projector2 in measurement:
                    if projector1 == projector2:
                        substitutions[projector1 * projector1] = projector1
                    else:
                        substitutions[projector1*projector2] = 0
                        substitutions[projector2*projector1] = 0
    # Projectors commute between parties in a partition
    for n1 in range(len(parties)):
        for n2 in range(n1+1, len(parties)):
            for measurement1 in parties[n1]:
                for measurement2 in parties[n2]:
                    for projector1 in measurement1:
                        for projector2 in measurement2:
                            substitutions[projector2*projector1] = \
                                projector1*projector2
    return substitutions
