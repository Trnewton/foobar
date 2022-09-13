from fractions import Fraction, gcd
from collections import defaultdict

cycle_idx_sym_cache = {1:[((1,), 1)]}
def cycle_index_symmetric(n):
    '''Compute the cycle index polynomial for the symetric group of degree n.

    Uses the reccurrence relation to compute the cycle index of the degree n
    symetric group. The polynomial is represented by a list of of tuples where
    each tuple represents a term of the polynomial. The first element of each term
    is a tuple that represents the powers of the variables and the second the
    coefficient of the term. For example:
        ((0, 2, 1), 5) = 5 * a_2^2 * a_3
        ((1, 0, 0, 1), 1) = a_1 * a_4

    Note: This function caches results.

    Args:
        n (int) : Degree of symetric group to compute cycle index for

    Returns:
        cycle_index (list of tuples) : List of polynomial terms as described
    '''
    # Check for cached result
    if n in cycle_idx_sym_cache:
        return cycle_idx_sym_cache[n]

    # Use default dict to combinded terms with same variables
    cycle_index = defaultdict(Fraction)
    cycle_index[tuple(0 for n in range(n-1)) + (1,)] = Fraction(1)

    # Use recursive equation to add lower degree cycle indexes to get desired degree
    for l in range(1, n):
        cycle_index_nl = cycle_index_symmetric(n-l)

        # Multiply each term of lower degree cycle index by variable from recurrence
        # equation and add to cycle index we are computing
        for term, coef in cycle_index_nl:
            # Generate new term for cycle index
            new_term = []

            # Increase power of variable if it is not a new variable
            for m, a_m in enumerate(term, start=1):
                new_term += [a_m if m != l else a_m+1]

            # Add new variable to term if it is new
            for m in range(len(term)+1, n+1):
                new_term += [0 if m != l else 1]

            # Add term and coefficient to polynomial
            cycle_index[tuple(new_term)] += coef

    # Convert default dict to list of terms and cache it
    cycle_index_list = cycle_index.items()
    cycle_index_list = [(idx, val/n) for idx, val in cycle_index_list]
    cycle_idx_sym_cache[n] = cycle_index_list

    return cycle_idx_sym_cache[n]




def solution(w, h, s):
    '''
    We can recognize that the set of star configurations defines a set of matrices
    which we will denote by X. We can then further realize that the actions of
    switching rows and columns of a matrix form a group that acts on the set
    X, we call this group G. We finally note that the the equivalent configurations
    defined in the problem are exactly those matrices that are congruent to each
    other as defined by the action of G. This means that the number of unique
    matrices is equal to the number of unique orbits in X under G. We can thus apply
    the Cauchy-Frobenius Lemma to compute the solution to the problem.

    To apply the desired lemma we would sum the number of fixed x in X for each
    g in G but we can simplify this. We can note that G is exactly the Cartesian
    product of the two symmetric permutation groups of degree equal to w and h. To
    compute the number of unique orbits, the cycle index polynomial of a permutation
    group is used. The cycle index of a permutation group of degree n represents
    the types and number of permutations in a permutation group. It works on the
    principle that all permutations can be decomposed into cycles. It is possible
    to represent such a cycle decomposition as a multivariate monomial where each
    variable represents a cycle length and the exponents the number of cycles of
    that length. The index cycle polynomial is then the average of all such monomials.
    Evaluating this polynomial at an integer q then tells you the number of ways
    of ordering n things with q possible values such that if they permute they are
    considered equal.

    We can compute the cycle index polynomial of G by computing the cycle indexes
    of the symmetric groups that construct G. These cycle indexes are easy to compute
    since they have a recurrence equation. We can then compute the product of these
    two polynomials by looking at how two terms a_p^n and b_q^m would multiply.
    We realizing that the Cartesian pair of cycles has order equal to lcm(p,q)
    and so in the product there will be n*m * (p*q/lcm(p,q)) = n*m*gcd(p,q)
    such cycles in the decomposition. Thus we can compute the cycle index for G.
    '''

    # Compute cycle indexes for each symetric group of interest
    cycle_idx_w = cycle_index_symmetric(w)
    cycle_idx_h = cycle_index_symmetric(h)

    # Compute product of cycle indexes evaluated at s
    cycle_idx_at_s = 0
    for term_n, coeff_n in cycle_idx_w:
        for term_m, coeff_m in cycle_idx_h:
            # Get new coefficient of product term
            term_tot = coeff_n * coeff_m

            # Loop over each varible in the terms, skipping those that are zero
            for i, var_n in enumerate(term_n):
                if var_n != 0:
                    for j, var_m in enumerate(term_m):
                        if var_m != 0:
                            # Compute resulting product term using theory of cycles
                            # in product spaces as described above
                            term_tot *= (s) ** (var_n * var_m * gcd(i+1, j+1))

            cycle_idx_at_s += term_tot

    return str(cycle_idx_at_s)

if __name__ == '__main__':
    cases = [[2,2,2], [2,3,4]]

    for case in cases:
        print(solution(*case))

    # n = Fraction(3)
    # cycle_idx = cycle_index_symmetric(n)
    # print([val for idx, val in cycle_idx])
    # print([idx for idx, val in cycle_idx])
