from fractions import Fraction
from collections import defaultdict
import numpy as np
from itertools import permutations


cycle_idx_sym_cache = {1:[((1,), 1)]}
def cycle_index_symmetric(n):
    if n in cycle_idx_sym_cache:
        return cycle_idx_sym_cache[n]

    cycle_index = defaultdict(Fraction)
    cycle_index[tuple(0 for n in range(int(n)-1)) + (1,)] = Fraction(1)
    for l in range(1, int(n)):
        cycle_index_nl = cycle_index_symmetric(n-l)
        for term, coef in cycle_index_nl:
            new_term = tuple(a_m if m+1 != l else a_m + 1 for m, a_m in enumerate(term)) \
            #    + tuple(0 if m+1 + len(term) != l  else 1 for m in range(len(term)-l))
            cycle_index[new_term] += coef

    cycle_index_list = cycle_index.items()
    cycle_index_list = [(idx, val/n) for idx, val in cycle_index_list]
    cycle_idx_sym_cache[n] = cycle_index_list

    return cycle_idx_sym_cache[n]

def cycle_index_prod_symetric(n, m):
    cycle_idx_n = cycle_index_symmetric(n)
    cycle_idx_m = cycle_index_symmetric(m)



def solution(w, h, s):
    w = Fraction(w)
    h = Fraction(h)
    s = Fraction(s)



if __name__ == '__main__':
    import time
    # cases = [[2,2,2], [2,3,4]]

    # for case in cases:
    #     print(solution(*case))

    n = Fraction(3)
    cycle_idx = cycle_index_symmetric(n)
    print([val for idx, val in cycle_idx])
    print([idx for idx, val in cycle_idx])
