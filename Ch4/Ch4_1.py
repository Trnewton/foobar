# -*- coding: utf-8 -*-
'''
This solution consists of two parts; first checking whether a pair of trainers will thumb
wrestle forever and second finding a maximal matching of trainers that will fight forever.

We can solve the first problem efficiently by realizing that the thumb war games state is 
completely determined by the number of bananas one of the trainers has and the total number
of bananas. This means we can analysis the system of difference  equations by looking at 
only one of the trainers banana stash. We can then perform a change of variables to 
normalize the resulting system and realize that we are dealing with the Dyadic 
Transformation.This means we can use the binary representation of the initial conditions
to determine if the system will orbit or not.

To solve the second problem we can realize that the possible pairings of trainer that 
loop forever can be thought of as edges in a graph. This then means the problem of 
distracting the greatest number of trainers is equivalent to finding the maximal cardinal
matching of the resulting graph. We use Edmonds Blossom algorithm to compute this.
'''

import numpy as np
from fractions import Fraction
from collections import defaultdict

def dyadic_orbit_check(a, b):
    '''Checks if two trainers will wrestle forever (true) or not (false).
    
    By associating the thumb war game with the Dyadic Transformation
    we can check for fixed points by seeing if the initial condition
    has a finite binary representation once normalized.
    
    Args:
        a (int): Number of bananas first trainer has
        b (int): Number of bananas second trainer has
    Returns:
        bool: True iff the pairing of thumb warriors dooesn't reach a fixed point    
    '''
    T = a + b
    f = Fraction(a, T)

    denom = f.denominator
    return (np.log2(denom)%1) != 0

def compute_adj_dict(banana_list):
    '''Computes the pairings of trainers that will fight forever.
    
    Args:
        banana_list (list): List of initial banana holdings by trainers
    Returns:
        adj_list (dict of int:lists): dict of list of trainer pairs that will fight forever
    '''
    adj_dict = defaultdict(list)
    for n, a in enumerate(banana_list):
        for m, b in enumerate(banana_list[n:]):
            if dyadic_orbit_check(a, b):
                adj_dict[n].append(n+m)
                adj_dict[n+m].append(n)
    return adj_dict
    
def maximal_matching(G, N):
    '''Uses Edmonds Blossom algorithm to compute maximum matching set of G.
    
    Reference: László Lovász and M. D. Plummer. Matching Theory, volume 121 of 
                Annals of Discrete Mathematics. North Holland, 1986.

    Args:
        G (dict of int: set): Ajacency map of G
        N (int): Number of nodes in G
    Returns:
        matching (list of tuples): List of the bidirectional matchings in the 
            maximum matching.
    '''

    unscanned = [n for n in range(N)]
    mu = np.arange(N)
    phi = np.arange(N)
    rho = np.arange(N)

    def is_matched(x):
        return mu[x] != x

    def is_blossom(x):
        return ~is_matched(x) or phi[mu[x]] != mu[x]

    def is_base(x):
        return is_blossom(x) and phi[x] == x
    
    def is_out_of_forest(x):
        return mu[x] != x and phi[mu[x]] == mu[x] and phi[x] == x
    
    def S(x):
        '''Creates M-alternating path for x.'''
        S_x = [x]
        while True:
            
            if mu[x] != x:
                x = mu[x] 
                S_x.append(x)
            else:
                return S_x
            
            if phi[x] != x:
                x = phi[x]
                S_x.append(x)
            else:
                return S_x

    has_blossom = True
    while has_blossom:
        # Check if blossom is unscanned
        has_blossom = False
        for x in unscanned:
            if is_blossom(x):
                has_blossom = True
                break

        if has_blossom:
            # Check if blossom has neighbour to work with
            break_outter = False
            for y in G[x]:
                if is_blossom(y) and rho[y] != rho[x]:
                    S_x = S(x)
                    S_y = S(y)

                    if not (set(S_x) & set(S_y)): # Found augmented path
                        # Augment matching
                        for i, x_2i in enumerate(S_x[2::2], start=1):
                            mu[S_x[2*i-1]] = x_2i
                            mu[x_2i] = S_x[2*i-1]

                        for i, y_2i in enumerate(S_y[2::2], start=1):
                            mu[S_y[2*i-1]] = y_2i
                            mu[y_2i] = S_y[2*i-1]

                        mu[x] = y
                        mu[y] = x

                        # Reset maps
                        phi = np.arange(N)
                        rho = np.arange(N)
                        unscanned = [n for n in range(N)]

                        # Break to outter most loop
                        break_outter = True
                        break
                    else: # Found blossom
                        # Find first base point shared by S(x) and S(y)
                        S_x_set = set(S_x)
                        S_y_set = set(S_y)
                        for u in S_x:
                            if u in S_y_set and is_base(u):
                                break
                        
                        # Update maps to shrink blossom
                        for i, x_2i in enumerate(S_x[2::2], start=1):
                            if x_2i in S_y_set:
                                break
                            phi[x_2i] = S_x[2*i-1]

                        for i, y_2i in enumerate(S_y[2::2], start=1):
                            if y_2i in S_x_set:
                                break
                            phi[y_2i] = S_y[2*i-1]
                        
                        phi[x] = y
                        phi[y] = x
                        
                        # Update base point
                        for z, rho_z in enumerate(rho):
                            if rho_z in S_x_set or rho_z in S_y_set:
                                rho[z] = u
                
                # Extend tree
                elif is_out_of_forest(y):
                    phi[y] = x

            if not break_outter: # Remove x from unscanned
                unscanned.remove(x)

    matching = []
    seen = defaultdict(list)
    for x, y in enumerate(mu):
        if x != y and not seen[x]:
            matching.append((x, y))
        seen[x] = True
        seen[y] = True

    return matching

def solution(banana_list):
    # Check edges
    if len(banana_list) < 2:
        return 1
    if len(banana_list) < 3:
        return 2*int(~dyadic_orbit_check(banana_list[0], banana_list[1]))
        
    adj_dict =  compute_adj_dict(banana_list)

    return len(banana_list) - 2 * len(maximal_matching(adj_dict, len(banana_list)))

if __name__ == '__main__':
    b_list = [[1,1],[1,7,3,21,13,19]]

    count = 0
    for b in b_list:
        print(solution(b))