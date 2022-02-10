from fractions import Fraction
from fractions import Fraction
from collections import defaultdict

import numpy as np
import networkx as nx
from pyvis.network import Network as pvnet

PRIMES = [
            2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 
            67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 
            137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197,
            199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 
            277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 
            359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 
            439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 
            521, 523, 541
        ]
NUM_P = len(PRIMES)

def dyadic_orbit_check(a, b):
    '''Checks if two trainers will wrestle forever (true) or not (false).
    
    By associating the thumb war game with the Dyadic Transformation
    we can check for fixed points by seeing if the initial condition
    has a finite binary representation once normalized.
    
    Args:
        a (int): Number of bananas first trainer has
        b (int): Number of bananas second trainer has
    Returns
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
        adj_list (list): list of trainer pairs that will fight forever
    '''
    adj_dict = defaultdict(list)
    for n, a in enumerate(banana_list):
        for m, b in enumerate(banana_list[n:]):
            if dyadic_orbit_check(a, b):
                adj_dict[n].append(n+m)
                adj_dict[n+m].append(n)
    return adj_dict
    
def compute_adj_mat(banana_list):
    N = len(banana_list)
    adj_mat = np.zeros((N,N))
    for n, a in enumerate(banana_list):
        for m, b in enumerate(banana_list[n:]):
            if dyadic_orbit_check(a,b):
                adj_mat[n,n+m] = 1
                adj_mat[n+m,n] = 1
    return adj_mat
    
def matching_number(adj, N):
    '''Determines matching number of stable pairings graph. 
    
    If matrix is set with only 1s and -1s then we will most likely get 
    duplicate rows and hence zero eigenvalues, so as a "fix" we randomize
    the enteries to make it unlikely that rows are linearly dependant. 
    Should implement proper fix....
    '''
    # Convert adj list into a skew-symetric matrix that realizes the given graph
    adj_mat = np.zeros((N,N))
    count = 0
    for n, m in adj:
        rand_int = np.random.randint(1,10)
        adj_mat[n,m] = PRIMES[count]
        adj_mat[m,n] = -PRIMES[count]
        count += 1
    
    # Compute matching number of graph by finding number of non-zero eigenvalues
    matching_number = np.sum(~np.isclose(np.linalg.eigvals(adj_mat), 0)) // 2
    return matching_number

def maximal_matching(G, N):
    '''
    Args:
        G (dict of int: set): Ajacency map of G
        N (int): Number of nodes in G
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

                        for z, rho_z in enumerate(rho):
                            if rho_z in S_x_set or rho_z in S_y_set:
                                rho[z] = u

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
    if len(banana_list) < 2:
        return 1
    if len(banana_list) < 3:
        return 2*int(~dyadic_orbit_check(banana_list[0], banana_list[1]))

    adj_dict = compute_adj_dict(banana_list)

    return len(banana_list) - 2 * len(maximal_matching(adj_dict, len(banana_list)))

def test_solution(banana_list):
    result = solution(banana_list)

    G = nx.from_numpy_array(compute_adj_mat(banana_list))
    nx_result = len(banana_list) - 2 * len(nx.max_weight_matching(G))

    # print(result)
    # print(nx_result)

    return result == nx_result

if __name__ == "__main__":
    # banana_list = [7, 3, 5, 3, 8, 6, 1]

    # print(compute_adj_mat(banana_list))

    # print(test_solution(banana_list))

    # quit()

    N = 1000
    max_trainer = 100
    max_banana = 1073741823
    b_list = []
    for n in range(N):
        b_list.append(np.random.randint(1, max_banana, size=np.random.randint(1, max_trainer)))

    count = 0
    for b in b_list:
        if not test_solution(b):
            print('BAD')
            count += 1
    
    print(count)