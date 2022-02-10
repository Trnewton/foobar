import math
import fractions
import time 
from collections import defaultdict

eucd_cache = {}
def euclidean_algo(a, b):
    # if (a, b) in eucd_cache:
    #     return eucd_cache[(a,b)]

    old_r,r = a,b
    old_s,s = 1,0
    old_t,t = 0,1

    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t 

    eucd_cache[(a,b)] = (old_r, old_s, old_t)

    return old_r, old_s, old_t

def euclid_dist(p, q=(0,0)):
    '''Computes Euclidean distance between points p and q.'''
    return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

def mirror_xy(deltas, n, m):
    # Compute x
    x = deltas[0] + 2 * (abs((n+1)//2) * deltas[2])
    x = x + 2 * ((n//2) * deltas[0]) if n >= 0 else -x - 2 * (abs(n+1)//2 * deltas[0])

    # Compute  y
    y = deltas[1] + 2 * (abs((m+1)//2) * deltas[3])
    y = y + 2 * ((m//2) * deltas[1]) if m >= 0 else -y - 2 * (abs(m+1)//2 * deltas[1])

    return (x, y)

def solution_slow_2(dimensions, your_position, trainer_position, distance):
    # check direct shot
    if euclid_dist(your_position, trainer_position) > distance:
        return 0

    w = fractions.Fraction(dimensions[0])
    h = fractions.Fraction(dimensions[1])
    my_x = fractions.Fraction(your_position[0])
    my_y = fractions.Fraction(your_position[1])
    tr_x = fractions.Fraction(trainer_position[0])
    tr_y = fractions.Fraction(trainer_position[1])

    # Store distances from trainer to walls 
    my_deltas = (my_x, my_y, w - my_x, h - my_y)
    tr_deltas = (tr_x, tr_y, w - tr_x, h - tr_y)

    # [left, below, right, above]
    hyper_square = [int(math.ceil(float(distance - delta)/dimensions[n % 2]))
                    for n, delta in enumerate(my_deltas)]

    hyper_width = hyper_square[0] + hyper_square[2] + 1
    hyper_height = hyper_square[1] + hyper_square[3] + 1
    hittable = [[False for m in range(hyper_height)] for n in range(hyper_width)]

    #### Iterate over possible shots to see if they are within range
    for n in range(-hyper_square[0], hyper_square[2]+1):
        for m in range(-hyper_square[1], hyper_square[3]+1):
            # Compute mirror targets x and y
            target = mirror_xy(tr_deltas, n, m)

            if euclid_dist(your_position, target) <= distance:
                hittable[n + hyper_square[0]][m + hyper_square[1]] = True
    
    
    # print 'Done distance check'

    #### Check the possible shots in range for self-hits and double hits on trainer
    count = 0
    for i, col in enumerate(hittable):
        n = i - hyper_square[0]
        for j, can_hit in enumerate(col):
            m = j - hyper_square[1]
            if can_hit:
                if n == 0 and m == 0:
                    count += 1
                    continue
                target = mirror_xy(tr_deltas, n, m)
                # To check to see if we are hitting ourselves or the trainer twice
                # we check to see either lay on the ray going to the mirror target
                dx = fractions.Fraction(target[0] - your_position[0])
                dy = fractions.Fraction(target[1] - your_position[1])

                if dx == 0 or dy == 0:
                    continue

                slope = fractions.Fraction(dy, dx)
                dx = slope.denominator if n > 0 else - slope.denominator
                dy = abs(slope.numerator) if m > 0 else - abs(slope.numerator)

                bad_flag = False
                # Check for collisions with self
                consts = [
                    0, # Even-Even self
                    2 * my_deltas[0] * dy, # odd-even self
                    (tr_deltas[1] - my_deltas[1]) * dx - (tr_deltas[0] - my_deltas[0]) * dy, # even-even trainer
                    (tr_deltas[1] - my_deltas[1]) * dx + (tr_deltas[0] + my_deltas[0]) * dy, # odd-even trainer
                    - 2 * my_deltas[1] * dx, # even-odd self
                    2 * my_deltas[0] * dy - 2 * my_deltas[1] * dx, # odd-odd self
                    -(tr_deltas[1] + my_deltas[1]) * dx - (tr_deltas[0] - my_deltas[0]) * dy, # even-odd trainer
                    -(tr_deltas[1] + my_deltas[1]) * dx + (tr_deltas[0] + my_deltas[0]) * dy, # odd-odd trainer
                    ]
                a = 2 * w * dy
                b = - 2 * h * dx
                d, x_star, y_star = euclidean_algo(a, b)                
                n_low, n_high = (0, n) if n > 0 else (n, 0)

                for i, c in enumerate(consts):
                    if c % d == 0:
                        p = c / d

                        x_0 = x_star * p
                        y_0 = y_star * p
                        n_j = 2 * x_0 - 1 * ((i % 2) == 1) # 1,3,5,7
                        m_j = 2 * y_0 - 1 * ((i // 4) == 1)

                        j = 0
                        while n_low < n_j:
                            j -= 1
                            n_j = 2 * (x_0 + (j if (b/d) > 0 else - j) * (b/d)) - 1 * ((i % 2) == 1)
                            m_j = 2 * (y_0 - (j if (b/d) > 0 else - j) * (a/d)) - 1 * ((i // 4) == 1)
                        
                        if ((i%4)//2 == 0): # self hit
                            while n_j <= n_high:
                                if 0 < abs(n_j) <= abs(n) and 0 < abs(m_j) <= abs(m):
                                    m_me_x, m_me_y  = mirror_xy(my_deltas, n_j, m_j)
                                    if your_position[0] < m_me_x < target[0] or  your_position[0] > m_me_x > target[0]:
                                        bad_flag = True
                                        break
                                j += 1
                                n_j = 2 * (x_0 + (j if (b/d) > 0 else - j) * (b/d)) - 1 * ((i % 2) == 1)
                                m_j = 2 * (y_0 - (j if (b/d) > 0 else - j) * (a/d)) - 1 * ((i // 4) == 1)
                        else: # double trainer
                            while n_j <= n_high:
                                if 0 <= abs(n_j) < abs(n) and 0 <= abs(m_j) < abs(m):
                                    m_tr_x, m_tr_y  = mirror_xy(tr_deltas, n_j, m_j)
                                    if your_position[0] < m_tr_x < target[0] or your_position[0] > m_tr_x > target[0]:
                                        bad_flag = True
                                        break
                                j += 1
                                n_j = 2 * (x_0 + (j if (b/d) > 0 else - j) * (b/d)) - 1 * ((i % 2) == 1)
                                m_j = 2 * (y_0 - (j if (b/d) > 0 else - j) * (a/d)) - 1 * ((i // 4) == 1)

                if not bad_flag:
                    count += 1

    return count

def solution(dimensions, your_position, trainer_position, distance):
    # check direct shot
    if euclid_dist(your_position, trainer_position) > distance:
        return 0

    t_0 = time.time()

    w = fractions.Fraction(dimensions[0])
    h = fractions.Fraction(dimensions[1])
    my_x = fractions.Fraction(your_position[0])
    my_y = fractions.Fraction(your_position[1])
    tr_x = fractions.Fraction(trainer_position[0])
    tr_y = fractions.Fraction(trainer_position[1])

    # Store distances from trainer to walls 
    my_deltas = (my_x, my_y, w - my_x, h - my_y)
    tr_deltas = (tr_x, tr_y, w - tr_x, h - tr_y)

    # [left, below, right, above]
    hyper_square = [int(math.ceil(float(distance - delta)/dimensions[n % 2]))
                    for n, delta in enumerate(my_deltas)]

    hyper_width = hyper_square[0] + hyper_square[2] + 1
    hyper_height = hyper_square[1] + hyper_square[3] + 1
    hittable = [[False for m in range(hyper_height)] for n in range(hyper_width)]

    #### Iterate over possible shots to see if they are within range
    for n in range(-hyper_square[0], hyper_square[2]+1):
        for m in range(-hyper_square[1], hyper_square[3]+1):
            # Compute mirror targets x and y
            target = mirror_xy(tr_deltas, n, m)

            if euclid_dist(your_position, target) <= distance:
                hittable[n + hyper_square[0]][m + hyper_square[1]] = True
    
    t_1 = time.time()
    print('sol2 t0:', t_1 - t_0)
    
    # print 'Done distance check'

    #### Check the possible shots in range for self-hits and double hits on trainer
    seen = {}
    count = 0
    for i, col in enumerate(hittable):
        n = i - hyper_square[0]
        for j, can_hit in enumerate(col):
            m = j - hyper_square[1]
            if can_hit:
                if n == 0 and m == 0:
                    count += 1
                    continue
                target = mirror_xy(tr_deltas, n, m)
                # To check to see if we are hitting ourselves or the trainer twice
                # we check to see either lay on the ray going to the mirror target
                dx = fractions.Fraction(target[0] - your_position[0])
                dy = fractions.Fraction(target[1] - your_position[1])

                if dx == 0 or dy == 0:
                    continue

                slope = fractions.Fraction(dy, dx)
                dx = slope.denominator if n > 0 else - slope.denominator
                dy = abs(slope.numerator) if m > 0 else - abs(slope.numerator)

                if (dx,dy) in seen:
                    continue
                seen[(dx,dy)] = True

                bad_flag = False
                # Check for collisions with self
                consts = [
                    fractions.Fraction(0), # Even-Even self
                    2 * my_deltas[0] * dy, # odd-even self
                    (tr_deltas[1] - my_deltas[1]) * dx - (tr_deltas[0] - my_deltas[0]) * dy, # even-even trainer
                    (tr_deltas[1] - my_deltas[1]) * dx + (tr_deltas[0] + my_deltas[0]) * dy, # odd-even trainer
                    - 2 * my_deltas[1] * dx, # even-odd self
                    2 * my_deltas[0] * dy - 2 * my_deltas[1] * dx, # odd-odd self
                    -(tr_deltas[1] + my_deltas[1]) * dx - (tr_deltas[0] - my_deltas[0]) * dy, # even-odd trainer
                    -(tr_deltas[1] + my_deltas[1]) * dx + (tr_deltas[0] + my_deltas[0]) * dy, # odd-odd trainer
                    ]
                a = 2 * w * dy
                b = - 2 * h * dx
                d, x_star, y_star = euclidean_algo(a, b)                
                n_low, n_high = (fractions.Fraction(0), fractions.Fraction(n)) if n > 0 else (fractions.Fraction(n), fractions.Fraction(0))
                m_low, m_high = (fractions.Fraction(0), fractions.Fraction(m)) if m > 0 else (fractions.Fraction(m), fractions.Fraction(0))

                for i, c in enumerate(consts):
                    if c % d == 0:
                        p = c / d

                        x_0 = x_star * p
                        y_0 = y_star * p
                        
                        t_n_low = (d/b) * (((n_low + (i % 2)) / 2) - x_0)
                        t_n_high = (d/b) * (((n_high + (i % 2)) / 2) - x_0)
    
                        t_m_low = - (d/a) * (((m_low + (i // 4)) / 2) - y_0)
                        t_m_high = - (d/a) * (((m_high + (i // 4)) / 2) - y_0)

                        int_n = [0,0]
                        if (d/b) > 0:
                            int_n[0] = t_n_low
                            int_n[1] = t_n_high
                        else:
                            int_n[0] = t_n_high
                            int_n[1] = t_n_low
                            
                        int_m = [0,0]
                        if (d/a) < 0:
                            int_m[0] = t_m_low
                            int_m[1] = t_m_high
                        else:
                            int_m[0] = t_m_high
                            int_m[1] = t_m_low

                        # Check endpoints 
                        if ((((d*b) > 0) != (n > 0)) != ((i%4)//2 == 0)):
                            # int_n[0] < t <= int_n[1]
                            b_1 = math.floor(int_n[1])
                            if (int_n[0] % 1) == 0:
                                a_1 = int_n[0] + 1
                            else:
                                a_1 = math.ceil(int_n[0])
                        else:
                            # int_n[0] <= t < int_n[1]
                            a_1 = math.ceil(int_n[0])
                            if (int_n[1] % 1) == 0:
                                b_1 = int_n[1] - 1
                            else:
                                b_1 = math.floor(int_n[1])
                        

                        if ((((d*a) > 0) != (m > 0)) != ((i%4)//2 == 0)):
                            # int_m[0] <= t < int_m[1]
                            a_2 = math.ceil(int_m[0])
                            if (int_m[1] % 1) == 0:
                                b_2 = int_m[1] - 1
                            else:
                                b_2 = math.floor(int_m[1])
                        else:
                            # int_m[0] < t <= int_m[1]
                            b_2 = math.floor(int_m[1])
                            if (int_m[0] % 1) == 0:
                                a_2 = int_m[0] + 1
                            else:
                                a_2 = math.ceil(int_m[0])


                        # Check for interval overlap
                        if (a_1 <= a_2 <= b_1 or a_1 <= b_2 <= b_1 or a_2 <= b_1 <= b_2) and (a_1 <= b_1) and (a_2 <= b_2):
                            bad_flag = True
                            break
                
                if not bad_flag:
                    count += 1
    t_2 = time.time()
    print('sol2 t1:', t_2 - t_1)

    return count

if __name__ == '__main__':
    import time 
    verbose = True
    test_2 = ([2, 100], [1, 99], [1, 1], 1000)

    start = time.time()
    sol_2 = solution(*test_2)
    end = time.time()
    print sol_2 
    print (end-start)

    # start = time.time()
    # sol_2 = solution_slow(*test_2)
    # end = time.time()
    # print sol_2 
    # print (end-start)

    quit()

    import time 
    import numpy as np
    np.random.seed(1)
    verbose = False
    dim_max = 1250
    dim_min = 10
    dist_max = 10000
    N = 100

    count = 0

    for n in range(N):
        if n % 10 == 0:
            print('n:', n)
            
        dim = [np.random.randint(dim_min, dim_max+1), np.random.randint(2, dim_max+1)]
        my = [np.random.randint(1, dim[0]), np.random.randint(1, dim[1])]
        tr = [np.random.randint(1, dim[0]), np.random.randint(1, dim[1])]
        while tr == my:
            tr = [np.random.randint(1, dim[0]), np.random.randint(1, dim[1])]

        dist = np.random.randint(1, min(max(dim)*20, dist_max))
        
        start_1 = time.time()
        sol_1 = solution(dim, my, tr, dist)
        end_1 = time.time()
        time_1 = end_1 - start_1

        start_2 = time.time()
        sol_2 = solution_slow(dim, my, tr, dist)
        end_2 = time.time()
        time_2 = end_2 - start_2

        if sol_1 != sol_2:
            print(dim, my, tr, dist)
            print(sol_1)
            print(sol_2)
            # break
        
        if time_1 > time_2:
            count += 1
            # print 'Slow'
            # print(dim, my, tr, dist)
        print(time_1)
        print(time_2)


    print count