from turtle import color
import numpy as np
import math
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
import fractions

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

def solution_slow(dimensions, your_position, trainer_position, distance):
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

    count = 0
    bad_bois = []
    # We iterate over possible shots to see if they are within range
    for n in range(-hyper_square[0], hyper_square[2]+1):
        for m in range(-hyper_square[1], hyper_square[3]+1):
            # Compute mirror targets x and y
            target = mirror_xy(tr_deltas, n, m)

            if euclid_dist(your_position, target) <= distance:
                if n == 0 and m == 0:
                    count += 1
                    continue

                # print('m, n:', m, n)
                # To check to see if we are hitting ourselves or the trainer twice
                # we check to see either lay on the ray going to the mirror target
                dx = fractions.Fraction(target[0] - your_position[0])
                dy = fractions.Fraction(target[1] - your_position[1])

                dx_ij = fractions.Fraction(my_deltas[2*(dx>0)])
                dy_ij = fractions.Fraction(my_deltas[1+2*(dy>0)])
                i = 0 
                j = 0
                good_flag = True
                while good_flag:
                    if i == n and j == m:
                        break
                    
                    # Compute second mirror targets x and y
                    m_tr_x, m_tr_y = mirror_xy(tr_deltas, i, j)
                    m_tr_dx = m_tr_x - your_position[0]
                    m_tr_dy = m_tr_y - your_position[1]
                    # Check if we have a collision
                    if dx == 0:
                        if m_tr_dx == 0:
                            good_flag = False
                            break
                    elif dy == 0:
                        if m_tr_dy == 0:
                            good_flag = False
                            break
                    else:
                        if euclid_dist(your_position, (m_tr_x, m_tr_y)) <= distance:
                            if (m_tr_dx * dy) == (m_tr_dy * dx) and abs(m_tr_dx) <= abs(dx):
                                good_flag = False
                                break
                    
                    # Compute next mirror square to check
                    if abs(dy * dx_ij) > abs(dy_ij * dx): # Beam hits top/bottom 
                        j += 1 if dy > 0 else -1
                        dx_ij -= (dy_ij * abs(dx/dy))
                        dy_ij = h
                    elif abs(dy * dx_ij) < abs(dy_ij * dx): # Beam hits side
                        i += 1 if dx > 0 else -1
                        dy_ij -= (dx_ij * abs(dy/dx))
                        dx_ij = w
                    else: # Beam hit corner
                        i += 1 if dx > 0 else -1
                        j += 1 if dy > 0 else -1
                        dx_ij = w
                        dy_ij = h
                    
                    # Compute mirror me x and y
                    m_my_x, m_my_y = mirror_xy(my_deltas, i, j)
                    m_my_dx = m_my_x - your_position[0]
                    m_my_dy = m_my_y - your_position[1]
                    
                    # Check if we have a collision
                    if dx == 0:
                        if m_my_dx == 0:
                            good_flag = False
                            break
                    elif dy == 0:
                        if m_my_dy == 0:
                            good_flag = False
                            break
                    else:
                        if euclid_dist(your_position, (m_my_x, m_my_y)) <= distance:
                            if (m_my_dx * dy) == (m_my_dy * dx) and abs(m_my_dx) <= abs(dx):
                                good_flag = False
                                break

                if good_flag:
                    count += 1
                else:
                    bad_bois.append((n,m))

    return count, bad_bois, max(hyper_square[0], hyper_square[2]), max(hyper_square[1], hyper_square[3])

def find_bad(N, M, dimensions, your_position, trainer_position, distance):
    w = fractions.Fraction(dimensions[0])
    h = fractions.Fraction(dimensions[1])
    my_x = fractions.Fraction(your_position[0])
    my_y = fractions.Fraction(your_position[1])
    tr_x = fractions.Fraction(trainer_position[0])
    tr_y = fractions.Fraction(trainer_position[1])

    # Store distances from trainer to walls 
    my_deltas = (my_x, my_y, w - my_x, h - my_y)
    tr_deltas = (tr_x, tr_y, w - tr_x, h - tr_y)

    count = 0
    bad_bois = []
    # We iterate over possible shots to see if they are within range
    for n in range(-N, N+1):
        for m in range(-M, M+1):
            # Compute mirror targets x and y
            target = mirror_xy(tr_deltas, n, m)
            if n == 0 and m == 0:
                count += 1
                continue

            # print('m, n:', m, n)
            # To check to see if we are hitting ourselves or the trainer twice
            # we check to see either lay on the ray going to the mirror target
            dx = fractions.Fraction(target[0] - my_x)
            dy = fractions.Fraction(target[1] - my_y)

            dx_ij = fractions.Fraction(my_deltas[2*(dx>0)])
            dy_ij = fractions.Fraction(my_deltas[1+2*(dy>0)])
            i = 0 
            j = 0
            good_flag = True
            while good_flag:
                if i == n and j == m:
                    break
                
                # Compute second mirror targets x and y
                m_tr_x, m_tr_y = mirror_xy(tr_deltas, i, j)
                m_tr_dx = m_tr_x - my_x
                m_tr_dy = m_tr_y - my_y
                # Check if we have a collision
                if dx == 0:
                    if m_tr_dx == 0:
                        good_flag = False
                        bad_bois.append((n,m,0))
                        break
                elif dy == 0:
                    if m_tr_dy == 0:
                        good_flag = False
                        bad_bois.append((n,m,0))
                        break
                else:
                    if (m_tr_dx * dy) == (m_tr_dy * dx) and 0 <= m_tr_dx / dx <= 1:
                        good_flag = False
                        bad_bois.append((n,m,1))
                        break
                
                # Compute next mirror square to check
                if abs(dy * dx_ij) > abs(dy_ij * dx): # Beam hits top/bottom 
                    j += 1 if dy > 0 else -1
                    dx_ij -= (dy_ij * abs(dx/dy))
                    dy_ij = h
                elif abs(dy * dx_ij) < abs(dy_ij * dx): # Beam hits side
                    i += 1 if dx > 0 else -1
                    dy_ij -= (dx_ij * abs(dy/dx))
                    dx_ij = w
                else: # Beam hit corner
                    i += 1 if dx > 0 else -1
                    j += 1 if dy > 0 else -1
                    dx_ij = w
                    dy_ij = h
                
                # Compute mirror me x and y
                m_my_x, m_my_y = mirror_xy(my_deltas, i, j)
                m_my_dx = m_my_x - my_x
                m_my_dy = m_my_y - my_y
                
                # Check if we have a collision
                if dx == 0:
                    if m_my_dx == 0:
                        good_flag = False
                        bad_bois.append((n,m,0))
                        break
                elif dy == 0:
                    if m_my_dy == 0:
                        good_flag = False
                        bad_bois.append((n,m,0))
                        break
                else:
                    if (m_my_dx * dy) == (m_my_dy * dx) and 0 <= m_my_dx / dx <= 1:
                        good_flag = False
                        bad_bois.append((n,m,-1))
                        break

            if good_flag:
                count += 1

    return count, bad_bois


def plot_mirror_lattice(problem, n, m):
    dims = problem[0]
    my_pos = problem[1]
    tr_pos = problem[2]
    ray_dist = problem[3]

    my_deltas = (my_pos[0], my_pos[1], dims[0] - my_pos[0], dims[1] - my_pos[1])
    tr_deltas = (tr_pos[0], tr_pos[1], dims[0] - tr_pos[0], dims[1] - tr_pos[1])

    fig = plt.figure(figsize=(16,16))

    for i in range(-n, n+2):
        plt.vlines(x=i*dims[0], ymin=-m*dims[1], ymax=(m+1)*dims[1])
        plt.annotate(f'{i-1}', ((i-0.5)*dims[0], -(m+0.25)*dims[1]), fontsize=20)
    for j in range(-m, m+2):
        plt.hlines(y=j*dims[1], xmin=-n*dims[0], xmax=(n+1)*dims[0])
        plt.annotate(f'{j-1}', (-(n+0.25)*dims[0], (j-0.5)*dims[1]), fontsize=20)

    for i in range(-n, n+1):
        for j in range(-m, m+1):
            m_my_pos = mirror_xy(my_deltas, i, j)
            m_tr_pos = mirror_xy(tr_deltas, i, j)

            plt.plot(m_my_pos[0], m_my_pos[1], 'g*')
            plt.plot(m_tr_pos[0], m_tr_pos[1], 'rs')


    plt.plot(0, 0, 'bo', markersize=20)

    plt.show()

def plot_mirror_lattice_slope(problem, n, m, slopes):

    colours = ['k', 'g', 'b']

    dims = problem[0]
    my_pos = problem[1]
    tr_pos = problem[2]
    ray_dist = problem[3]

    my_deltas = (my_pos[0], my_pos[1], dims[0] - my_pos[0], dims[1] - my_pos[1])
    tr_deltas = (tr_pos[0], tr_pos[1], dims[0] - tr_pos[0], dims[1] - tr_pos[1])

    fig = plt.figure(figsize=(16,16))

    for i in range(-n, n+2):
        plt.vlines(x=i*dims[0], ymin=-m*dims[1], ymax=(m+1)*dims[1])
        plt.annotate(f'{i-1}', ((i-0.5)*dims[0], -(m+0.25)*dims[1]), fontsize=10)
    for j in range(-m, m+2):
        plt.hlines(y=j*dims[1], xmin=-n*dims[0], xmax=(n+1)*dims[0])
        plt.annotate(f'{j-1}', (-(n+0.25)*dims[0], (j-0.5)*dims[1]), fontsize=10)

    for i in range(-n, n+1):
        for j in range(-m, m+1):
            m_my_pos = mirror_xy(my_deltas, i, j)
            m_tr_pos = mirror_xy(tr_deltas, i, j)

            plt.plot(m_my_pos[0], m_my_pos[1], 'g*')
            plt.plot(m_tr_pos[0], m_tr_pos[1], 'rs')


    plt.plot(0, 0, 'bo', markersize=10)
    ray_range = plt.Circle(my_pos, ray_dist, fill=False)
    plt.gca().add_patch(ray_range)
    plt.axis('off')

    for k, l, type in slopes:
        target = mirror_xy(tr_deltas, k, l)
        plt.plot((my_pos[0], target[0]), (my_pos[1], target[1]))

        rect = Rectangle((k*dims[0], l*dims[1]), dims[0], dims[1], alpha=0.2, color=colours[type])
        plt.gca().add_patch(rect)

    plt.show()

def get_rate(m, n, dimensions, your_position, trainer_position):

    w = fractions.Fraction(dimensions[0])
    h = fractions.Fraction(dimensions[1])
    my_x = fractions.Fraction(your_position[0])
    my_y = fractions.Fraction(your_position[1])
    tr_x = fractions.Fraction(trainer_position[0])
    tr_y = fractions.Fraction(trainer_position[1])

    # Store distances from trainer to walls 
    my_deltas = (my_x, my_y, w - my_x, h - my_y)
    tr_deltas = (tr_x, tr_y, w - tr_x, h - tr_y)

    target = mirror_xy(tr_deltas, n, m)
        
    dx = fractions.Fraction(target[0] - your_position[0])
    dy = fractions.Fraction(target[1] - your_position[1])

    slope = fractions.Fraction(dy, dx)
    dx = slope.denominator
    dy = slope.numerator

    dx_rates_my = [dx / (2*w), (dx + 2 * (1 - 2 * (n>0)) * my_deltas[2*(n<0)]) / (2*w)]
    dy_rates_my = [dy / (2*h), (dy + 2 * (1 - 2 * (m>0)) * my_deltas[1 + 2*(m<0)]) / (2*h)]

    # Check for double collisions with trainer
    dx_rates_tr = [
        (dx - tr_deltas[0] + my_deltas[0]) / (2*w),
        (dx - tr_deltas[0] + my_deltas[0] + 2 * (1 - 2 * (n>0)) * tr_deltas[2*(n<0)]) / (2*w),
        ]
    dy_rates_tr = [
        (dy - tr_deltas[1] + my_deltas[1]) / (2*h),
        (dy - tr_deltas[1] + my_deltas[1] + 2 * (1 - 2 * (m>0)) * tr_deltas[1 + 2*(m<0)]) / (2*h),
        ]

    return dx_rates_my, dy_rates_my, dx_rates_tr, dy_rates_tr

if __name__ == '__main__':
    test = ([4, 6], [3, 1], [3, 2], 3)
    n = 1
    m = 1

    sol, bad_bois = find_bad(n, m, *test)

    for boi in bad_bois:
        print(boi)
        rates = get_rate(boi[0], boi[1], *test[:-1])
        print('\t', rates[:2])
        print('\t', rates[2:])

    # plot_mirror_lattice(test, n, m)
    plot_mirror_lattice_slope(test, n, m, 
        bad_bois
        )