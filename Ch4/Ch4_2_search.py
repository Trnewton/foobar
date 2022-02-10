'''
Passes:
1
2
4
5
6/7
'''


import math
import fractions

def euclid_dist(p, q=(0,0)):
    '''Computes Euclidean distance between points p and q.'''
    return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

def mirror_xy(deltas, n, m):
    # Compute x
    x = deltas[0] + 2 * (abs((n+1)//2) * deltas[2])
    x = x + 2 * (n//2) * deltas[0] if n >= 0 else -x - 2 * (abs(n+1)//2 * deltas[0])

    # Compute  y
    y = deltas[1] + 2 * (abs((m+1)//2) * deltas[3])
    y = y + 2 * ((m//2) * deltas[1]) if m >= 0 else -y - 2 * (abs(m+1)//2 * deltas[1])

    return (x, y)

def solution(dimensions, your_position, trainer_position, distance):
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
                            if (m_tr_dx * dy) == (m_tr_dy * dx) and (m_tr_dx) <= (dx):
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
                            if (m_my_dx * dy) == (m_my_dy * dx) and (m_my_dx) <= (dx):
                                good_flag = False
                                break

                if good_flag:
                    count += 1

    return count

if __name__ == '__main__':
    verbose = True

    test_cases = [
            ([[3,2], [1,1], [2,1], 4], 7),
            # ([[300,275], [150,150], [185,100], 500], 9)
        ]

    for test, result in test_cases:
        print(test)
        print(result)
        print(solution(*test))
        print('')