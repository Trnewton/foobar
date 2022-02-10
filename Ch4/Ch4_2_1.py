'''
Passes:
1
2
4
5
6/7
'''


import math


def euclid_dist(p, q=(0,0)):
    '''Computes Euclidean distance between points p and q.'''
    return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

def mirror_xy(deltas, n, m):
    # Compute x
    x = deltas[0] + 2 * (abs((n+1)//2) * deltas[2] + abs(abs(n+1)-0.5)//2 * deltas[0])
    x = x if n >= 0 else -x

    # Compute  y
    y = deltas[1] + 2 * (abs((m+1)//2) * deltas[3] + abs(abs(m+1)-0.5)//2 * deltas[1])
    y = y if m >= 0 else -y

    return (x, y)

def solution(dimensions, your_position, trainer_position, distance):
    # check direct shot
    if euclid_dist(your_position, trainer_position) > distance:
        return 0


    w = dimensions[0]
    h = dimensions[1]
    my_x = your_position[0]
    my_y = your_position[1]
    tr_x = trainer_position[0]
    tr_y = trainer_position[1]

    # Store distances from trainer to walls 
    my_deltas = (my_x, my_y, w - my_x, h - my_y)
    tr_deltas = (tr_x, tr_y, w - tr_x, h - tr_y)

    # [left, below, right, above]
    hyper_square = [int(math.ceil(float(distance - delta)/dimensions[n % 2]))
                    for n, delta in enumerate(my_deltas)]

    count = 0
    for n in range(-hyper_square[0], hyper_square[2]+1):
        for m in range(-hyper_square[1], hyper_square[3]+1):
            # Compute target x and y
            target_x, target_y = mirror_xy(tr_deltas, n, m)

            if euclid_dist(your_position, (target_x, target_y)) <= distance:
                dx = target_x - your_position[0]
                dy = target_y - your_position[1]

                good_flag = True
                for i in range(0, n + 1 - 2 * (n < 0), 1 - 2 * (n < 0)):
                    if not good_flag:
                        break
                    for j in range(0, m + 1 - 2 * (m < 0), 1 - 2 * (m < 0)):
                        # if n == 2 and m == -1:
                        #     print(n, m, i, j)
                        # Compute mirror me x and y
                        m_my_x, m_my_y = mirror_xy(my_deltas, i, j)
                        m_my_dx = m_my_x - your_position[0]
                        m_my_dy = m_my_y - your_position[1]

                        if m_my_dx == 0:
                            if dx == 0:
                                good_flag = False
                                break   
                        elif dy/dx == (m_my_dy/m_my_dx):
                            # bad
                            good_flag = False
                            break

                        if i == n and j == m:
                            break
                        
                        # Compute mirror trainer x and y
                        m_tr_x, m_tr_y = mirror_xy(tr_deltas, i, j)
                        m_tr_dx = m_tr_x - your_position[0]
                        m_tr_dy = m_tr_y - your_position[1]   
                        
                        if m_tr_dx == 0:
                            if dx == 0:
                                good_flag = False
                                break    
                        elif dy/dx == (m_tr_dy/m_tr_dx):
                            good_flag = False
                            break

                        if m_my_dx == 0:
                            if dx == 0:
                                good_flag = False
                                break   
                        elif dy/dx == (m_my_dy/m_my_dx):
                            # bad
                            good_flag = False
                            break
                    
                if good_flag:
                    # print(n, m)
                    count += 1

    return count

if __name__ == '__main__':
    verbose = True

    test_cases = [
            ([[3,2], [1,1], [2,1], 100], 7),
            # ([[300,275], [150,150], [185,100], 500], 9)
        ]

    for test, result in test_cases:
        print(test)
        print(result)
        print(solution(*test))
        print('')