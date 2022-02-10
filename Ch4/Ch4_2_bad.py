import math
import copy 

def euclid_dist(p, q=(0,0)):
    '''Computes Euclidean distance between points p and q.'''
    return math.sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)
        

def solution(dimensions, your_position, trainer_position, distance):
    # check direct shot
    if euclid_dist(your_position, trainer_position) > distance:
        return 0

    count = 1

    w = dimensions[0]
    h = dimensions[1]
    my_x = your_position[0]
    my_y = your_position[1]
    tr_x = trainer_position[0]
    tr_y = trainer_position[1]

    # Create queue for breath-first search of bounce tree
    queue = []

    walls = [0, 1, 2, 3] # L, B, R, T

    # Store distances from trainer to walls 
    corners = (0, 0, w, h)
    my_deltas = (my_x, my_y, w - my_x, h - my_y)
    tr_deltas = (tr_x, tr_y, w - tr_x, h - tr_y)

    # Add initial bounce to queue
    for wall_1 in walls:
        if my_deltas[(wall_1 + 1) % 4] != tr_deltas[(wall_1 + 1) % 4]:
            count += 1

        if wall_1 % 2 == 0:
            shot_path = [abs(my_deltas[wall_1] + tr_deltas[wall_1]), abs(my_y - tr_y)]
        else:
            shot_path = [abs(my_x - tr_x), abs(my_deltas[wall_1] + tr_deltas[wall_1])]
            
        if euclid_dist(shot_path) <= distance:
            for wall_2 in walls:
                if wall_1 != wall_2:
                    # Create ray path for checking hittability
                    ## we want to store the bounds of the wall_1 -> wall_2 light ray and their slope
                    source_1 = [0,0]
                    source_1[wall_2 % 2] += corners[wall_2]
                    source_1[(wall_2 + 1) % 2] += corners[(wall_2 + 1) % 4]
                    source_2 = [0,0]
                    source_2[wall_2 % 2] += corners[wall_2]
                    source_2[(wall_2 - 1) % 2] += corners[(wall_2 - 1) % 4]

                    # TODO: Split into zigzag and corner case for now but should be unifiable
                    if (wall_1 % 2) == (wall_2 % 2): # zigzag
                        ray_1 = [0, 0]
                        ray_1[(wall_2 + 1) % 2] += my_deltas[(wall_2 + 1) % 4]
                        ray_1[wall_2 % 2] += corners[2 + (wall_2 % 2)] + my_deltas[wall_2]

                        ray_2 = [0, 0]
                        ray_2[(wall_2 + 1) % 2] += my_deltas[(wall_2 - 1) % 4]
                        ray_2[wall_2 % 2] += corners[2 + (wall_2 % 2)] + my_deltas[wall_2]
                    else: # corner
                        # One bounce ray
                        ray_1 = [0, 0]
                        ray_1[(wall_2 + 1) % 2] += corners[(1 + wall_2) % 4] + my_deltas[(wall_2 + 1) % 2]
                        ray_1[wall_2 % 2] += my_deltas[wall_2 % 2]

                        # Two bounce ray
                        ray_2 = [0, 0]
                        ray_2[(wall_2 + 1) % 2] += corners[(3 + wall_2) % 4] + my_deltas[(wall_2 + 1) % 4]
                        ray_2[wall_2 % 2] += my_deltas[wall_2]
                    
                    queue.append((shot_path, [source_1, ray_1], [source_2, ray_2], [wall_1, wall_2]))            

    print(count)

    # Breath-first search of bounce tree
    while queue: 
        node = queue.pop()
        
        # check if distance if shorter than max distance
        wall_2 = node[-1][-1]
        shot_path = copy.copy(node[0])
        shot_path[wall_2 % 2] += tr_deltas[wall_2]
        dist = euclid_dist(shot_path)

        if dist <= distance:
            if verbose:
                print(node)
                print(dist)
            # add possible next bounces
            wall_1 = node[-1][-2]
            max_ray = node[1]

            # check if bounce can hit target
            if True:
                count += 1

            # Construct next possible bounces with correct rays

            #TODO: We assume that a zigzag can hit all three walls but this is not the case
            if (wall_1 % 2) == (wall_2 % 2): # Last bounce was zigzag
                for wall_new in walls:
                    if wall_new != wall_2:
                        max_ray_new = copy.copy(max_ray)
                        max_ray_new[wall_new % 2] += w if wall_new % 2 == 0 else h
                        queue.append((shot_path, max_ray_new, node[-1] + [wall_new]))
            else: # Last bounce was corner
                # check if we can hit both walls
                # L = 
                if True:
                    # add next corner bounce with direction of bounces
                    wall_new = (wall_2 + wall_2 - wall_1) % 4
                    max_ray_new = max_ray 
                    max_ray_new[wall_new % 2] += w if wall_new % 2 == 0 else h
                    queue.append((shot_path, max_ray, node[-1] + [wall_new]))
                    # add zigzag bounce
                    wall_new = (wall_2 + 2) % 4
                    max_ray_new = max_ray 
                    max_ray_new[wall_new % 2] += w if wall_new % 2 == 0 else h
                    queue.append((shot_path, max_ray, node[-1] + [wall_new]))
                else:
                    pass

    return count


if __name__ == '__main__':
    verbose = True

    test_cases = [
            # ([[3,2], [1,1], [2,1], 4], 7),
            ([[300,275], [150,150], [185,100], 500], 9)
        ]

    for test, result in test_cases:
        print(test)
        print(result)
        print(solution(*test))
        print('')