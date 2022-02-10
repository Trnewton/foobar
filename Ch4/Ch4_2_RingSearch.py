from math import sqrt, ceil
from fractions import gcd


def euclid_dist(p, q=(0,0)):
    '''Computes Euclidean distance between points p and q.'''
    return sqrt((p[0]-q[0])**2 + (p[1]-q[1])**2)

def mirror_xy(deltas, n, m):
    '''Computes position of object with in mirror universe in the (n,m) mirror world.'''
    # Compute x
    x = deltas[0] + 2 * (abs((n+1)//2) * deltas[2])
    x = x + 2 * ((n//2) * deltas[0]) if n >= 0 else -x - 2 * (abs(n+1)//2 * deltas[0])

    # Compute  y
    y = deltas[1] + 2 * (abs((m+1)//2) * deltas[3])
    y = y + 2 * ((m//2) * deltas[1]) if m >= 0 else -y - 2 * (abs(m+1)//2 * deltas[1])

    return (x, y)

def solution(dimensions, your_position, trainer_position, distance):
    '''Computes number of way that we can hit a trainer.

    We use the observation that any ray reflection can be equivalently thought of
    as the ray traveling through the wall into a world that is mirrored along the
    wall from the current world. This means we can work in a lattice universe of
    mirrors of our original world. We  then first compute in which worlds our
    beam is long enough to hit the trainer. We then check if we have hit ourself
    before the trainer or 'hit' the trainer twice.

    To do the first part we can simply use the mirror lattice analogy to
    compute the distance the ray travels to hit the trainer and see if it is shorter
    than the allowed length. The easiest way to do this would to be to iterate
    over all possible worlds, but we can do better! By realizing that we only need
    to check the worlds at the end of the rays beam length we can iterate over
    the ellipsoid path of worlds that could possibly be out of range.

    The second part is done by noting that the path the beam leaving us travels
    is uniquely defined by the angle it leaves with. This means that instead of
    solving a set Diophantine equation to detect collisions we can simply check
    to see if two beams have the same angle. To make this efficient we can use
    the fact that each angle is defined by the x and y component of the right
    angle triangle it defines, as long as we also take into consideration the
    direction of the components (sign). Finally we must make sure we iterate
    over worlds in an order so that our we have already visited all worlds that
    our beam might stray though. We do this by following the quadrant oreder
    1, 2, 4, 3.
    '''
    # check direct shot
    if euclid_dist(your_position, trainer_position) > distance:
        return 0

    # Store distances to walls for the trainer and us
    my_deltas = (
        your_position[0], # Distance from me to left wall
        your_position[1], # Distance from me to bottom wall
        dimensions[0] - your_position[0], # Distance from me to right wall
        dimensions[1] - your_position[1] # Distance from me to top wall
        )
    tr_deltas = (
        trainer_position[0], # Distance from trainer to left wall
        trainer_position[1], # Distance from trainer to bottom wall
        dimensions[0] - trainer_position[0], # Distance from trainer to right wall
        dimensions[1] - trainer_position[1]# Distance from trainer to top wall
        )

    # Number of mirror worlds on each side of center in mirror universe
    # [left, below, right, above]
    hyper_square = [int(ceil(float(distance - delta)/dimensions[n % 2]))
                    for n, delta in enumerate(my_deltas)]


    #### Iterate over possible shots to see if they are within range
    # Initiate bounds of ellipse in each x slice of worlds to be empty
    ellipse_bounds = [[1,-1] for m in range(hyper_square[1] + hyper_square[3] + 1)]
    # start at bottom of ellipse and iterate clockwise until we get to the top
    n = 0
    m = -hyper_square[1]
    while n < 1 and m <= hyper_square[3]:
        mirror_tr = mirror_xy(tr_deltas, n, m)
        # If the current world has a hittable trainer we can move further out
        # from the center if not we move towards the center, depending on quadrant
        if euclid_dist(your_position, mirror_tr) <= distance:
            ellipse_bounds[m][0] = n
            n -= 1 if m <= 0 else 0
            m += 0 if m <= 0 else 1
        else:
            ellipse_bounds[m][0] = n + 1
            n += 0 if m <= 0 else 1
            m += 1 if m <= 0 else 0

    # start at top of ellipse and iterate clockwise until we get to the bottom
    n = 0
    m = hyper_square[3]
    while n > -1 and m >= -hyper_square[1]:
        mirror_tr = mirror_xy(tr_deltas, n, m)
        # If the current world has a hittable trainer we can move further out
        # from the center if not we move towards the center, depending on quadrant
        if euclid_dist(your_position, mirror_tr) <= distance:
            ellipse_bounds[m][1] = n
            n += 1 if m >= 0 else 0
            m -= 0 if m >= 0 else 1
        else:
            ellipse_bounds[m][1] = n - 1
            n -= 0 if m >= 0 else 1
            m -= 1 if m >= 0 else 0

    # Now to check hittability we just check if the world is inside the ellipse
    def hittable(n, m):
        '''Checks if a mirror world trainer is within range.'''
        return ellipse_bounds[m][0] <= n <= ellipse_bounds[m][1]


    #### Check the possible shots in range for self-hits and double hits on trainer
    def check_for_collisions(n, m, seen):
        '''Checks for any collisions in ray path from me to trainer in mirror world (n,m).'''
        collision = True
        if hittable(n, m):
            # Find position of trainer in mirror worlds
            mirror_tr = mirror_xy(tr_deltas, n, m)
            dx_tr = mirror_tr[0] - your_position[0]
            dy_tr = mirror_tr[1] - your_position[1]
            # Compute angle/slope of ray needed to hit trainer
            g_tr = abs(gcd(dx_tr, dy_tr)) # abs to preserve sign of dx and dy
            slope_tr = (dx_tr / g_tr, dy_tr / g_tr)

            # Find position of me in mirror worlds
            mirror_me = mirror_xy(my_deltas, n, m)
            dx_me = mirror_me[0] - your_position[0]
            dy_me = mirror_me[1] - your_position[1]
            # Compute angle/slope of ray needed to hit me
            g_me = max(abs(gcd(dx_me, dy_me)), 1) # abs to preserve sign of dx and dy
            slope_me = (dx_me / g_me, dy_me / g_me)

            # Add slopes to dict of seen slopes so we can check for future collisions
            # and check that we have not got a collision
            if abs(dx_me) <= abs(dx_tr): # We could hit ourself in this mirror world
                if slope_me not in seen:
                    seen[slope_me] = True

                if slope_tr not in seen:
                    seen[slope_tr] = True
                    collision = False
            else: # We could NOT hit ourself in this mirror world
                if slope_tr not in seen:
                    seen[slope_tr] = True
                    collision = False

                if slope_me not in seen:
                    seen[slope_me] = True

        return collision, seen

    seen = {}
    count = 0

    # We iterate over the quadrants of the mirror universe
    # 0 = top right, 1 = top left, 2 = bottom right, 3 = bottom left
    # we follow this order so that we check closest worlds first
    for Q in [0,1,2,3]:
        # Compute the direction of the diagonal traversals, start point,
        # and bounds of quadrant
        direction = [1 - 2 * (Q%2),  1 - 2 * (Q//2)]
        n = - (Q%2)
        m = - (Q//2)
        n_bounds = [n, (1 - 2 * (Q%2)) * hyper_square[2 * (1 - (Q%2))]]
        m_bounds = [m, (1 - 2 * (Q//2)) * hyper_square[1 + 2 * (1 - (Q//2))]]

        # If quadrant is empty continue
        if abs(n_bounds[0]) > abs(n_bounds[1]) or abs(m_bounds[0]) > abs(m_bounds[1]):
            continue

        # Check if starting world of quadrant has collision
        collision, seen = check_for_collisions(n, m, seen)
        count += 0 if collision else 1

        # While we are not at the farthest point of the current mirror universe quadrant
        while abs(n) < abs(n_bounds[1]) or abs(m) < abs(m_bounds[1]):
            # Compute the next direction to traverse in the mirror universe
            # which is either towards the x- or y- axis in anti-parallel
            # directions. Also compute which bounds to check for loop
            if n == n_bounds[1]:
                diag_dir = [-direction[0], direction[1]]
                mask = [False, True, True, False]
            elif m == m_bounds[1]:
                diag_dir = [direction[0], -direction[1]]
                mask = [True, False, False, True]
            elif n == n_bounds[0]:
                diag_dir = [direction[0], -direction[1]]
                mask = [True, False, False, False]
            elif m == m_bounds[0]:
                diag_dir = [-direction[0], direction[1]]
                mask = [False, False, True, False]

            # Iterate along diagonal component of mirror worlds
            while (n!=n_bounds[0] or mask[0]) and (n!=n_bounds[1] or mask[1]) \
                and (m!=m_bounds[0] or mask[2]) and (m!=m_bounds[1] or mask[3]):
                n += diag_dir[0]
                m += diag_dir[1]

                # Check if current world has a collision
                collision, seen = check_for_collisions(n, m, seen)
                count += 0 if collision else 1

            # We've hit a boundary of the mirror universe so we move along the boundry
            if (m == m_bounds[0] and abs(n) < abs(n_bounds[1])) or m == m_bounds[1]:
                n += direction[0]
            elif (n == n_bounds[0] and abs(m) < abs(m_bounds[1])) or n == n_bounds[1]:
                m += direction[1]

            # Check if current world has a collision
            collision, seen = check_for_collisions(n, m, seen)
            count += 0 if collision else 1

    return count

if __name__ == '__main__':
    import time
    verbose = True
    test_2 = ([4, 6], [3, 1], [3, 2], 3)

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