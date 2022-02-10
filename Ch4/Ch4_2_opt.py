def opt_searc():

    ring_width = hyper_square[0] + hyper_square[2] + 1
    ring_height = hyper_square[1] + hyper_square[3] + 1

    ring_count = 0
    while True:
        if verbose:
            print(hyper_square)
        # Check upper and lower sides
        dy_low = my_deltas[1] + tr_deltas[1] + 2 * ((hyper_square[1]//2)*tr_deltas[3] + ((hyper_square[1]+1)//2)*tr_deltas[1])
        dy_high = my_deltas[3] + tr_deltas[1] + 2 * ((hyper_square[3]//2)*tr_deltas[1] + ((hyper_square[3]+1)//2)*tr_deltas[3])
        for n in range(-hyper_square[0], hyper_square[2]+1):
            # compute dx
            if n >= 0:
                dx = abs(my_x - tr_x) + 2*(((n+1)//2)*tr_deltas[2] + (n//2)*tr_deltas[0])
            else: # n < 0
                dx = my_deltas[0] + tr_deltas[0] + 2 * ((abs(n-1)//2)*tr_deltas[2] + abs(n//2)*tr_deltas[0])

            # lower side            
            if euclid_dist(your_position, (dx, dy_low)) <= distance:
                if verbose:
                    print('a', n, -hyper_square[1])
                    # print(dx, dy_low)
                ring_count += 1

            # upper side
            if euclid_dist(your_position, (dx, dy_high)) <= distance:
                if verbose:
                    print('b', n, hyper_square[3])
                    # print(dx, dy_high)
                ring_count += 1

        # Check left and right sides
        dx_left = my_deltas[0] + tr_deltas[0] + 2 * (((hyper_square[0]-1)//2)*tr_deltas[2] + (hyper_square[0]//2)*tr_deltas[0])
        dx_right = abs(my_x - tr_x) + 2*(((hyper_square[2]+1)//2)*tr_deltas[2] + (hyper_square[2]//2)*tr_deltas[0])
        for n in range(1-hyper_square[1], hyper_square[3]):
            # compute dy
            if n == 0:
                dy = abs(my_y - tr_y)
            elif n > 0:
                dy = my_deltas[3] + tr_deltas[1] + 2 * ((n//2)*tr_deltas[1] + ((n-1)//2)*tr_deltas[3])
            else: # n < 0
                dy = my_deltas[1] + tr_deltas[1] + 2 * (abs(n//2)*tr_deltas[3] + abs((n+1)//2)*tr_deltas[1])

            # left side            
            if euclid_dist(your_position, (dx_left, dy)) <= distance:
                if verbose:
                    print('c', -hyper_square[0], n)
                    # print(dx_left, dy)
                ring_count += 1

            # right side
            if euclid_dist(your_position, (dx_right, dy)) <= distance:
                if verbose:
                    print('d', hyper_square[2], n)
                    # print(dx_right, dy)
                ring_count += 1
            

        if ring_count == (2 * ring_height + 2 * (ring_width-1)):
            # We are done
            print('done')
            count += ring_width * ring_height
            ring_height = 0
            ring_width = 0
            break
        else:
            count += ring_count
            ring_count = 0
            ring_height -= 2
            ring_width -= 2
            hyper_square = [dim-1 for dim in hyper_square]