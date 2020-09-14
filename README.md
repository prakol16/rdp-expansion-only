# Expansion only RDP
An implementation of the Douglas-Peucker algorithm, modified so that the simplified polygon contains the original.

Usage:

    rdp_closed(np.array([[0, 0], [-0.1, 1], [0, 2], [4, -1], [0, 0]]), 0.2)
    --> [[-0.1  ,  0.025],
         [-0.1  ,  2.075],
         [ 4.   , -1.   ],
         [-0.1  ,  0.025]]
         
