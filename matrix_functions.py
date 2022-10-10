#Code made for learning matrix management

import numpy as np

# Creates a square matrix and fills with the required value
def matrix_creator(n, value):
    my_matrix = np.zeros((n, n))
    for i in range(0, len(my_matrix)):
        my_matrix[i] = value
    return(my_matrix)
stiffness = matrix_creator(2,1)
print(stiffness)
