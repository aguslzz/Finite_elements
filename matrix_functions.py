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

def gen_kmat(keq_arr, nodes):
    k_mat = np.zeros((nodes, nodes))
    for i in range(len(k_mat) - 1):
        k_mat[i, i] += keq_arr[i]
        k_mat[i, i + 1] -= keq_arr[i]
        k_mat[i + 1, i] -= keq_arr[i]
        k_mat[i + 1, i + 1] += keq_arr[i]
    return k_mat