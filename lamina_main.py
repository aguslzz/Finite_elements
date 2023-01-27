#We will solve a mechanical resistance problem analyzing an iron laminate as a compound of springs, using hooke´s law and basic linear algebra

import numpy as np
from matrix_functions import *

# elements = float(input("Input the number of elements that you want to divide your laminate in (discretization): "))
elements = 4
# width = float(input("Input the top width of your laminate(in): "))
topwidth = 2
# width = float(input("Input the bottom width of your laminate(in): "))
bottomwidth = 1
# height = float(input("Input the height of your laminate(in): "))
height = 10
# thickness = float(input("Input the thickness of your laminate(in): "))
thickness = 0.125
# load = float(input("Input the magnitude of the load(lb): "))
load = 1000
# E = float(input("Input the elasticity modulus(lb/in^2): "))
my_E = 10.4e6

# Functions
def matrix_creator(n, value):                                   #creates a square matrix and fills with the required value
    my_matrix = np.zeros((n, n))
    for i in range(0, len(my_matrix)):
        my_matrix[i] = value
    return(my_matrix)
stiffness = matrix_creator(2,1)
print(stiffness)

def gen_kmat(keq_arr, nodes):                                   #generates k matrix as required in the exercise, Milo´s creation
    k_mat = np.zeros((nodes, nodes))
    for i in range(len(k_mat) - 1):
        k_mat[i, i] += keq_arr[i]
        k_mat[i, i + 1] -= keq_arr[i]
        k_mat[i + 1, i] -= keq_arr[i]
        k_mat[i + 1, i + 1] += keq_arr[i]
    return k_mat

# Equivalent area for each element
division = height/elements    #lenght of each element 
element_areas = np.zeros(elements+1)
for i in range(0, len(element_areas)):
    element_areas[i] = ((topwidth + (((bottomwidth-topwidth)/height)*(i*division)))*thickness)

# Stiffness constant for each element
elasticitys = np.zeros(elements)
for i in range(0, (len(element_areas)-1)):
    elasticitys[i] = (((element_areas[i]+element_areas[i+1])*my_E)/(2*division))

# Stiffness matrix for each element
stiff1 = matrix_creator(5,0)
for i in range(0,2):
    stiff1[0][i] = elasticitys[0]
    stiff1[i][0] = elasticitys[0]
    stiff1[i][i] = elasticitys[0]
stiff1[0][1] = elasticitys[0]*(-1)
stiff1[1][0] = elasticitys[0]*(-1)

stiff2 = matrix_creator(5,0)
for i in range(1,3):
    stiff2[0][i] = elasticitys[1]
    stiff2[i][0] = elasticitys[1]
    stiff2[i][i] = elasticitys[1]
stiff2[0,:] = 0
stiff2[:,0] = 0
stiff2[1][2] = elasticitys[1]*(-1)
stiff2[2][1] = elasticitys[1]*(-1)

stiff3 = matrix_creator(5,0)
for i in range(2,4):
    stiff3[0][i] = elasticitys[2]
    stiff3[i][0] = elasticitys[2]
    stiff3[i][i] = elasticitys[2]
stiff3[0,:] = 0
stiff3[:,0] = 0
stiff3[1,:] = 0
stiff3[:,1] = 0
stiff3[2][3] = elasticitys[2]*(-1)
stiff3[3][2] = elasticitys[2]*(-1)

stiff4 = matrix_creator(5,0)
for i in range(3,5):
    stiff4[4][i] = elasticitys[3]
    stiff4[i][4] = elasticitys[3]
    stiff4[i][i] = elasticitys[3]
stiff4[3][4] = elasticitys[3]*(-1)
stiff4[4][3] = elasticitys[3]*(-1)

# Object stiffness matrix
stiffness = matrix_creator(5,0)
stiffness = stiff1 + stiffness
stiffness = stiff2 + stiffness
stiffness = stiff3 + stiffness
stiffness = stiff4 + stiffness
# Another much simplier way jasjasj
stif = gen_kmat(elasticitys, elements+1)

# Creation of the missing matrixes 
loads = np.zeros(elements+1)      #Loads in each node
loads[-1] = load                  #Load is applied in the last node

# Proving the diference between solutions without contour conditions being applied
notsolution = (np.linalg.inv(stif))@loads   #Nodal displacements

# Solving and applying contour conditions
stif[0,:] = 0
stif[0,0] = 1
solution = (np.linalg.inv(stif))@loads   #Nodal displacements
solution[np.abs(solution) < 10**(-10)] = 0
print(solution)

# Postprocessing phase, stress in each element
stress = np.zeros(elements)
for i in range(len(stress)):
    stress[i] = ((my_E)*((solution[i+1]-solution[i])/division))
print(stress)
# Another way of calculating the stress in each element
prom_areas = np.zeros(elements)
for i in range(len(prom_areas)):
    prom_areas[i] = ((element_areas[i]+element_areas[i+1])/2)
stress2 = np.zeros(elements)
for i in range(len(stress2)):
    stress2[i] = loads[-1]/prom_areas[i]
print(stress2)

# Postprocessing phase, reactions in the object
reactions = stiffness@solution-loads
reactions[np.abs(reactions) < 10**(-10)] = 0


