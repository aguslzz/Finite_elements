#Code made for learning matrix management
import numpy as np
from matrix_functions import *
from math import cos, sin, pi, radians, degrees

nodes = 5
elements = 6                         # Link element type as shown in the problem graphs
# section_A = float(input("Input the area of the transversal section of the links (in^2)"))
section_A = 8
# load = float(input("Input the magnitude of the load(lb): "))
loads = np.zeros(nodes*2)
loads[7] = -500
loads[9] = -500
# E = float(input("Input the elasticity modulus(lb/in^2): "))
my_E = 1.90e6
# l1 = float(input("Input the lenght of links 1, 3, 4 and 6 (in): "))
l1 = 36 
# l1 = float(input("Input the lenght of links 2 and 5 (in): "))
l2 = 50.9
# Calculus of the equivalent k in each element k1 (for elements 1, 3, 4 and 6) and k2 (for elements 2 and 5)
k1 = float((((section_A)*(my_E))/l1))
k2 = float((((section_A)*(my_E))/l2))

# Global system variation among the local system
def K_e(theta):
    stheta= sin(theta)
    ctheta = cos(theta)
    K = np.array([[ctheta*ctheta, stheta*ctheta,(-1)*(ctheta*ctheta), (-1)*(stheta*ctheta)],
    [stheta*ctheta, stheta*stheta, (-1)*(stheta*ctheta), (-1)*(stheta*stheta)],
    [(-1)*(ctheta*ctheta), (-1)*(stheta*ctheta), ctheta*ctheta, stheta*ctheta],
    [(-1)*(stheta*ctheta), (-1)*(stheta*stheta), stheta*ctheta, stheta*stheta]])
    K[K==-0] = 0
    return K

# Stiffness matrixes 
big_k = np.zeros((nodes*2,nodes*2))          #stiffness global matrix
k_1 = k1*K_e(0)                              #stiffness matrix for elements 1, 3 and 6
k_4 = k1*K_e(radians(90))                    #stiffness matrix for element 4
k_4[np.abs(k_4) < 10**(-10)] = 0
k_2 = k2*K_e(radians(135))                   #stiffness matrix for element 2
k_5 = k2*K_e(radians(45))                    #stiffness matrix for element 5
k__5 = k_5[0][0]

for i in range(0, len(k_1)):                #filling element 1 stiffness
    big_k[0][i+0] += k_1[0][i]                
    big_k[1][i+0] += k_1[1][i]
    big_k[2][i+0] += k_1[2][i]
    big_k[3][i+0] += k_1[3][i]

for i in range(0, len(k_1)):                #filling element 3 stiffness
    big_k[4][i+4] += k_1[0][i]
    big_k[5][i+4] += k_1[1][i]
    big_k[6][i+4] += k_1[2][i]
    big_k[7][i+4] += k_1[3][i]

for i in range(0, len(k_1)):                #filling element 6 stiffness
    big_k[6][i+6] += k_1[0][i]
    big_k[7][i+6] += k_1[1][i]
    big_k[8][i+6] += k_1[2][i]
    big_k[9][i+6] += k_1[3][i]

big_k[3][3] += k1                            #filling element 4 stiffness
big_k[3][7] += (-1)*(k1)
big_k[7][3] += (-1)*(k1)
big_k[7][7] += k1

for i in range(0, len(k_1)):                #filling element 2 stiffness
    big_k[2][i+2] += k_2[0][i]
    big_k[3][i+2] += k_2[1][i]
    big_k[4][i+2] += k_2[2][i]
    big_k[5][i+2] += k_2[3][i]

big_k[2][2] += k__5                          #filling element 5 stiffness
big_k[2][3] += k__5
big_k[3][2] += k__5
big_k[3][3] += k__5
big_k[8][2] += (-1)*(k__5)                          
big_k[8][3] += (-1)*(k__5)  
big_k[9][2] += (-1)*(k__5)  
big_k[9][3] += (-1)*(k__5)  
big_k[2][8] += (-1)*(k__5)                          
big_k[2][9] += (-1)*(k__5)  
big_k[3][8] += (-1)*(k__5)  
big_k[3][9] += (-1)*(k__5)  
big_k[8][8] += k__5                          
big_k[8][9] += k__5
big_k[9][8] += k__5
big_k[9][9] += k__5
big_k[np.abs(big_k) < 10**(-10)] = 0

# Applying boundary conditions
big_k = np.delete(big_k, 0, 0)
big_k = np.delete(big_k, 0, 0)
big_k = np.delete(big_k, 2, 0)
big_k = np.delete(big_k, 2, 0)
big_k = np.delete(big_k, 0, 1)
big_k = np.delete(big_k, 0, 1)
big_k = np.delete(big_k, 2, 1)
big_k = np.delete(big_k, 2, 1)
loads = np.delete(loads, 0, 0)
loads = np.delete(loads, 0, 0)
loads = np.delete(loads, 2, 0)
loads = np.delete(loads, 2, 0)

# Solution
solution = (np.linalg.inv(big_k))@loads          #Global nodal displacement for elements 2,4 and 5. 1 and 3 dont have displacement because they are supported
print(solution)

# Postprocessing phase
element5 = np.zeros(4)                           #Element displacements
element5[0] = solution[0]
element5[1] = solution[1]
element5[2] = solution[4]
element5[3] = solution[5]

element2 = np.zeros(4)
element2[0] = solution[0]
element2[1] = solution[1]
element2[2] = 0
element2[3] = 0

element4 = np.zeros(4)
element4[0] = solution[0]
element4[1] = solution[1]
element4[2] = solution[2]
element4[3] = solution[3]

def K_einv(theta):
    stheta= sin(theta)
    ctheta = cos(theta)
    Kinv = np.array([[ctheta, stheta, 0, 0],
    [(-1)*(stheta), ctheta, 0, 0],
    [0, 0, ctheta, stheta],
    [0, 0, (-1)*(stheta), ctheta]])
    Kinv[Kinv==-0] = 0
    return Kinv

local5 = K_einv(radians(45))@element5                   #Local element and nodal displacement
local2 = K_einv(radians(135))@element2    
local4 = K_einv(radians(90))@element4    

stress5_x = (my_E)*((local5[0]-local5[2])/l2)           #Nodal stress in x and y direction
stress5_y = (my_E)*((local5[1]-local5[3])/l2)
stress2_x = (my_E)*((local2[0]-local2[2])/l2)
stress2_x = (my_E)*((local2[1]-local2[3])/l2)
stress4_x = (my_E)*((local4[0]-local4[2])/l2)
stress4_x = (my_E)*((local4[1]-local4[3])/l2)

print(local5[0])
print(local5[2])
print(stress5_x)
