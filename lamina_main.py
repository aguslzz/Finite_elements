#Code made for learning matrix management
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


# Equivalent area for each element
division = height/elements    #lenght of each element 
element_areas = np.zeros(elements+1)
for i in range(0, len(element_areas)):
    element_areas[i] = ((topwidth + (((bottomwidth-topwidth)/height)*(i*division)))*thickness)
    i += 1
print(element_areas)

# Stiffness constant for each element
elasticitys = np.zeros(elements)
for i in range(0, (len(element_areas)-1)):
    elasticitys[i] = (((element_areas[i]+element_areas[i+1])*my_E)/(2*division))
print(elasticitys)

# Stiffness matrix for each element

    













