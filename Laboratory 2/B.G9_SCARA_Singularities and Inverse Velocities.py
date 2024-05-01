B. Singularities and Inverse Velocities python file

import numpy as np
import sympy as sp
 
# link lengths in mm
a1 = float(input("a1 = "))
a2 = float(input("a2 = "))
a3 = float(input("a3 = "))
a4 = float(input("a4 = "))
a5 = float(input("a5 = "))

T1 = float(input("T1 = "))
T2 = float(input("T2 = "))
d3 = float(input("d3 = "))

# degree to radian
T1 = (T1/180.0)*np.pi
T2 = (T2/180.0)*np.pi

# Parametric Table (theta, alpha, r, d)
PT = [[T1,(0.0/180.0)*np.pi,a2,a1],
   [T2,(180.0/180.0)*np.pi,a4,a3],
          [(0.0/180.0)*np.pi,(0.0/180.0)*np.pi,0,a5+d3]] 

# HTM formulae
i = 0
H0_1 = [[np.cos(PT[i][0]),-np.sin(PT[i][0])*np.cos(PT[i][1]),np.sin(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.cos(PT[i][0])],
            [np.sin(PT[i][0]),np.cos(PT[i][0])*np.cos(PT[i][1]),-np.cos(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.sin(PT[i][0])],
            [0,np.sin(PT[i][1]),np.cos(PT[i][1]),PT[i][3]],
            [0,0,0,1]]

i = 1
H1_2 = [[np.cos(PT[i][0]),-np.sin(PT[i][0])*np.cos(PT[i][1]),np.sin(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.cos(PT[i][0])],
            [np.sin(PT[i][0]),np.cos(PT[i][0])*np.cos(PT[i][1]),-np.cos(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.sin(PT[i][0])],
            [0,np.sin(PT[i][1]),np.cos(PT[i][1]),PT[i][3]],
            [0,0,0,1]]

i = 2
H2_3 = [[np.cos(PT[i][0]),-np.sin(PT[i][0])*np.cos(PT[i][1]),np.sin(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.cos(PT[i][0])],
            [np.sin(PT[i][0]),np.cos(PT[i][0])*np.cos(PT[i][1]),-np.cos(PT[i][0])*np.sin(PT[i][1]), PT[i][2]*np.sin(PT[i][0])],
            [0,np.sin(PT[i][1]),np.cos(PT[i][1]),PT[i][3]],
            [0,0,0,1]]
    
H0_1 = np.matrix(H0_1)
H1_2 = np.matrix(H1_2)
H2_3 = np.matrix(H2_3)

H0_2 = np.dot(H0_1,H1_2)
H0_3 = np.dot(H0_2,H2_3)

# Jacobian Matrix Program

#1. Linear / Translational Vectors
Z_1 = [[0],[0],[1]] # The [0,0,1] vector

#Row 1 to 3, Column 1
J1a = [[1,0,0],
           [0,1,0],
           [0,0,1]] #R0_0
J1a = np.dot(J1a,Z_1)

J1b_1 =H0_3[0:3,3:] 
J1b_1 = np.matrix(J1b_1)

J1b_2 = [[0],[0],[0]] # position vector at the base
J1b_2 = np.matrix(J1b_2)

J1b = J1b_1 - J1b_2

J1 = [[(J1a[1,0]*J1b[2,0])-(J1a[2,0]*J1b[1,0])],
          [(J1a[2,0]*J1b[0,0])-(J1a[0,0]*J1b[2,0])],
          [(J1a[0,0]*J1b[1,0])-(J1a[1,0]*J1b[0,0])]]
J1 = np.matrix(J1)

#Row 1 to 3, Column 2
J2a = H0_1[0:3,0:3]
J2a = np.dot(J2a,Z_1)

J2b_1 = H0_3[0:3,3:]
J2b_1 = np.matrix(J2b_1)

J2b_2 = H0_1[0:3,3:]
J2b_2 = np.matrix(J2b_2)

J2b = J2b_1 - J2b_2

J2 = [[(J2a[1,0]*J2b[2,0])-(J2a[2,0]*J2b[1,0])],
          [(J2a[2,0]*J2b[0,0])-(J2a[0,0]*J2b[2,0])],
          [(J2a[0,0]*J2b[1,0])-(J2a[1,0]*J2b[0,0])]]
J2 = np.matrix(J2)

#Row 1 to 3, Column 3
J3 = H0_2[0:3,0:3]
J3 = np.dot(J3,Z_1)
J3 = np.matrix(J3)

#2. Rotation / Orientation Vectors
 
#Row 4 to 6, Column 1
J4 = J1a
J4 = np.matrix(J4)

#Row 4 to 6, Column 2
J5 = J2a
J5 = np.matrix(J5)

#Row 4 to 6, Column 3
J6 = [[0],[0],[0]]
J6 = np.matrix(J6)

#3. Concatenated Jacobian Matrix
JM1 = np.concatenate((J1,J2,J3),1)
JM2 = np.concatenate((J4,J5,J6),1)

J = np.concatenate((JM1,JM2),0)
J = np.matrix(J)

## Singularity
D_J = np.linalg.det(JM1)
print("D_J = ",D_J)

## Inverse Velocity
I_V = np.linalg.inv(JM1)
print("I_V = ",I_V)
