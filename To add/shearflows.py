# -*- coding: utf-8 -*-
"""
Created on Sun Feb 24 11:02:15 2019

@author: Roy
"""

import math
import numpy as np
import random

#S_y/I_xzz
#I_xx=161250000 #mm^4
I_zz=1.565571e+09
I_yy=9.719095e+07
CGz=263.08492817627456
spar_height=205
large_pitch=81.33511390749766
small_pitch=8.512094985784858
stiffener_pitch=89.84720889328251
stiffener_area=40.56
C_a = 605 #mm
t_1=1.1
t_2=2.8

theta = 28
P = 97.4e+3
Ry = 100e+3
Rz = 100e+3
#Sum of forces in y and z
Sy = Ry-P*math.sin(math.radians(theta))
Sz = Rz+P*math.cos(math.radians(theta))

random.seed(9000)
boomarea=[]
for i in range (0,17):
    boomarea.append(random.uniform(100000,20000))
    

#print (boomarea)
angle=math.degrees(0.2012196329874848)
#print (angle)
z_pos= [560.2076385798177, 560.2076385798177, 470.6229157394533, 
         470.6229157394533, 381.0381928990888, 381.0381928990888, 
         291.45347005872435, 291.45347005872435, 201.86874721835989, 
         201.86874721835989, 112.28402437799541, 112.28402437799541, 
         100.0, 100.0, 37.51424233216027, 37.51424233216027, 0]
#print(z_pos)
y_pos= [8.869774538649953, -8.869774538649953, 26.609323615949847, 
        -26.609323615949847, 44.348872693249746, -44.348872693249746, 
        62.08842177054964, -62.08842177054964, 79.82797084784953, 
        -79.82797084784953, 97.56751992514943, -97.56751992514943, 
        100.0, -100.0, 78.07387583997614, -78.07387583997614, 0]

#print(y_pos)

boomarea2 = []
b1 = (t_1*stiffener_pitch)/6.0*(2.0+y_pos[1]/y_pos[0]) + stiffener_area + (t_1*stiffener_pitch)/6.0*(2.0+y_pos[2]/y_pos[0])
boomarea2.append(b1)
boomarea2.append(b1)
for i in range(0,4):
    j = 2+2*i
    b2 = (t_1*stiffener_pitch)/6.0*(2.0+y_pos[j-2]/y_pos[j]) + stiffener_area + (t_1*stiffener_pitch)/6.0*(2.0+y_pos[j+2]/y_pos[j])   
    boomarea2.append(b2)
    boomarea2.append(b2)    


q_bottom=[0]
q_prev = float(0)
for i in range(0,6):
    j = 11-(2*i)
    q=(-boomarea[j]*Sy/I_zz)*y_pos[j] - (boomarea[j]*Sz/I_yy)*(z_pos[j]-CGz) + q_prev
    q_prev = q
    q_bottom.append(q)
print(q_bottom)
print(q_bottom[-1])


q_top=[]
q_prev = q_bottom[-1]
for i in range(0,6):
    j = (2*i)
    q=(-boomarea[j]*Sy/I_zz)*y_pos[j] - (boomarea[j]*Sz/I_yy)*(z_pos[j]-CGz) + q_prev
    q_prev = q
    q_top.append(q)
print(q_top)


q_spar=[]
q=(-boomarea[13]*Sy/I_zz)*y_pos[13] - (boomarea[13]*Sz/I_yy)*(z_pos[13]-CGz)
q_spar.append(q)
print(q_spar)


q_curve=[0]
q16=(-boomarea[15]*Sy/I_zz)*y_pos[15] - (boomarea[15]*Sz/I_yy)*(z_pos[15]-CGz)
q_curve.append(q16)
q17=(-boomarea[16]*Sy/I_zz)*y_pos[16] - (boomarea[16]*Sz/I_yy)*(z_pos[16]-CGz)+q16
q_curve.append(q17)
q15=(-boomarea[14]*Sy/I_zz)*y_pos[14] - (boomarea[14]*Sz/I_yy)*(z_pos[14]-CGz)+q17
q_curve.append(q15)
print(q_curve)


d=spar_height*math.sin(math.radians(90-angle))
A1 = 1/4.0*math.pi*spar_height*spar_height*0.5
A2 = 0.5*spar_height*(C_a - 0.5*spar_height)

#Moment around boom 14 by top CCW positive
M1 = sum(q_top[0:5])*d*stiffener_pitch + q_top[5]*small_pitch*d

#Moment by q12
M2 = q_bottom[-1]*(y_pos[0]-y_pos[1])*(z_pos[0]-z_pos[14])

#Moment by curved part
M3 = -q_curve[1]*(y_pos[16]-y_pos[15])*(z_pos[13]-z_pos[15]) + q_curve[1]*(z_pos[15]-z_pos[16])*(y_pos[15]-y_pos[13])
M4 = -q_curve[2]*(y_pos[14]-y_pos[16])*(z_pos[13]-z_pos[16]) - q_curve[2]*(z_pos[14]-z_pos[16])*(y_pos[16]-y_pos[13])
M5 = -q_curve[3]*(y_pos[12]-y_pos[14])*(z_pos[13]-z_pos[14]) - q_curve[3]*(z_pos[12]-z_pos[14])*(y_pos[14]-y_pos[13])

#Moments due to q_s0
M6 = -2.0*A1 #times q_s01 
M7 = -2.0*A2 #times q_s02

#Moments by external forces
M = -Rz*0.5*spar_height - P*math.cos(math.radians(theta))*spar_height + P*math.sin(math.radians(theta))*(z_pos[13]-z_pos[16])

#Compute the angle of twist for curved cell (clockwise positive)
k1 = sum(q_curve[1:3])*stiffener_pitch/t_1 + q_curve[3]*large_pitch/t_1 - q_spar[0]*spar_height/t_2
l1 = 2.0*(stiffener_pitch+large_pitch)/t_1 + spar_height/t_2 #times qs01
m1 = -spar_height/t_2 #times qs02

#Compute the angle of twist for triangular cell (clockwise positive)
k2 = -sum(q_bottom)*stiffener_pitch/t_1 + -sum(q_top[0:5])*stiffener_pitch/t_1 + -q_top[5]*small_pitch/t_1 + q_spar[0]*spar_height/t_2
l2 = -spar_height/t_2 #times qs01
m2 = (2.0*small_pitch + 11.0*stiffener_pitch)/t_1 + spar_height/t_2

A = np.array([[M6,M7],[(l1-l2),(m1-m2)]])
B = np.array([[(M-(M1+M2+M3+M4+M5))],[(k2-k1)]])
print (np.linalg.solve(A, B ))
q_s0=np.linalg.solve(A,B)
q_s01=q_s0[0][0]
q_s02=q_s0[1][0]
