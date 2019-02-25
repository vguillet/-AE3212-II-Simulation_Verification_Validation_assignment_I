import numpy as np
import matplotlib.pyplot as plt
from src.class_SLICE import Slice

from config import *
#get_ipython().run_line_magic('matplotlib','qt')

# CONSTANTS #
cross_section = Slice()

theta_max = theta*np.pi/180    # [rad]
P = p

I1_zz = cross_section.I_yy    # [mm**4]
I1_yy = cross_section.I_zz    # [mm**4]

sc = 106.12931274441615     # [mm] shear center distance from LE
hinge = ha/2
act = hinge*np.sqrt(2)

arm_a = act*np.sin(np.pi/4-theta_max)-(sc - hinge)*np.sin(theta_max) 
arm_q = (ca/4-sc)*np.cos(theta_max) 
arm_y = (sc-hinge)*np.cos(theta_max)
arm_z = (sc-hinge)*np.sin(theta_max)

h_arm_a = act*np.sin(np.pi/4-theta_max)
h_arm_q = (ca/4-hinge)*np.cos(theta_max)

def Y_force(n_steps=1000, plot=True, info=True, ret=False):

    # INITIATION #
    dx = la/(n_steps)
    x_pos = np.arange(0, la+1, dx)
    n_x1 = int(n_steps*x1/la)+1
    n_x2 = int(n_steps*x2/la)+1
    n_x3 = int(n_steps*x3/la)+1
    xa1 = x2 - xa/2
    xa2 = x2 + xa/2
       
    F1 = 25000 # [N] - INITIAL GUESS
    correction = 0
    error = 100
    itteration = 0
    itteration_factor = 0.03 # 0.03 works very well
    y0 = 0
    a0 = 0
    
    # ITTERATION IN Y DIRECTION #
    while abs(error) >= 10**-10: # [mm]
        itteration += 1
        
        shear = [0]*len(x_pos)
        moment = [0]*len(x_pos)
        y_disp = [0]*len(x_pos)
        angle = [0]*len(x_pos)
        
        y_disp[0] = y0
        angle[0] = a0
        
        F1 += correction
        F3 = (F1*(x2-x1) + q*la*(la/2-x2))/(x3-x2)
        F2 = F1 + F3 - q*la
        
        for i in range(1, n_steps+1):
                        
            # SHEAR INTEGRATION #        
            shear[i]=shear[i-1] + q*dx        
            if i == n_x1:
                shear[i] = shear[i] - F1
            if i == n_x2:
                shear[i] = shear[i] + F2
            if i == n_x3:
                shear[i] = shear[i] - F3
                
            # MOMENT INTEGRATION #    
            moment[i] = moment[i-1] - shear[i]*dx
            
            # ANGLE INTEGRATION #
            angle[i] = angle[i-1] + moment[i]*dx/(E*I1_zz)
            
            # DISPLACEMENT INTEGRATION #
            y_disp[i] = y_disp[i-1] + angle[i]*dx
        
        # DISPLACEMENT AND ANGLE ADJUSTMENT #
        x2_adjustment = - y_disp[n_x2]
        x3_adjustment = - y_disp[n_x3] + d3 +y_disp[n_x2]
        angle_adjustment = x3_adjustment/(x3-x2)
        for i in range(0, n_steps+1):
            y_disp[i] = y_disp[i] + x2_adjustment + x3_adjustment*(x_pos[i]-x2)/(x3-x2)
            angle[i] = angle[i] + angle_adjustment
            
        correction = F1*(d1 - y_disp[n_x1])*itteration_factor
        y0 = y_disp[0]
        a0 = angle[0]
        error = y_disp[n_x1] - d1
    
    if info:    
        print("Calculated forces in y direction in "+str(itteration)+" itterations, final error: "+str(error)+" mm")
        #print("F1y = "+str(F1)+" F2y = "+str(F2)+" F3y = "+str(F3))
    if plot:
        par = shear
        axis = [0]*len(x_pos)       
        h = np.arange(min(par),max(par),(max(par)-min(par))/10)
        h1x = [x1]*len(h)
        h2x = [x2]*len(h)
        h3x = [x3]*len(h)
        h1a = [xa1]*len(h)
        h2a = [xa2]*len(h)
        
        plt.figure()        
        plt.plot(x_pos, par)
        plt.plot(h1x, h, h2x, h, h3x, h, c="g")
        plt.plot(h1a, h, h2a, h, c="y")
        plt.plot(x_pos, axis, c="r")
        plt.savefig("fig",dpi=1000)
        plt.show()
    if ret:
        return([F1, 0, -F2, 0, F3])
        
    

def Z_force(FY, n_steps=1000, plot=True, info=True, ret=False):
    
    # INITIATION #
    dx = la/(n_steps)
    x_pos = np.arange(0, la+1, dx)
    n_x1 = int(n_steps*x1/la)+1
    n_x2 = int(n_steps*x2/la)+1
    n_x3 = int(n_steps*x3/la)+1
    xa1 = x2 - xa/2
    xa2 = x2 + xa/2
    n_a1 = int(n_steps*xa1/la)+1
    n_a2 = int(n_steps*xa2/la)+1
    
    F1 = 25000. # [N] - INITIAL GUESS
    
    error = 100
    correction = 0
    itteration = 0
    itteration_factor = 1 # 0.35 works well
    
    z0 = 0
    a0 = 0
    
    # ITTERATION IN Z DIRECTION #
    while abs(error) >= 10**-10: # [mm]
        itteration += 1
        F1 += correction
        
        shear = [0]*len(x_pos)
        moment = [0]*len(x_pos)
        disp = [0]*len(x_pos)
        angle = [0]*len(x_pos)
        
        disp[0] = z0
        angle[0] = a0
        
        A = np.matrix([ [1,      1,      1         ],
                        [0,      0,      h_arm_a     ],                   
                        [(x2-x1), (x3-x1), (xa1-x1)] ])
        
        B = np.matrix([ [-P-F1],
                        [-P*h_arm_a - q*la*h_arm_q],
                        [-P*(xa2-x1)] ])
        C = np.linalg.solve(A,B)
        
        
        F2 = C[0,0]
        F3 = C[1,0]
        R = C[2,0]
        
        for i in range(1, n_steps+1):
            
            # SHEAR INTEGRATION #
            shear[i] = shear[i-1]                
            if i == n_x1:
                shear[i] = shear[i] + F1
            if i == n_a1:
                shear[i] = shear[i] + R
            if i == n_x2:
                shear[i] = shear[i] + F2
            if i == n_a2:
                shear[i] = shear[i] + P
            if i == n_x3:
                shear[i] = shear[i] + F3
                
            # MOMENT INTEGRATION #    
            moment[i] = moment[i-1] + shear[i]*dx
            
            # ANGLE INTEGRATION #
            angle[i] = angle[i-1] + moment[i]*dx/(E*I1_yy)
            
            # DISPLACEMENT INTEGRATION #
            disp[i] = disp[i-1] + angle[i]*dx
        
        # DISPLACEMENT AND ANGLE ADJUSTMENT #
        x2_adjustment = - disp[n_x2]
        x3_adjustment = - disp[n_x3] + disp[n_x2]
        angle_adjustment = x3_adjustment/(x3-x2)
        for i in range(0, n_steps+1):
            disp[i] = disp[i] + x2_adjustment + x3_adjustment*(x_pos[i]-x2)/(x3-x2)
            angle[i] = angle[i] + angle_adjustment
        
        correction = - F1*(disp[n_x1])*itteration_factor
        z0 = disp[0]
        a0 = angle[0]
        error = disp[n_x1]

    if info:    
        print("Calculated forces in z direction in "+str(itteration)+" itterations, final error: "+str(error)+" mm")
        #print("F1z = "+str(F1)+" R = "+str(R)+" F2z = "+str(F2)+" P = "+str(P)+" F3z = "+str(F3))  
    if plot:
        par = disp
        axis = [0]*len(x_pos)
        h = np.arange(min(par),max(par),(max(par)-min(par))/10) 
        h1x = [x1]*len(h)
        h2x = [x2]*len(h)
        h3x = [x3]*len(h)
        h1a = [xa1]*len(h)
        h2a = [xa2]*len(h)
        
        plt.figure()        
        plt.plot(x_pos, par)
        plt.plot(h1x, h, h2x, h, h3x, h, c="g")
        plt.plot(h1a, h, h2a, h, c="y")
        plt.plot(x_pos, axis, c="r")
        plt.savefig("fig",dpi=1000)
        plt.show()
    if ret:
        return([F1, R, F2, P, F3])

# print(Y_force(plot=True, info=True, ret=True))
#FY = [898063.8378819056, 0, -1560749.1474143378, 0, 677421.7095343104]
#print(Z_force(FY, plot=True, info=True, ret=True))
