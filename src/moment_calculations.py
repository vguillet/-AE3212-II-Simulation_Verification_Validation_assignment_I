import numpy as np
import matplotlib.pyplot as plt
from src.class_SLICE import Slice
#get_ipython().run_line_magic('matplotlib','qt')

from config import *

# CONSTANTS #
cross_section = Slice()

theta_max = theta*np.pi/180 # [rad]

hinge = ha/2
act = hinge*np.sqrt(2)

I1_zz = cross_section.I_yy    # [mm**4]
I1_yy = cross_section.I_zz    # [mm**4]
J = I1_yy + I1_zz

sc = 106.12931274441615 # [mm] shear center distance from LE

arm_a = act*np.sin(np.pi/4-theta_max)-(sc - hinge)*np.sin(theta_max) 
arm_q = (ca/4-sc)*np.cos(theta_max) 
arm_y = (sc-hinge)*np.cos(theta_max)
arm_z = (sc-hinge)*np.sin(theta_max)

force_name = ["F1", "R", "F2", "P", "F3"]

def moment_calc(force_y, force_z, n_steps=1000):

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
    
    n_list = [n_x1, n_a1, n_x2, n_a2, n_x3]
    
    theta = [0]*len(x_pos)
    torque = [0]*len(x_pos)
    
    for i in range(1, n_steps+1):
                
        # INTEGRATE MOMENTS TO GET TORQUES#
        torque[i] = torque[i-1] + q*arm_q*dx
        
        if n_list.count(i) == 1:
            if n_list.index(i) == 1 or n_list.index(i) == 3:
                torque[i] = torque[i] + FZ[n_list.index(i)]*arm_a
            else:
                torque[i] = torque[i] + FY[n_list.index(i)]*arm_y - FZ[n_list.index(i)]*arm_z

        #INTEGRATE TORQUES TO GET THETA #
        theta[i] = theta[i-1] - torque[i]*dx/(G*J)
    
    theta_adjustment = theta[n_a1] - theta_max
    theta_deg = theta
    for i in range(0, n_steps+1):
        theta[i] = theta[i] - theta_adjustment
        theta_deg[i] = theta[i]*180/np.pi
    
    par = theta
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
    plt.show()

FZ = [20553.09815723586, -112366.7238326473, 5540.080519864504, 97400.0, -11126.45484445306]
FY = [898063.8378819056, 0, -1560749.1474143378, 0, 677421.7095343104]   
    
#moment_calc(FY,FZ, n_steps=1000)

#test (P+R)*arm_a + sum(FY)*arm_y + q*la*arm_q - z*arm_z
    
    
    
    