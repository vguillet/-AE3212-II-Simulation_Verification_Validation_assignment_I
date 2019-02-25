import numpy as np
import matplotlib.pyplot as plt
from src.class_SLICE import Slice
#get_ipython().run_line_magic('matplotlib','qt')

from src.tools_plot import *

from config import *


def moment_calc(FY, FZ, x_pos, feature_lst):

    # INITIATION #
    n_list = feature_lst
    n_steps = len(x_pos)

    # CONSTANTS #
    cross_section = Slice()

    theta_max = theta*np.pi/180     # [rad]
    hinge = ha/2
    act = hinge*np.sqrt(2)

    I1_zz = cross_section.I_yy    # [mm**4]
    I1_yy = cross_section.I_zz    # [mm**4]
    J = I1_yy + I1_zz
    sc = cross_section.shear_center_u # [mm] shear center distance from LE

    arm_a = act*np.sin(np.pi/4-theta_max)-(sc - hinge)*np.sin(theta_max)
    arm_q = (ca/4-sc)*np.cos(theta_max)
    arm_y = (sc-hinge)*np.cos(theta_max)
    arm_z = (sc-hinge)*np.sin(theta_max)

    n_a1 = n_list[1]

    theta_lst = [0]*len(x_pos)
    torque = [0]*len(x_pos)

    for i in range(1, n_steps):
        dx = x_pos[i-1]

        # INTEGRATE MOMENTS TO GET TORQUES#
        torque[i] = torque[i-1] + q*arm_q*dx

        if n_list.count(i) == 1:
            if n_list.index(i) == 1 or n_list.index(i) == 3:
                torque[i] = torque[i] + FZ[n_list.index(i)]*arm_a
            else:
                torque[i] = torque[i] + FY[n_list.index(i)]*arm_y - FZ[n_list.index(i)]*arm_z

        # INTEGRATE TORQUES TO GET THETA #
        theta_lst[i] = theta_lst[i-1] - torque[i]*dx/(G*J)

    theta_adjustment = theta_lst[n_a1] + theta_max
    theta_deg = theta_lst
    for i in range(0, n_steps):
        theta_lst[i] = theta_lst[i] - theta_adjustment
        theta_deg[i] = -theta_lst[i]*180/np.pi

    return torque, theta_deg
