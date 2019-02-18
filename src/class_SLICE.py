from src.class_BOOM import Boom
import numpy as np
from math import *
from config import *

import matplotlib.pyplot as plt


class Slice:
    def __init__(self, label, x_location):
        # TODO: Add material properties in inputs
        self.label = label
        self.x_location = x_location

        # --- Setup boom structure
        self.booms = self.construct_boom_structure()
        # --- Calculate centroid
        self.z_centroid = self.calc_centroid()

        # print(self.z_centroid)

    @staticmethod
    def construct_boom_structure():
        booms = []

        r = ha/2
        stiffener_pitch = (pi*r + 2*sqrt(r**2+(ca-r)**2))/nst
        angle_tail = atan(r/(ca-r))

        stiffener_count = 0

        # Adding slanted edge booms
        boom_z = ca-stiffener_pitch/2
        while boom_z > r:
            stiffener_count += 2
            boom_y = (ca - boom_z) * tan(angle_tail)
            booms.append(Boom(boom_y, boom_z, "Stringer", stiffener_count-1))
            booms.append(Boom(-boom_y, boom_z, "Stringer", stiffener_count))
            boom_z -= stiffener_pitch

        # Adding sheet booms
        stiffener_count += 2
        booms.append(Boom(r, r, "Spar", stiffener_count-1))
        booms.append(Boom(-r, r, "Spar", stiffener_count))

        # Adding circular section booms
        boom_z = 0
        while stiffener_count != nst+1:
            stiffener_count += 2
            boom_z += stiffener_pitch

            booms.append(Boom((r * sin(boom_z / r)), r - (r * cos(boom_z / r)), "Stringer", stiffener_count-1))
            booms.append(Boom(-(r * sin(boom_z / r)), r - (r * cos(boom_z / r)), "Stringer", stiffener_count))

        # Add leading edge boom
        booms.append(Boom(0, 0, "Stringer", nst+2))

        return booms

    def calc_centroid(self):

        # ca = 605.0  # [mm]     airfoil cord
        # ha = 205.0  # [mm]     airfoil height
        # tsk = 1.1  # [mm]      skin thickness
        # tst = 1.2  # [mm]      stringer thickness
        # hst = 16.0  # [mm]     height of stringer
        # wst = 19.0  # [mm]     width of stringer

        # List z boom positions
        z_pos = []
        for boom in self.booms:
            if boom.type == "Stringer":
                z_pos.append(boom.z_location)
            elif boom.type == "Spar":
                z_spar = boom.z_location

        a_st = wst * tst + (hst - tst) * tst  # area of one stiffener

        # Areas and centroids of parts of cross section from leading edge
        a_circ = np.pi * ha * tsk / 2  # area of circular section
        z_circ = ha / 2 * (1 - 4 / (3 * np.pi))  # centroid of circular section
        a_straight = 2 * tsk * np.sqrt((ca - ha / 2) ** 2 + (ha / 2) ** 2)  # area of both straight parts
        z_straight = ha / 2 + (ca - ha / 2) / 2  # centroid of both straight parts
        a_stiffeners = a_st * len(z_pos)  # area of all stiffeners
        z_stiffeners = sum(z_pos) / len(z_pos)  # centroid of all stiffeners
        a_spar = ha * tsk

        z_centroid = (z_circ * a_circ + z_straight * a_straight + z_stiffeners * a_stiffeners + a_spar * z_spar) / \
                     (a_circ + a_straight + a_stiffeners + a_spar)

        return z_centroid

    def print_boom_structure(self):
        z = []
        y = []
        label = []

        for i in range(len(self.booms)):
            z.append(self.booms[i].z_location)
            y.append(self.booms[i].y_location)
            label.append(self.booms[i].label)

        plt.scatter(z, y)
        plt.axis("equal")
        for i, txt in enumerate(label):
            plt.annotate(txt, (z[i], y[i]))
        plt.show()


class Simple_slice(Slice):
    def __init__(self, label, x_location):
        Slice.__init__(self, label, x_location)

        self.x_load = 0
        self.y_load = None
        self.z_load = None

    def __repr__(self):
        return "Simple slice " + str(self.label)


class Hinged_slice(Slice):
    def __init__(self, label, x_location):
        Slice.__init__(self, label, x_location)

        self.x_load = 0
        self.y_load = None
        self.z_load = None

    def __repr__(self):
        return "Hinged slice " + str(self.label)
