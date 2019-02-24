from src.class_BOOM import Boom
import numpy as np
from math import *
from config import *

import matplotlib.pyplot as plt


class Slice:
    def __init__(self, label=None, x_location=0):
        # TODO: Add material properties in inputs
        self.label = label
        self.x_location = x_location

        self.x_internal_load = 0
        self.y_internal_load = 0
        self.z_internal_load = 0

        self.x_internal_moment = 0
        self.y_internal_moment = 0
        self.z_internal_moment = 0

        # --- Setup boom structure
        self.construct_boom_structure()
        # --- Calculate centroid
        self.calc_centroid()
        # --- Calculate moment of inertia
        self.calc_mom_inertia()

    def construct_boom_structure(self):
        #  ______________________________ Creating boom structure
        booms = []

        r = ha/2
        self.stiffener_pitch = (pi*r + 2*sqrt(r**2+(ca-r)**2))/nst
        self.angle_tail = atan(r/(ca-r))

        stiffener_count = 0

        # --- Adding slanted edge booms
        boom_z = ca-self.stiffener_pitch/2
        while boom_z > r:
            stiffener_count += 2
            boom_y = (ca - boom_z) * tan(self.angle_tail)
            booms.append(Boom(boom_y, boom_z, "Stringer", stiffener_count-1))
            booms.append(Boom(-boom_y, boom_z, "Stringer", stiffener_count))
            boom_z -= self.stiffener_pitch

        # --- Adding sheet booms
        stiffener_count += 2
        booms.append(Boom(r, r, "Spar", stiffener_count-1))
        booms.append(Boom(-r, r, "Spar", stiffener_count))

        # --- Adding circular section booms
        circ_booms = []
        boom_z = 0
        front_stiffeners = nst + 3

        while stiffener_count != nst+1:
            stiffener_count += 2
            front_stiffeners -= 2
            boom_z += self.stiffener_pitch

            circ_booms.append(Boom(-(r * sin(boom_z / r)), r - (r * cos(boom_z / r)), "Stringer", front_stiffeners))
            circ_booms.append(Boom((r * sin(boom_z / r)), r - (r * cos(boom_z / r)), "Stringer", front_stiffeners-1))

        for boom in reversed(circ_booms):
            booms.append(boom)


        # --- Add leading edge boom
        booms.append(Boom(0, 0, "Stringer", nst+2))

        self.booms = booms

        # ______________________________ Creating variables required
        # --- Calculate spar-stringer skin length (trailing edge side)
        self.booms_z_positions = []
        self.booms_y_positions = []

        for boom in booms:
            self.booms_z_positions.append(boom.z_location)
            self.booms_y_positions.append(boom.y_location)

        for i in range(len(booms)):
            if booms[i].type == "Spar":
                spar_label = booms[i].label
                prev_boom_label = booms[i-2].label

                spar_loc = self.booms_z_positions[i]
                prev_boom_loc = self.booms_z_positions[i-2]

        # print("prev_boom_loc", prev_boom_loc, prev_boom_label)
        # print("spar_loc", spar_loc, spar_label)

        self.small_pitch = (prev_boom_loc-spar_loc)/cos(self.angle_tail)
        self.large_pitch = self.stiffener_pitch-self.small_pitch

        # --- List stringers z locations
        self.z_positions_stringers = []
        self.y_positions_stringers = []

        for boom in booms:
            if boom.type == "Stringer":
                self.z_positions_stringers.append(boom.z_location)
                self.y_positions_stringers.append(boom.y_location)
        return

    def calc_centroid(self):

        # ca = 605.0  # [mm]     airfoil cord
        # ha = 205.0  # [mm]     airfoil height
        # tsk = 1.1  # [mm]      skin thickness
        # tst = 1.2  # [mm]      stringer thickness
        # hst = 16.0  # [mm]     height of stringer
        # wst = 19.0  # [mm]     width of stringer

        # List booms z position
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

        self.z_centroid = z_centroid
        return

    def calc_mom_inertia(self):
        """"Enter lists with the z and y co-ords of the stiffeners and the angle\
        of the skin in the straight section."""

        stif_z = self.z_positions_stringers
        stif_y = self.y_positions_stringers

        angle = self.angle_tail

        ail_h = ha  # height of the aileron     [mm]
        ail_w = ca  # chord of the aileron      [mm]

        # we'll start off by finding moment of inertia for a stiffener
        w = wst  # width of stiffener        [mm]
        h = hst  # height of stiffener       [mm]
        t = tst  # thickness of stiffener    [mm]

        A_base = w * t  # area of the horiz. part   [mm^2]
        A_top = (h - t) * t  # area of the vert. part    [mm^2]
        A_stif = A_base + A_top  # area of stiffener         [mm^2]
        self.area_stiffener = A_stif

        # here are the z co-ordinates of the stiffeners             [mm]
        #    stif_z  = [560.2076385798177, 560.2076385798177, 470.6229157394533, \
        #               470.6229157394533, 381.0381928990888, 381.0381928990888, \
        #               291.45347005872435, 291.45347005872435, 201.86874721835989, \
        #               201.86874721835989, 112.28402437799541, 112.28402437799541, \
        #               37.51424233216027, 37.51424233216027, 0]

        # same in the y-axis                                        [mm]
        #    stif_y  = [8.869774538649953, -8.869774538649953, 26.609323615949847, \
        #               -26.609323615949847, 44.348872693249746, -44.348872693249746, \
        #               62.08842177054964, -62.08842177054964, 79.82797084784953, \
        #               -79.82797084784953, 97.56751992514943, -97.56751992514943, \
        #               78.07387583997614, -78.07387583997614, 0]

        #    angle   = 0.2012196329874848    # angle of skin             [rad]

        c_y = 0  # centroid y co-ord (symm)  [mm]
        c_z = self.z_centroid  # centroid z co-ord         [mm]
        stif_zz = 0  # initial value mom inert   [mm^4]
        stif_yy = 0  # initial value mom inert   [mm^4]

        # now calculate the Steiner terms for all stiffeners
        for i in range(len(stif_z)):
            stein_z = A_stif * (stif_z[i] - c_y) ** 2
            stif_zz = stif_zz + stein_z

        for i in range(len(stif_y)):
            stein_y = A_stif * (stif_y[i] - c_z) ** 2
            stif_yy = stif_yy + stein_y

        # print('---------------------------------------------------------------------')
        # print('Moment of inertia due to the Steiner term of stiffeners:')
        # print('I_zz =', '%e' % stif_zz, 'mm^4, and I_yy =', '%e' % stif_yy, 'mm^4')
        # print('---------------------------------------------------------------------')

        skin_z = (ail_w - (ail_h / 2)) / 2  # skin centroid         [mm]
        skin_t = 1.1  # skin thickness        [mm]
        skin_l = np.sqrt((ail_w - (ail_h / 2)) ** 2 + (ail_h / 2) ** 2)  # skin length           [mm]

        # add the moment of inertia for the straight part of the skin
        skin_zz = ((skin_t * skin_l ** 3) * np.cos(angle) ** 2) / 12
        skin_yy = ((skin_t * skin_l ** 3) * np.sin(angle) ** 2) / 12
        skin_stein = skin_t * skin_l * (skin_z - c_z) ** 2

        # print('Moment of inertia due to the straight part of the skin:')
        # print('I_zz =', '%e' % (skin_zz + skin_stein), 'mm^4, and I_yy =', '%e' % skin_yy, 'mm^4')
        # print('---------------------------------------------------------------------')

        z_out = (ail_h / 2) - (4 * ail_h / 2) / (3 * np.pi)  # centroid outer        [mm]
        z_in = (ail_h / 2) - (4 * (ail_h - skin_t) / 2) / (3 * np.pi)  # centroid inner        [mm]
        A_out = (np.pi * (ail_h / 2) ** 2) / 2  # area outer circ       [mm^2]
        A_in = (np.pi * ((ail_h - skin_t) / 2) ** 2) / 2  # area inner circ       [mm^2]
        circ_z = (A_out * z_out + A_in * z_in) / (A_out + A_in)  # centroid z            [mm]

        # add the moment of inertia for the circular part of the skin
        circ_zz = A_out * (ail_h / 2) ** 2 / 4 - A_in * ((ail_h - skin_t) / 2) ** 2 / 4
        circ_yy = A_out * (ail_h / 2) ** 2 / 4 - A_in * ((ail_h - skin_t) / 2) ** 2 / 4
        circ_stein = (A_in + A_out) * (circ_z - c_z) ** 2

        # print('Moment of inertia due to the circular part of the skin:')
        # print('I_zz =', '%e' % (circ_zz + circ_stein), 'mm^4, and I_yy =', '%e' % circ_yy, 'mm^4')
        # print('---------------------------------------------------------------------')

        spar_z = ail_h / 2  # centroid of the spar      [mm]
        spar_h = ail_h  # height of the spar        [mm]
        spar_t = 2.8  # thickness of spar         [mm]

        # finally, add the moment of inertia for the spar
        spar_zz = (spar_h * spar_t ** 3) / 12
        spar_yy = (spar_t * spar_h ** 3) / 12
        spar_stein = spar_h * spar_t * (spar_z - c_z) ** 2

        # print('Moment of inertia due to the spar:')
        # print('I_zz =', '%e' % (spar_zz + spar_stein), 'mm^4, and I_yy =', '%e' % spar_yy, 'mm^4')
        # print('---------------------------------------------------------------------')

        # then add all the moment of inertia together
        I_zz = stif_zz + skin_zz + skin_stein + circ_zz + circ_stein + spar_zz + spar_stein
        I_yy = stif_yy + skin_yy + circ_yy + spar_yy

        # print('Total moment of inertia:')
        # print('I_zz =', '%e' % I_zz, 'mm^4, and I_yy =', '%e' % I_yy, 'mm^4')
        # print('---------------------------------------------------------------------')

        # now switch to the new axis system
        phi = np.radians(28)  # aileron rot angle         [rad]
        I_u = (I_zz + I_yy) / 2 + ((I_zz - I_yy) / 2) * np.cos(2 * phi)
        I_v = (I_zz + I_yy) / 2 - ((I_zz - I_yy) / 2) * np.cos(2 * phi)

        # print('Total moment of inertia in the rotated axis system:')
        # print('I_u =', '%e' % I_u, 'mm^4, and I_v =', '%e' % I_v, 'mm^4')
        # print('---------------------------------------------------------------------')

        self.I_u = I_zz
        self.I_v = I_yy

        self.I_zz = I_u
        self.I_yy = I_v

        self.polar_I_uv = self.I_u + self.I_v
        self.polar_I_zy = self.I_zz + self.I_yy
        return

    def plot_boom_structure(self):
        label_stringer = []

        z_spar = []
        y_spar = []
        label_spar = []

        for i in range(len(self.booms)):
            if self.booms[i].type == "Stringer":
                label_stringer.append(self.booms[i].label)
            elif self.booms[i].type == "Spar":
                z_spar.append(self.booms[i].z_location)
                y_spar.append(self.booms[i].y_location)
                label_spar.append(self.booms[i].label)

        plt.scatter(self.z_positions_stringers, self.y_positions_stringers, color='b', label="Stringers")
        plt.scatter(z_spar, y_spar, color='r', label="Spar")

        for i, txt in enumerate(label_stringer):
            plt.annotate(txt, (self.z_positions_stringers[i], self.y_positions_stringers[i]))

        for i, txt in enumerate(label_spar):
            plt.annotate(txt, (z_spar[i], y_spar[i]))

        i = 1
        upper_half_z = [ca]
        upper_half_y = [0]

        lower_half_z = [ca]
        lower_half_y = [0]

        while i < len(self.booms_z_positions):
            upper_half_y.append(self.booms_y_positions[i-1])
            upper_half_z.append(self.booms_z_positions[i-1])

            lower_half_y.append(self.booms_y_positions[i])
            lower_half_z.append(self.booms_z_positions[i])

            i += 2

        upper_half_y.append(0)
        upper_half_z.append(0)

        lower_half_y.append(0)
        lower_half_z.append(0)

        plt.plot(upper_half_z, upper_half_y, color="g")
        plt.plot(lower_half_z, lower_half_y, color="g")

        plt.legend()
        plt.grid()
        plt.axis("equal")
        plt.show()


class Simple_slice(Slice):
    def __init__(self, label, x_location):
        Slice.__init__(self, label, x_location)

        self.x_load = 0
        self.y_load = 0
        self.z_load = 0

    def __repr__(self):
        return "Simple slice " + str(self.label)


class Hinged_slice(Slice):
    def __init__(self, label, x_location, x_load=0, y_load=0, z_load=0):
        Slice.__init__(self, label, x_location)

        self.x_load = x_load
        self.y_load = y_load
        self.z_load = z_load

    def __repr__(self):
        return "Hinged slice " + str(self.label)
