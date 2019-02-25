from src.tools_plot import *

import src.force_calculations_v2 as fc
import src.moment_calculations as mc

from src.class_RIB import *
from src.class_SLICE import *

from config import *
# _________________________________ Calculate reaction forces iteratively


def reaction_forces():
    FY = fc.Y_force(n_steps=10000, plot=False, info=False, ret=True)
    FZ = fc.Z_force(FY, n_steps=10000, plot=False, info=False, ret=True)
    return FY, FZ


fy, fz = reaction_forces()

fx = [0, 0, 0, 0, 0]

print("\n=========================================================")
print("                 Reaction forces found!")
print("=========================================================\n")

# _________________________________ Define basic model generation parameters
slice_interval = la/nb_of_slices

# _________________________________ Gen. initial aileron discretisation
model = list()

# --------- Adding ribs and hinged slice to model and setting actuator status
model.append(Hinged_rib("A", x1, d1, x_load=fx[0], y_load=fy[0], z_load=fz[0]))
model.append(Actuator_rib("B", x2-(xa/2), x_load=fx[1], y_load=fy[1], z_load=fz[1]))

model.append(Hinged_slice(1, x2, x_load=fx[2], y_load=fy[2], z_load=fz[2]))

model.append(Actuator_rib("C", x2+(xa/2), x_load=fx[3], y_load=fy[3], z_load=fz[3]))
model.append(Hinged_rib("D", x3, d3, x_load=fx[4], y_load=fy[4], z_load=fz[4]))

# --------- Adding simple slices to model
loc = -slice_interval
for i in range(nb_of_slices):
    loc += slice_interval
    model.append(Simple_slice(i+1, loc))

# --------- Sort model points according to x position:
n_times = len(model)
# Traverse all list elements
for i in range(n_times):
    for j in range(0, n_times-i-1):
        # traverse the list from 0 to n-i-1 (e.g number of possible swap left)
        if model[j].x_location > model[j+1].x_location:
            model[j], model[j+1] = model[j+1], model[j]

# --------- List x positions
x_lst = []
for i in range(len(model)):
    x_lst.append(model[i].x_location)

# _________________________________ Calculate internal forces in every slice
# --- Calc. y internal forces
for i in range(len(model)-1):
    # print(model[i].y_load)
    model[i+1].y_internal_load = -(-model[i].y_internal_load + model[i+1].y_load + -q*(model[i+1].x_location - model[i].x_location))

# --- Calc. z internal forces
for i in range(len(model)-1):
    # print(model[i].z_load)
    model[i+1].z_internal_load = -(-model[i].z_internal_load + model[i+1].z_load)
    
# --- Calc. internal torque
from src.moment_calculations import moment_calc

feature_index = []
x_separation = []

for i in range(len(model)-1):
    x_separation.append(model[i+1].x_location-model[i].x_location)

for i in range(len(model)):
    if type(model[i]) is not Simple_slice:
        feature_index.append(i)

torque_lst, theta_lst = moment_calc(fy, fz, x_separation, feature_index)




# --- List results
internal_y = []
internal_z = []

for slice_point in model:
    internal_y.append(slice_point.y_internal_load)
    internal_z.append(slice_point.z_internal_load)


# ============================================ Print/Plot functions
# --------- Print model layout
print("Model structure:", model)

print("\n---------------------------------------------------------\n")

test_slice = Slice("Test Slice", None)

# --------- Plot aileron cross-section and print positions of booms
print("The following excludes Spar positions:")
print("Stringer z-position ("+length_unit+"):", test_slice.z_positions_stringers)
print("Stringer y-position ("+length_unit+"):", test_slice.y_positions_stringers)

print("\nThe following contains all booms location (including spar booms)")
print("Booms z-position ("+length_unit+"):", test_slice.booms_z_positions)
print("Booms y-position ("+length_unit+"):", test_slice.booms_y_positions)

print("\n---------------------------------------------------------\n")

print("~~~~~~~~~~~~~~Physical properties:~~~~~~~~~~~~~~")
print("Stiffener pitch ("+length_unit+"):", test_slice.stiffener_pitch)
print("Small pitch ("+length_unit+"):", test_slice.small_pitch)
print("Large pitch ("+length_unit+"):", test_slice.large_pitch)

print("\nTrailing edge angle ("+angle_unit+"):", test_slice.angle_tail)
print("Trailing edge angle (deg):", test_slice.angle_tail*180/pi)

print("\n~~~~~~~~~~~~~~Geometry properties:~~~~~~~~~~~~~~")
print("Centroid of section ("+length_unit+"):", test_slice.u_centroid)
print("Area stiffener ("+length_unit+"^2):", test_slice.area_stiffener)

print("\nMoment of inertia (reference axis system):")
print("Moment of inertia z ("+length_unit+"^4):", "%e" % test_slice.I_zz)
print("Moment of inertia y ("+length_unit+"^4):", "%e" % test_slice.I_yy)

print("\nMoment of inertia (aileron axis system):")
print("Moment of inertia u ("+length_unit+"^4):",  "%e" % test_slice.I_u)
print("Moment of inertia v ("+length_unit+"^4):",  "%e" % test_slice.I_v)

print("\nPolar moment of inertia:")
print("Polar moment of inertia J ("+length_unit+"^4):",  "%e" % test_slice.polar_I_zy)

print("\nShear center u:", test_slice.shear_center_u)


# plot_internal_forces(internal_y, "y")
# plot_internal_forces(internal_z, "z")
# plot_internal_forces(torque_lst, "Torque")
# plot_internal_forces(theta_lst, "Theta")
# test_slice.plot_boom_structure()

# --------- Plot aileron model structure
# plot_3d_aileron(model)

# --------- Plot leading/trailing edge deflections
# plot_le_te_deflection(le_deflection, te_deflection, x_lst)




