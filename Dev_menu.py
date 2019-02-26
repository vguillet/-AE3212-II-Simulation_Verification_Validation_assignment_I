from src.tools_plot import *

import src.force_calculations_v2 as fc

from src.class_SLICE import *

from config import *
# ============================================ Calculate reaction forces iteratively


def reaction_forces():
    FY = fc.Y_force(n_steps=10000, plot=False, info=False, ret=True)
    FZ = fc.Z_force(FY, n_steps=10000, plot=False, info=False, ret=True)
    return FY, FZ

fy, fz = reaction_forces()
# print("fy", fy)
# print("fz", fz)
fx = [0, 0, 0, 0, 0]

print("\n=========================================================")
print("                 Reaction forces found!")
print("=========================================================\n")

# ============================================ Define basic model generation parameters
slice_interval = la/(nb_of_slices-1)

# _________________________________ Gen. initial aileron discretisation
model = list()

# --------- Adding ribs and hinged slice to model and setting actuator status
model.append(Rib("A", x1, deflection=d1, x_load=fx[0], y_load=fy[0], z_load=fz[0]))
model.append(Rib("B", x2-(xa/2), x_load=fx[1], y_load=fy[1], z_load=fz[1]))

model.append(Hinged_slice(1, x2, x_load=fx[2], y_load=fy[2], z_load=fz[2]))

model.append(Rib("C", x2+(xa/2), x_load=fx[3], y_load=fy[3], z_load=fz[3]))
model.append(Rib("D", x3, deflection=d3, x_load=fx[4], y_load=fy[4], z_load=fz[4]))

# --------- Adding simple slices to model
loc = 0
for i in range(nb_of_slices):
    model.append(Simple_slice(i+1, loc))
    loc += slice_interval

# --------- Sort model points according to x position:
n_times = len(model)
# Traverse all list elements
for i in range(n_times):
    for j in range(0, n_times-i-1):
        # traverse the list from 0 to n-i-1 (e.g number of possible swap left)
        if model[j].x_location > model[j+1].x_location:
            model[j], model[j+1] = model[j+1], model[j]

# --------- List all x positions
x_lst = []
for i in range(len(model)):
    x_lst.append(model[i].x_location)

# _________________________________ Calculate internal forces in every slice
features = []
feature_index = []
x_separation = []

count = 0
for i in range(len(model)-1):
    count += 1
    x_separation.append(model[i+1].x_location-model[i].x_location)

for i in range(len(model)):
    if type(model[i]) is not Simple_slice:
        features.append(model[i])
        feature_index.append(i)


print(features)

# --- Calc. y/z internal forces and z internal moment
for i in range(len(model)-1):
    model[i+1].y_internal_load = -(-model[i].y_internal_load + model[i+1].y_load + -q*(model[i+1].x_location - model[i].x_location))
    model[i+1].z_internal_load = -(-model[i].z_internal_load + model[i+1].z_load)

    fy1 = 0
    fz1 = 0

    fy2 = 0
    fz2 = 0

    fy3 = 0
    fz3 = 0

    fy4 = 0
    fz4 = 0

    fy5 = 0
    fz5 = 0

    if model[i+1].x_location > features[0].x_location:
        fy1 = features[0].y_load
        fz1 = features[0].z_load

    if model[i+1].x_location > features[1].x_location:
        fy2 = features[1].y_load
        fz2 = features[1].z_load

    if model[i+1].x_location > features[2].x_location:
        fy3 = features[2].y_load
        fz3 = features[2].z_load

    if model[i+1].x_location > features[3].x_location:
        fy4 = features[3].y_load
        fz4 = features[3].z_load

    if model[i+1].x_location > features[4].x_location:
        fy5 = features[4].y_load
        fz5 = features[4].z_load

    model[i+1].y_internal_moment = -(fy1*(model[i+1].x_location - features[0].x_location)
                                     + fy2 * (model[i + 1].x_location - features[1].x_location)
                                     + fy3 * (model[i + 1].x_location - features[2].x_location)
                                     + fy4 * (model[i + 1].x_location - features[3].x_location)
                                     + fy5 * (model[i + 1].x_location - features[4].x_location)
                                     - (q*model[i+1].x_location)*(model[i+1].x_location/2))

    model[i+1].z_internal_moment = -(fz1*(model[i+1].x_location - features[0].x_location)
                                     + fz2 * (model[i + 1].x_location - features[1].x_location)
                                     + fz3 * (model[i + 1].x_location - features[2].x_location)
                                     + fz4 * (model[i + 1].x_location - features[3].x_location)
                                     + fz5 * (model[i + 1].x_location - features[4].x_location))

# --- Calc. internal torque
from src.moment_calculations import moment_calc

torque_lst, theta_lst = moment_calc(fy, fz, x_separation, feature_index)

for i in range(len(model)-1):
    model[i+1].x_torque = torque_lst[i]
    model[i+1].displacement_theta = theta_lst[i]


# --- List results
f_internal_y = []
f_internal_z = []
m_internal_y = []
m_internal_z = []
x_torque = []

for slice_point in model:
    f_internal_y.append(slice_point.y_internal_load)
    f_internal_z.append(slice_point.z_internal_load)

    m_internal_y.append(slice_point.y_internal_moment)
    m_internal_z.append(slice_point.z_internal_moment)

    x_torque.append(slice_point.x_torque)

# ============================================ Calc. internal shear
for i in range(len(model)-1):
    if type(model[i+1]) is Rib:
        model[i+1].calc_shear_flow()

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

print("\nShear center u ("+length_unit+"^4):", test_slice.shear_center_u)


# plot_internal_forces(f_internal_y, "y load")
# plot_internal_forces(f_internal_z, "z load")

plot_internal_forces(m_internal_y, "y moment")
plot_internal_forces(m_internal_z, "z moment")

# plot_internal_forces(x_torque, "Torque")
# plot_internal_forces(theta_lst, "Theta")

# test_slice.plot_boom_structure()

# --------- Plot aileron model structure
# plot_3d_aileron(model)

# --------- Plot leading/trailing edge deflections
# plot_le_te_deflection(le_deflection, te_deflection, x_lst)




