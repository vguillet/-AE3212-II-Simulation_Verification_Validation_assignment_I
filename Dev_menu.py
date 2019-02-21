from src.tools_plot import *
from src.class_RIB import *
from src.class_SLICE import *

# _________________________________ Define basic model parameters
from config import *

# --------- Auto parameters
slice_interval = la/nb_of_slices

# _________________________________ Gen. initial aileron discretisation
model = list()

# --------- Adding ribs to model and setting actuator status
model.append(Hinged_rib("A", x1, d1))
model.append(Actuator_rib("B", x2-(xa/2), True))
model.append(Actuator_rib("C", x2+(xa/2), False))
model.append(Hinged_rib("D", x3, d3))

# --------- Adding hinged slice
model.append(Hinged_slice(1, x2))

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

# _________________________________ Calculate reaction forces interatively

# ============================================ Print/Plot functions
# --------- Print model layout
print(model)

test_slice = Slice("Test Slice", None)

# --------- Plot aileron cross-section and print positions of booms
print("\nThe following excludes Spar positions:")
print("Stringer z-position ("+length_unit+"):", test_slice.z_positions_stringers)
print("Stringer y-position ("+length_unit+"):", test_slice.y_positions_stringers)

print("\nThe following contains all booms location (including spar booms)")
print("Booms z-position ("+length_unit+"):", test_slice.booms_z_positions)
print("Booms y-position ("+length_unit+"):", test_slice.booms_y_positions)

print("\nStiffener pitch ("+length_unit+"):", test_slice.stiffener_pitch)
print("Small pitch ("+length_unit+"):", test_slice.small_pitch)
print("Large pitch ("+length_unit+"):", test_slice.large_pitch)

print("\nTrailing edge angle ("+angle_unit+"):", test_slice.angle_tail)
print("Trailing edge angle (deg):", test_slice.angle_tail*180/pi)

print("\nCentroid of section ("+length_unit+"):", test_slice.z_centroid)
print("Area stiffener ("+length_unit+"^2):", test_slice.area_stiffener)

print("\nMoment of inertia (reference axis system):")
print("Moment of inertia z ("+length_unit+"^4):", "%e" % test_slice.I_zz)
print("Moment of inertia y ("+length_unit+"^4):", "%e" %  test_slice.I_yy)

print("\nMoment of inertia (aileron axis system):")
print("Moment of inertia u ("+length_unit+"^4):",  "%e" % test_slice.I_u)
print("Moment of inertia v ("+length_unit+"^4):",  "%e" % test_slice.I_v)

print("\nPolar moment of inertia:")
print("Polar moment of inertia J ("+length_unit+"^4):",  "%e" % test_slice.polar_I_zy)


test_slice.plot_boom_structure()

# --------- Plot aileron model structure
# plot_3d_aileron(model)

# --------- Plot leading/trailing edge deflections
# plot_le_te_deflection(le_deflection, te_deflection, x_lst)
