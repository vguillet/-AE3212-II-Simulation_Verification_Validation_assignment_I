from src.class_RIB import *
from src.class_SLICE import *

# _________________________________ Define basic model parameters
from config import *

# --------- Auto parameters
slice_interval = la/nb_of_slices

# _________________________________ Gen. initial aileron discretisation
model = []

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

print(model)


# ____________________________________________ Print functions

model[0].print_boom_structure()

# 3D plot of aileron discretisation
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = []
y = []
z = []
for i in range(len(model)):
    if type(model[i]) is Simple_slice:
        for boom in model[i].booms:
            x.append(model[i].x_location)
            y.append(boom.y_location)
            z.append(boom.z_location)

ax.scatter(x, z, y, c='r', marker='o')
ax.set_xlabel('X')
ax.set_ylabel('Z')
ax.set_zlabel('Y')

# plt.gca().set_aspect('equal')
plt.show()
