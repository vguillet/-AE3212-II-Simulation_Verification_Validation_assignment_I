import matplotlib.pyplot as plt
from src.class_SLICE import *
from src.class_RIB import *


def plot_3d_aileron(model):
    # 3D plot of aileron discretisation
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    x_ss = []
    y_ss = []
    z_ss = []
    for i in range(len(model)):
        if type(model[i]) is Simple_slice:
            for boom in model[i].booms:
                x_ss.append(model[i].x_location)
                y_ss.append(boom.y_location)
                z_ss.append(boom.z_location)

    ax.scatter(x_ss, z_ss, y_ss, c='r', marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Z')
    ax.set_zlabel('Y')

    plt.grid()

    # plt.gca().set_aspect('equal')
    plt.show()


def plot_le_te_deflection(le_deflection, te_deflection, x_position):
    plt.plot(x_position, le_deflection, label="Leading edge deflection")
    plt.plot(x_position, te_deflection, label="Trailing edge deflection")

    plt.title("Leading and trailing edge deflection over the aileron span")
    plt.xlabel("x position (m)")
    plt.ylabel("Deflection (m)")

    plt.legend()
    plt.grid()

    plt.show()
