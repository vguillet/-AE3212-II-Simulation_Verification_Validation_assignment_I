from config import *


class Rib:
    def __init__(self, label, x_location):
        # TODO: Add material properties in inputs
        self.label = label
        self.x_location = x_location

        self.x_internal_load = 0
        self.y_internal_load = 0
        self.z_internal_load = 0

        self.x_internal_moment = 0
        self.y_internal_moment = 0
        self.z_internal_moment = 0

        self.x_torque = 0
        self.displacement_theta = 0

    def calc_shear(self):
        # TODO: Implement shear calc
        return


class Hinged_rib(Rib):
    def __init__(self, label, x_location, deflection, x_load=0, y_load=0, z_load=0):
        Rib.__init__(self, label, x_location)

        self.deflection = deflection

        self.x_load = x_load
        self.y_load = y_load
        self.z_load = z_load

    def __repr__(self):
        return "Hinged rib " + str(self.label)


class Actuator_rib(Rib):
    def __init__(self, label, x_location, x_load=0, y_load=0, z_load=0):
        Rib.__init__(self, label, x_location)

        self.x_load = x_load
        self.y_load = y_load
        self.z_load = z_load

    def __repr__(self):
        return "Actuator rib " + str(self.label)
