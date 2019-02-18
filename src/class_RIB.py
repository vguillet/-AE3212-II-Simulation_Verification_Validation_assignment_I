

class Rib:
    def __init__(self, label, x_location):
        # TODO: Add material properties in inputs
        self.label = label
        self.x_location = x_location

    def calc_shear(self):
        # TODO: Implement shear calc
        return


class Hinged_rib(Rib):
    def __init__(self, label, x_location, deflection):
        Rib.__init__(self, label, x_location)

        self.deflection = deflection
        # TODO: Implement y-load based on deflection calc

    def __repr__(self):
        return "Hinged rib " + str(self.label)


class Actuator_rib(Rib):
    def __init__(self, label, x_location, actuator_jammed=False):
        Rib.__init__(self, label, x_location)

        if actuator_jammed:
            self.z_load = 0
        else:
            self.z_load = 97400

    def __repr__(self):
        return "Actuator rib " + str(self.label)