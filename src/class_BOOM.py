

class Boom:
    def __init__(self, y_location, z_location, type, label=None):
        self.label = label
        self.type = type

        self.y_location = y_location
        self.z_location = z_location

    def calc_boom_area(self):
        # TODO: Implement boom area/stresses based on force and MOI calc.
        return
