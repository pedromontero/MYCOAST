from math import pi, sqrt, atan2, sin, cos


class Vector:
    def __init__(self, wind=False):
        self.wind = wind

    def set_uv(self, u, v):
        self.u = u
        self.v = v

    def get_modtheta(self):
        """
        From u,v velocity components return module and bearing of current

        :param u: x-component of a current
        :param v: y- component of a current
        :param wind: for current or wind(True)
        :return: module, direction
        """
        if self.wind:
            coef = -1
        else:
            coef = 1
        deg2rad = 180. / pi
        mod = sqrt(self.u * self.u + self.v * self.v)
        theta = atan2(coef * self.u, coef * self.v) * deg2rad
        return mod, theta
