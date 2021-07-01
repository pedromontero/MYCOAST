from math import pi, atan2, sin, cos
import numpy as np


class Vector:

    def __init__(self, u=None, v=None, mod=None, theta=None, wind=False):
        self.wind = wind

        self.u = u
        self.v = v
        self.mod = mod
        self.theta = theta

        self.x = None
        self.y = None

        if self.mod is None or self.theta is None:
            self.mod, self.theta = self.get_modtheta()

        if self.u is None or self.v is None:
            self.u, self.v = self.get_uv()

    def set_point(self, x, y):
        self.x = x
        self.y = y

    def get_modtheta(self):
        """
        From u,v velocity components return module and bearing of current

        """
        if self.wind:
            coef = -1
        else:
            coef = 1
        deg2rad = 180. / pi

        mod = pow((pow(self.u, 2) + pow(self.v, 2)), .5)
        theta = atan2(coef * self.u, coef * self.v) * deg2rad
        return mod, theta

    def get_uv(self):
        """
        From module, bearing of a current return u,v velocity components

        """
        if self.wind:
            coef = -1
        else:
            coef = 1

        rad2deg = pi / 180.
        u = coef * self.mod * sin(self.theta * rad2deg)
        v = coef * self.mod * cos(self.theta * rad2deg)
        return u, v


class NumpyVector(Vector):

    def get_modtheta(self):
        """
        From u,v velocity components return module and bearing of current

        """
        if self.wind:
            coef = -1
        else:
            coef = 1
        deg2rad = 180. / pi

        mod = pow((pow(self.u, 2) + pow(self.v, 2)), .5)
        theta = np.arctan2(coef * self.u, coef * self.v) * deg2rad
        return mod, theta
