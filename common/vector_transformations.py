#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**vector_transformations.py**

* *Purpose:* module with several transformations for vectors:
* - uv2modtheta, modtheta2uv
* *python version:* 3.8
* *author:* Pedro Montero
* *license:* INTECMAR
* *requires:* math
* *date:* 2021/01/27
* *version:* 1.0.0
* *date version* 2021/01/27

"""

from math import pi, sqrt, atan2, sin, cos


def uv2modtheta(u, v, wind=False):
    """
    From u,v velocity components return module and bearing of current

    :param u: x-component of a current
    :param v: y- component of a current
    :param wind: for current or wind(True)
    :return: module, direction
    """
    if wind:
        coef = -1
    else:
        coef = 1
    deg2rad = 180./pi
    mod = sqrt(u * u + v * v)
    theta = atan2(coef * u, coef * v) * deg2rad
    return mod, theta


def modtheta2uv(mod, theta, wind=False):
    """
    From module, bearing of a current return u,v velocity components

    :param mod: module of a current
    :param theta: angle of a current
    :param wind: for current or wind(True)
    :return: u,v components of the current
    """
    if wind:
        coef = -1
    else:
        coef = 1

    rad2deg = pi / 180.
    u = coef * mod * sin(theta*rad2deg)
    v = coef * mod * cos(theta*rad2deg)
    return u, v
