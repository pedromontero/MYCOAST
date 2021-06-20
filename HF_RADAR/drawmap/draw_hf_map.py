#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
draw_hf_map.py

@Purpose: draw the current map from HF radar (nc file)
@version: 0.1

@python version: 3.4
@author: bvila
@license: INTECMAR
@requires: netCDF4, matplotlib, numpy, toolkits.basemap

@date 2015/10/09

@history:

"""


import netCDF4
import numpy as np
import datetime


from math import  pi

from BoundaryBox import BoundaryBox


def drawcurrents(lats, lons, ust, vst, mod, total, title, time, boundary_box):
    """

    :return:
    """

    import matplotlib.pyplot as plt
    import matplotlib.spines as spn
    from mpl_toolkits.basemap import Basemap

    fig = plt.figure(figsize=(15, 15))
    ax = fig.add_subplot(1, 1, 1)
    for child in ax.get_children():
        if isinstance(child, spn.Spine):
            child.set_color('#eeeeee')

    dx = 0.1
    middle_lon = boundary_box.middle_lon()
    middle_lat = boundary_box.middle_lat()
    m = Basemap(llcrnrlon=boundary_box.lon_min-dx,
                llcrnrlat=boundary_box.lat_min-dx,
                urcrnrlon=boundary_box.lon_max+dx,
                urcrnrlat=boundary_box.lat_max+dx,
                resolution='i', projection='tmerc', lon_0=middle_lon, lat_0=middle_lat)

    m.drawcoastlines()
    m.fillcontinents(color='grey', lake_color='aqua')

    if total:
        lon, lat = np.meshgrid(lons, lats)
        x, y = m(lon, lat)
    else:
        x, y = m(lons, lats)


    # draw coloured vectors.
    cs = m.quiver(x, y, ust, vst, mod, clim=[0, 0.7], scale=5)

    # add colorbar.
    cbar = m.colorbar(cs, location='bottom', pad="5%")
    cbar.ax.tick_params(labelsize=9)
    cbar.set_label('sea water velocity (m/s)', fontsize=9)

    name_fig = 'imaxe1.png'
    fig.savefig(name_fig, dpi=300, facecolor='w', edgecolor='w', format='png',
                transparent=False, pad_inches=0.1, bbox_inches='tight')

    plt.clf()
    plt.close('all')

    return




def unix_time(dt):
    """ Seconds since 01_01_1970."""
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = dt - epoch
    return delta.total_seconds()


def getvar_standardname(f, nome_standards):
    """Return values using the CF standard name of a variable in a netCDF file."""
    for var in f.variables:
        for atributo in (f.variables[var].ncattrs()):
            if atributo == 'standard_name':
                nome_atributo = (getattr(f.variables[var], 'standard_name'))
                for nome_standar in nome_standards:
                    if nome_atributo == nome_standar:
                        return f.variables[var]
    print('standard_name = {0} not found'.format(nome_standar))


def getvar_longname(f, nome_longs):
    """Return values using the CF long name of a variable in a netCDF file."""
    for var in f.variables:
        for atributo in (f.variables[var].ncattrs()):
            if atributo == 'long_name':
                nome_atributo = (getattr(f.variables[var], 'long_name'))
                for nome_long in nome_longs:
                    if nome_atributo == nome_long:
                        return f.variables[var]
    print('long_name = {0} not found'.format(nome_long))






def main():
    #file_in = r'http://150.145.136.27:8080/thredds/dodsC/Ibiza_NRT212/2020/2020_02/2020_02_12/HFR-Ibiza-Total_2020_02_12_1700.nc'
    file_in = '../codar2nc/data/HFR-Galicia-VILA_2021_04_26_0600.nc'
    print('vou a ler {0}'.format(file_in))

    f = netCDF4.Dataset(file_in)

    nc_attrs = f.ncattrs()

    detailed_text = 'NetCDF Global Attributes:\n\n'
    for nc_attr in nc_attrs:
        value = '%s' % repr(f.getncattr(nc_attr), )
        spam = f'- {nc_attr}: {value}; \n'
        detailed_text += spam

    # Radial or Total file
    total = False
    if f.getncattr('data_type') == 'HF radar total data':
        total = True

    # Extension

    boundary_box = BoundaryBox()
    boundary_box.lat_min = float(f.getncattr('geospatial_lat_min'))
    boundary_box.lat_max = float(f.getncattr('geospatial_lat_max'))
    boundary_box.lon_min = float(f.getncattr('geospatial_lon_min'))
    boundary_box.lon_max = float(f.getncattr('geospatial_lon_max'))

    # Variables with time

    times_in = getvar_standardname(f, ['time'])
    tempos = netCDF4.num2date(times_in[:], units=times_in.units)
    lat_in = getvar_standardname(f, ['latitude'])[:]
    lon_in = getvar_standardname(f, ['longitude'])[:]

    u_in = getvar_standardname(f, ['surface_eastward_sea_water_velocity',
                                   'eastward_sea_water_velocity'])[:]
    # if u_in is None:
    # u_in = getvar_standardname(f, 'eastward_sea_water_velocity')[:]
    v_in = getvar_standardname(f, ['surface_northward_sea_water_velocity',
                                   'northward_sea_water_velocity'])[:]
    # if v_in is None:
    # v_in = getvar_standardname(f, 'northward_sea_water_velocity')[:]
    print(v_in.shape)
    print(lat_in.shape)

    tempo = netCDF4.num2date(times_in[:], units=times_in.units)[0]
    times = unix_time(tempo)
    print(tempo)

    water_u = u_in[:]
    water_v = v_in[:]
    water_u = water_u[0][0]
    water_v = water_v[0][0]
    print(water_u.size, water_v.size, lon_in.size)

    mod = pow((pow(water_u, 2) + pow(water_v, 2)), .5)
    dir = (180 * np.arctan2(water_u, water_v)) / pi

    f.close()

    drawcurrents(lat_in, lon_in, water_u, water_v, mod, total, 'title', tempo, boundary_box)



if __name__ == '__main__':
    main()
