#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**json_struc_to_nc.py**

* *Purpose:* Read a json file with a struct and return a ncdf file with this structure

* *python version:* 3.8
* *author:* Pedro Montero
* *license:* INTECMAR
* *requires:*
* *date:* 2021/01/04
* *version:* 0.0.5
* *date version* 2021/01/04

"""
import sys
import json
import numpy

from collections import OrderedDict
import netCDF4


def json_struc_to_nc(input_json):
    """
    Read a json file with a structure and save it in a netcdf file
    :param input_json: name of a json file with the options
    :return:
    """
    try:
        with open(input_json, 'r') as f:
            inputs = json.load(f, object_pairs_hook=OrderedDict)
            nc_file = inputs['nc_file']
            json_file = inputs['json_file']
    except IOError:
        sys.exit('An error occurred trying to read the file.')
    except KeyError:
        sys.exit('An error with a key')
    except ValueError:
        sys.exit('Non-numeric data found in the file.')
    except Exception as err:
        print(err)
        sys.exit(f'Error with the input {input_json}')

    # Read json
    with open(json_file) as f:
        json_data = json.load(f)

    with netCDF4.Dataset(nc_file, "w", format="NETCDF4") as nc:  # TODO: tipo de netcdf debe de ser leido y no impuesto
        for dimension in json_data['dimensions']:

            nc.createDimension(dimension['name'], dimension['size'])
        for gl_att in json_data['global attributes']:

            nc.setncattr(gl_att['name'], gl_att['value'])
        for var in json_data['variables']:
            print('....', var['name'], var['type'], (* var['dimensions'],))
            fill_value = None
            for attr_var in var['attributes']:
                if attr_var['name'] == '_FillValue':
                    fill_value = attr_var['value']

            nc_var = nc.createVariable(var['name'],
                                       var['type'],
                                       (* var['dimensions'],),
                                       fill_value=fill_value)

            for attr_var in var['attributes']:
                print(attr_var['name'], attr_var['value'])
                if attr_var['name'] != '_FillValue':
                    if attr_var['type'] == 'int8':
                        attr_var['value'] = numpy.int8(attr_var['value'])
                    nc_var.setncattr(attr_var['name'], attr_var['value'])






if __name__ == '__main__':
    input_json = 'json_struct_to_nc.json'
    json_struc_to_nc(input_json)
