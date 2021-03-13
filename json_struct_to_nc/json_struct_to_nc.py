#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**json_struct_to_nc.py**

* *Purpose:* Read a json file with a struct and return a ncdf file with this structure

* *python version:* 3.8
* *author:* Pedro Montero
* *license:* INTECMAR
* *requires:*
* *date:* 2021/01/04
* *version:* 0.9.0
* *date version* 2021/01/07

"""
import sys
import json
import numpy

from collections import OrderedDict
import netCDF4


def json_struct_to_nc(input_json):
    """
    Read a json file with a structure and save it in a netcdf file
    :param input_json: name of a json file with the options
    :return: ncfile, name of netCDF
    """
    try:
        with open(input_json, 'r') as f:
            inputs = json.load(f, object_pairs_hook=OrderedDict)
            nc_file = inputs['nc_file']
            json_file = inputs['json_struct_file']
    except IOError:
        sys.exit('json_struct_to_nc: An error occurred trying to read the file.')
    except KeyError:
        sys.exit('json_struct_to_nc: An error with a key')
    except ValueError:
        sys.exit('json_struct_to_nc: Non-numeric data found in the file.')
    except Exception as err:
        print(err)
        sys.exit(f'json_struct_to_nc: Error with the input {input_json}')

    # Read json
    with open(json_file) as f:
        json_data = json.load(f)

    # Write netCDF
    with netCDF4.Dataset(nc_file, "w", format="NETCDF4") as nc:  # TODO: tipo de netcdf debe de ser leido y no impuesto
        # Dimensions
        for dimension in json_data['dimensions']:
            nc.createDimension(dimension['name'], dimension['size'])
        # Global attributes
        for gl_att in json_data['global attributes']:
            nc.setncattr(gl_att['name'], gl_att['value'])
        # Variables
        for var in json_data['variables']:
            fill_value = var['fill_value']
            nc_var = nc.createVariable(var['name'],
                                       var['type'],
                                       (* var['dimensions'],),
                                       fill_value=fill_value)
            # Attributes of a variable
            for attr_var in var['attributes']:
                if attr_var['type'] == 'ndarray':
                    attr_var['value'] = numpy.dtype(attr_var['type_element']).type(attr_var['value'])
                else:
                    attr_var['value'] = numpy.dtype(attr_var['type']).type(attr_var['value'])
                nc_var.setncattr(attr_var['name'], attr_var['value'])

            # Value
            if 'value' in var:
                #nc_var[:] = var['value']*numpy.ones(nc_var.shape, dtype=var['type'])
                nc_var[:] = numpy.full(nc_var.shape, var['value'], dtype=var['type'])


    return nc_file


if __name__ == '__main__':
    input_json = 'json_struct_to_nc.json'
    json_struct_to_nc(input_json)
