#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**nc_struc_to_json.py**

* *Purpose:* Read a netcdf and return a json file with the structure

* *python version:* 3.8
* *author:* Pedro Montero
* *license:* INTECMAR
* *requires:*
* *date:* 2020/12/07
* *version:* 1.0.0
* *date version* 2021/01/04

"""
import sys
import json
import numpy

from collections import OrderedDict
import netCDF4


def to_json(obj):
    """
    Transform numpy data type to a json compatible datatype
    :param obj: The object whose datatype will be converted
    :return:
    """
    if isinstance(obj, numpy.integer):
        return int(obj)
    elif isinstance(obj, numpy.floating):
        return float(obj)
    elif isinstance(obj, numpy.ndarray):
        return obj.tolist()
    elif isinstance(obj, bytes):
        return str(obj)
    else:
        return obj


def nc_struct_to_json(input_json):
    """
    Read a nc and save its structure in a json file
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

    json_dict = {}

    #  dimensions
    dims = []
    with netCDF4.Dataset(nc_file) as nc:
        for dim in nc.dimensions:
            dic_prov = {'name': nc.dimensions[dim].name, 'size': nc.dimensions[dim].size}
            dims.append(dic_prov)
    json_dict['dimensions'] = dims

    #  global attributes
    global_attr = []
    with netCDF4.Dataset(nc_file) as nc:
        for gl_att in nc.ncattrs():
            attribute = nc.getncattr(gl_att)
            type_attr = type(attribute).__name__
            dic_prov = {'name': gl_att, 'value': to_json(attribute), 'type': type_attr}
            global_attr.append(dic_prov)

    json_dict['global attributes'] = global_attr

    # variables
    variables = []
    with netCDF4.Dataset(nc_file) as nc:
        for var_name in nc.variables:
            var = nc.variables[var_name]

            # dimensions of a variable
            dimensions = []
            for i, dimension in enumerate(var.dimensions):
                dimensions.append({'name': dimension,
                                   'value': var.shape[i]})

            # attributes of a variable
            attributes = []
            for v_attr in var.ncattrs():

                attr = var.getncattr(v_attr)
                attributes.append({'name': v_attr,
                                   'value': to_json(attr),
                                   'type': type(attr).__name__})
            variable = {'name': var_name,
                        'type': str(var.datatype),
                        'dimensions': dimensions,
                        'attributes': attributes}
            variables.append(variable)

    json_dict['variables'] = variables

    #  write json file
    with open(json_file, 'w') as file:
        json.dump(json_dict, file, indent=4)


if __name__ == '__main__':
    input_json = 'nc_struct_to_json.json'
    nc_struct_to_json(input_json)
