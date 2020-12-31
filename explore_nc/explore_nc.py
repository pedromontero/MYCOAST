#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**explore_nc.py**

* *Purpose:* Read a netcdf

* *python version:* 3.8
* *author:* Pedro Montero
* *license:* INTECMAR
* *requires:*
* *date:* 2020/12/07
* *version:* 0.0.5
* *date version* 2020/12/10

"""
import sys
import json
import numpy

from collections import OrderedDict
import netCDF4


class MyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, numpy.integer):
            return int(obj), 'int'
        elif isinstance(obj, numpy.floating):
            return float(obj), 'float'
        elif isinstance(obj, numpy.ndarray):
            return obj.tolist(), 'list'
        else:
            return super(MyEncoder, self).default(obj)


def explore_nc(input_json):
    """
    read a nc and return variables names and attributes
    :param input_json:
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
    global_attr = {}

    with netCDF4.Dataset(nc_file) as nc:
        for gl_att in nc.ncattrs():
            print(gl_att, nc.getncattr(gl_att))
            global_attr[gl_att] = nc.getncattr(gl_att)
            print(global_attr[gl_att])

    json_dict['global attributes'] = global_attr
    with open(json_file, 'w') as file:
        json.dump(json_dict, file, indent=4, cls=MyEncoder)


if __name__ == '__main__':
    input_json = 'explore_nc.json'
    explore_nc(input_json)
