#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**platform.py**

* *Purpose:* Read data from a platform and a json structure file and create a common netcdf

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
from collections import OrderedDict
from datetime import datetime, timedelta
import numpy
import pandas as pd
import netCDF4

from db import consulta_pd
from WP6.json_struct_to_nc.json_struct_to_nc import json_struct_to_nc

from common import uv2modtheta


def date_range(start_date, end_date, dt):
    """
    Return a list of dates >= start_date and < end_date and dt =dt in seconds
    :param start_date:datetime, start date
    :param end_date: datetime, end date
    :param dt: int, the time step in seconds
    :return: a list the dates
    """
    dates = []
    delta = timedelta(seconds=dt)
    while start_date < end_date:
        dates.append(start_date)
        start_date += delta
    return dates


def read_data(db_json_file, start, end, ln_station,  ln_parameter):
    """
    Make  a query to the Ocean Platform database of a parameter of one station between two dates
     (>= start and < end) and return the result in a pandas dataframe

    :param db_json_file: str, a json file name with the DB connection parameters
    :param start: datetime, the start date
    :param end: datetime, the end date
    :param ln_station: int, index of the station
    :param ln_parameter: int, index of the parameter
    :return: a pandas dataframe with the time, value of the parameter and its QC flag
    """

    try:
        with open(db_json_file, 'r') as f:
            con_parameters = json.load(f, object_pairs_hook=OrderedDict)
    except IOError:
        sys.exit(f'read_data: An error happened trying to read the file {db_json_file}')
    except KeyError:
        sys.exit('read_data:An error with a key')
    except ValueError:
        sys.exit('read_data:Non-numeric data found in the file.')
    except Exception as err:
        print(err)
        sys.exit(f'read_data: Error with the input {db_json_file}')

    str_start = start.strftime("%Y-%m-%d %H:%M:%S")
    str_end = end.strftime("%Y-%m-%d %H:%M:%S")

    cons = "SELECT TOP (1000) [InstanteLectura] ,[Valor] ,[LnCodigoValidacion] " \
           "FROM [tHstUltimosValidados] INNER JOIN [tCfgCnxEstaciones] " \
           "ON  [tHstUltimosValidados].lnEstacion = [tCfgCnxEstaciones].IdEstacion " \
           "WHERE lnEstacion = ? " \
           "AND InstanteLectura >= CAST(? AS date) " \
           "AND InstanteLectura < CAST(? AS date) " \
           "AND lnParametro = ?"
    params = (ln_station, str_start, str_end, ln_parameter,)
    return consulta_pd(con_parameters, cons, params)

def platform_nc(input_json_file):
    """

    :param input_json_file:
    :return:
    """

    try:
        with open(input_json_file, 'r') as f:
            inputs = json.load(f, object_pairs_hook=OrderedDict)
            db_json_file = inputs['db_con']
            start_date = datetime.strptime(inputs['start'], "%Y-%m-%dT%H:%M:%S")
            end_date = datetime.strptime(inputs['end'], "%Y-%m-%dT%H:%M:%S")
            dt = inputs['dt']
            struct_json_file = inputs['struct_json']
    except IOError:
        sys.exit('platform_nc: An error trying to read the file.')
    except KeyError:
        sys.exit('platform_nc: An error with a key')
    except ValueError:
        sys.exit('platform_nc: Non-numeric data found in the file.')
    except Exception as err:
        print(err)
        sys.exit(f'platform_nc: Error with the input {input_json_file}')

    nc_file_name = json_struct_to_nc(struct_json_file)

    variables = [{'ln': 83, 'name_nc': 'DRYT', 'depth': -4.5, 'DM': 'R'},
                 {'ln': 81, 'name_nc': 'WSPD', 'depth': -7, 'DM': 'R'},
                 {'ln': 82, 'name_nc': 'WDIR', 'depth': -7, 'DM': 'R'},
                 {'ln': 86, 'name_nc': 'RELH', 'depth': -4.5, 'DM': 'R'},
                 {'ln': 10003, 'name_nc': 'GSPD', 'depth': -7, 'DM': 'R'},
                 {'ln': 10015, 'name_nc': 'GDIR', 'depth': -7, 'DM': 'R'},
                 {'ln': 10018, 'name_nc': 'DEWT', 'depth': -4.5, 'DM': 'R'},
                 {'ln': 20003, 'name_nc': 'TEMP', 'depth': 0.5, 'DM': 'R'},
                 {'ln': 20019, 'name_nc': 'TEMP', 'depth': 3.5, 'DM': 'R'},
                 {'ln': 20005, 'name_nc': 'PSAL', 'depth': 0.5, 'DM': 'R'},
                 {'ln': 20018, 'name_nc': 'PSAL', 'depth': 3.5, 'DM': 'R'},
                 {'ln': 20004, 'name_nc': 'CNDC', 'depth': 0.5, 'DM': 'R'},
                 {'ln': 20022, 'name_nc': 'CNDC', 'depth': 3.5, 'DM': 'R'},
                 {'ln': 20000, 'name_nc': 'HCSP', 'depth': 8, 'DM': 'R'},  # TODO: change x, y velocity to direction, speed
                 {'ln': 20002, 'name_nc': 'HCDT', 'depth': 8, 'DM': 'R'}]

    aux_vars = [{'ln': 20804, 'name_nc': 'HC-X', 'depth': 8, 'DM': 'R'},  # TODO: change x, y velocity to direction, speed
                {'ln': 20803, 'name_nc': 'HC-Y', 'depth': 8, 'DM': 'R'}]

    ln_estacion = 15001  # Cortegada

    first = True
    for variable in variables:
        df = read_data(db_json_file, start_date, end_date, ln_estacion,  variable['ln'])
        name = variable['name_nc'] + '_' + str(variable['depth'])
        flag_name = name + '_QC'
        df.rename(columns={'InstanteLectura': 'TIME',
                           'Valor': name,
                           'LnCodigoValidacion': flag_name},
                  inplace=True)

        if first:
            result = df
            first = False
        else:
            result = pd.merge(result, df,  how='outer', on='TIME')


        if variable['name_nc'] == 'HCSP':
            name_module = name
            name_module_flag = flag_name
        if variable['name_nc'] == 'HCDT':
            name_dir = name
            name_dir_flag = flag_name

    # Cambio de u, v a mod,dir
    first = True
    for aux_var in aux_vars:
        df = read_data(db_json_file, start_date, end_date, ln_estacion, aux_var['ln'])
        name = aux_var['name_nc'] + '_' + str(aux_var['depth'])
        flag_name = name + '_QC'
        df.rename(columns={'InstanteLectura': 'TIME',
                           'Valor': name,
                           'LnCodigoValidacion': flag_name},
                  inplace=True)

        if first:
            result_aux = df
            first = False
        else:
            result_aux = pd.merge(result_aux, df, how='outer', on='TIME')
    module = []

    direction = []
    for index, row in result_aux.iterrows():
        u = row[1]
        v = row[3]
        m, d = uv2modtheta(u, v)
        module.append(m)
        direction.append(d)

    result[name_module] = module
    result[name_module_flag] = result_aux[2]
    result[name_dir] = direction
    result[name_dir_flag] = result_aux[2]




    with netCDF4.Dataset(nc_file_name, mode='a') as nc_file:
        # TIME
        dates = date_range(start_date, end_date, dt)
        times = (nc_file["TIME"])
        times[:] = netCDF4.date2num(dates, units=times.units, calendar=times.calendar)
        qc_times = (nc_file["TIME_QC"])
        length, = times.shape
        qc_times[:] = numpy.ones(length, dtype="i1")

        for variable in variables:
            nc_name = variable['name_nc']
            qc_nc_name = nc_name + '_QC'
            dm_nc_name = nc_name + '_DM'

            name = variable['name_nc'] + '_' + str(variable['depth'])
            flag_name = name + '_QC'

            # Find out the depth
            depths = (nc_file['DEPH'][0])
            v_depth = variable['depth']
            i_depth = numpy.where(depths == v_depth)[0][0]

            nc_var = (nc_file[nc_name])
            values = numpy.array(result[name])

            nc_var[:, i_depth] = values

            nc_var_qc = (nc_file[qc_nc_name])
            values_qc = numpy.array(result[flag_name])
            nc_var_qc[:, i_depth] = values_qc

            nc_var_dm = (nc_file[dm_nc_name])
            nc_var_dm[:, i_depth] = variable['DM']

            print(f'Variable {nc_name} at {v_depth} saved in file {nc_file_name}')




if __name__ == '__main__':
    input_json_file = 'platform_nc.json'
    platform_nc(input_json_file)
