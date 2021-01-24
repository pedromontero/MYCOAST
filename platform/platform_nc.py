#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**platform.py**

* *Purpose:* Read data from a platform and a json structure file and create a mycoast netcdf

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
from json_struct_to_nc.json_struct_to_nc import json_struct_to_nc


def date_range(start_date, end_date, dt):
    """
    Return a list of dates >= start_date and < end_date and dt =dt in seconds
    :param start_date:
    :param end_date:
    :param dt:
    :return:
    """
    dates = []
    delta = timedelta(seconds=dt)
    while start_date < end_date:
        dates.append(start_date)
        start_date += delta
    return dates


def read_data(db_json_file, inicio, fin, ln_estacion,  ln_parametro):
    """

    :param db_json_file:
    :return:
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



    str_inicio = inicio.strftime("%Y-%m-%d %H:%M:%S")
    str_fin = fin.strftime("%Y-%m-%d %H:%M:%S")

    cons = "SELECT TOP (1000) [InstanteLectura] ,[Valor] ,[LnCodigoValidacion] " \
           "FROM [tHstUltimosValidados] INNER JOIN [tCfgCnxEstaciones] " \
           "ON  [tHstUltimosValidados].lnEstacion = [tCfgCnxEstaciones].IdEstacion " \
           "WHERE lnEstacion = ? " \
           "AND InstanteLectura >= CAST(? AS date) " \
           "AND InstanteLectura < CAST(? AS date) " \
           "AND lnParametro = ?"
    params = (ln_estacion, str_inicio, str_fin, ln_parametro,)
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

    variables = [{'name': 'TEMPERATURE_SUP', 'flag': 'QC_TEMPERATURE', 'ln': 20003, 'depth': 0.5},
                 {'name': 'SALINITY', 'flag': 'QC_SALINITY', 'ln': 20005, 'depth': 0.5}]



    ln_estacion = 15001  # Cortegada

    # TIME
    dates = date_range(start_date, end_date, dt)
    nc_file = netCDF4.Dataset(nc_file_name, mode='a')
    times = (nc_file["TIME"])
    times[:] = netCDF4.date2num(dates, units=times.units, calendar=times.calendar)
    qc_times = (nc_file["TIME_QC"])
    length, = times.shape
    qc_times[:] = numpy.ones(length, dtype="i1")
    nc_file.close()



    first = True
    for variable in variables:

        df = read_data(db_json_file, start_date, end_date, ln_estacion,  variable['ln'])

        df.rename(columns={'InstanteLectura': 'TIME',
                           'Valor': variable['name'],
                           'LnCodigoValidacion': variable['flag']},
                  inplace=True)

        if first:
            result = df
            first = False
        else:
            result = pd.merge(result, df,  how='outer', on='TIME')
    nc_file = netCDF4.Dataset(nc_file_name, mode='a')

    # Find out the depth

    depths = (nc_file['DEPH'][0])
    print(depths)
    v_depth = variable['depth']
    i_depth = numpy.where(depths == v_depth)[0][0]
    print(f'La profundidad de {v_depth} es {i_depth}')


    temps = (nc_file["TEMP"])
    temp = numpy.array(round(result['TEMPERATURE_SUP']/0.01), dtype='i2')
    temp =result['TEMPERATURE_SUP']/0.001
    temp = round(temp).astype('int')
    temp = numpy.array(temp, dtype='i2')
    print(temp)
    temps[:, i_depth] = temp
    nc_file.close()




if __name__ == '__main__':
    input_json_file = 'platform_nc.json'
    platform_nc(input_json_file)
