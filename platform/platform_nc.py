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
from datetime import datetime
from db import consulta_pd
import pandas as pd


def read_data(db_json_file, inicio, fin, ln_estacion,  ln_parametro):
    """

    :param db_json_file:
    :return:
    """

    try:
        with open(db_json_file, 'r') as f:
            con_parameters  = json.load(f, object_pairs_hook=OrderedDict)
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
    except IOError:
        sys.exit('An error occured trying to read the file.')
    except KeyError:
        sys.exit('An error with a key')
    except ValueError:
        sys.exit('Non-numeric data found in the file.')
    except Exception as err:
        print(err)
        sys.exit(f'Error with the input {input_json_file}')

    variables = [{'name': 'TEMPERATURE', 'flag': 'QC_TEMPERATURE', 'ln': 20003},
                 {'name': 'SALINITY', 'flag': 'QC_SALINITY', 'ln': 20005}]

    inicio = datetime(2020, 5, 15, 0, 0, 0)
    fin = datetime(2020, 5, 16, 0, 0, 0)

    ln_estacion = 15001  # Cortegada

    first = True
    for variable in variables:

        df = read_data(db_json_file, inicio, fin, ln_estacion,  variable['ln'])

        df.rename(columns={'InstanteLectura': 'TIME',
                           'Valor': variable['name'],
                           'LnCodigoValidacion': variable['flag']},
                  inplace=True)

        if first:
            result = df
            first = False
        else:
            result = pd.merge(result, df,  how='outer', on='TIME')
    print(result)



if __name__ == '__main__':
    input_json_file = 'platform_nc.json'
    platform_nc(input_json_file)
