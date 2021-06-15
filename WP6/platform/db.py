#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**db.py**

* *Purpose:* Queries to MS SQL

* *python version:* 3.8
* *author:* Pedro Montero
* *license:* INTECMAR
* *requires:*

* *date:* 2020/02/13

* *version:* 0.0.1
* *version date:*

"""
import pyodbc
import pandas as pd
import pandas.io.sql as psql


def consulta(con_pm, consulta, parametros_consulta):
    """
    Fai unha consulta xeral a BD

    :param con_pm: list with conexions parameters
    :param consulta: string da consulta SQL
    :param parametros_consulta: parametros da consulta
    :return:
    """
    try:
        cnxn = pyodbc.connect("DRIVER=" + con_pm['driver'] +
                              ";SERVER=" + con_pm['server'] +
                              ";DATABASE=" + con_pm['database'])
        cursor = cnxn.cursor()
        cursor.execute(consulta, parametros_consulta)
        resultado = []
        for tupla in cursor.fetchall():
            resultado.append(tupla)

        cursor.close()
        cnxn.close()
        return resultado

    except pyodbc.Error as e:
        print(e)


def consulta_pd(con_pm, consulta, parametros_consulta):
    """
    Fai unha consulta xeral a BD

    :param con_pm: list with conexions parameters
    :param consulta: string da consulta SQL
    :param parametros_consulta: parametros da consulta
    :return:
    """
    try:
        cnxn = pyodbc.connect("DRIVER=" + con_pm['driver'] +
                              ";SERVER=" + con_pm['server'] +
                              ";DATABASE=" + con_pm['database'])
        sql_query = pd.read_sql_query(consulta, cnxn, params=parametros_consulta)
        cnxn.close()
        return sql_query

    except pyodbc.Error as e:
        print(e)
