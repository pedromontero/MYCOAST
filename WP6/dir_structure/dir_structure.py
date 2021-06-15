#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
**dir_structure.py**

* *Purpose:* Create a structure of folders as depicted in CMEMS InSituTAC PUM:
             PUM: Product User Manual for multiparameter Copernicus In Situ Tac (PUM)
             Copernicus Marine In Situ Tac Data Management Team (2020).
             Product User Manual for multiparameter Copernicus In Situ TAC (PUM).
             https://doi.org/10.13155/43494




* *python version:* 3.8
* *author:* Pedro Montero
* *license:* INTECMAR - MYCOAST
* *requires:*

* *date:*

* *version:* 1.0.0
* *version date:*

"""
import os
from datetime import datetime, timedelta


def latest(root):
    """
    Creates or updates the 'latest' directory, removing the obsolete first day folder

    :param root: str, initial path
    :return:
    """
    path_latest = os.path.join(root, 'latest')
    today = datetime.today()
    one_month = today + timedelta(days=-32)
    day = one_month

    # Create one last month folders if they don´t exist
    while day <= today:
        folder_day = day.strftime("%Y%m%d")
        path_folder_day = os.path.join(path_latest, folder_day)
        print(path_folder_day)
        os.makedirs(path_folder_day, exist_ok=True)
        day += timedelta(days=1)

    # remove first day folder
    first_day = one_month + timedelta(days=-1)
    folder_first_day = first_day.strftime("%Y%m%d")
    path_folder_first_day = os.path.join(path_latest, folder_first_day)
    # Check if exist
    if os.path.isdir(path_folder_first_day):
        # Check if empty
        if not os.listdir(path_folder_first_day):
            print(f'{path_folder_first_day} empty.\n It will be removed')
            try:
                os.rmdir(path_folder_first_day)
            except OSError as e:   # if failed, report it back to the user
                print("Error: %s - %s." % (e.filename, e.strerror))
        else:
            print(f'{path_folder_first_day} not empty.\n It will not be removed')
    else:
        print(f'{path_folder_first_day} not exist.\n It will not be removed')


def monthly(root, types):
    """
    Creates or updates the monthly directory
    :param root: str, initial path
    :param types: list of platforms
    :return:
    """

    today = datetime.today()
    current_year = today.year
    current_month = today.month
    initial_year = current_year - 5
    ym_end = current_year*100 + current_month

    path_monthly = os.path.join(root, 'monthly')
    for pla in types:
        for y in range(initial_year, current_year+1):
            for m in range(1, 13):
                ym_folder = y * 100 + m
                if ym_folder < ym_end:
                    str_ym_folder = str(ym_folder)
                    full_path = os.path.join(path_monthly, pla, str_ym_folder)
                    os.makedirs(full_path)


def dir_structure(root):
    """
    Create a structure of folders similar to PUM description

    :param root: str, initial path
    :return:
    """
    types = ['BO', 'CT', 'DB', 'DC', 'FB', 'GL', 'HF', 'ML', 'MO', 'PF', 'RF', 'SD', 'SM', 'TG', 'TS', 'XB', 'XX']
    if not os.path.exists(root):
        print(f'The root path: {root} not exist')
        print('Create it and rerun the script')
        quit()
    else:
        print(f'The root path will be {root}')

    # latest
    latest(root)

    # monthly
    monthly(root, types)


def main():
    root = os.path.join('', '../../data/out/inicio')
    dir_structure(root)


if __name__ == '__main__':
    main()
