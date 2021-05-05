# -*- coding: utf-8 -*-

import matplotlib
from matplotlib import pyplot as plt

import numpy as np
import pandas as pd

from mpl_toolkits.basemap import Basemap

import xarray as xr

import re

from collections import OrderedDict

from datetime import datetime, timedelta

from scipy.spatial import cKDTree, KDTree

from pyproj import Proj

import argparse

from glob import glob

import json

import os

import logging

logging.basicConfig(format='%(asctime)s - %(levelname)s - %(message)s', level=logging.INFO)

deg2rad = np.pi/180 

dtypes        = {"TIME"        : 'float64',
                 "DEPH"        : 'float32',
                 "BEAR"        : 'float32',
                 "RNGE"        : 'float32',
                 "LONGITUDE"   : 'float32',
                 "LATITUDE"    : 'float32',
                 "XDST"        : 'int32',
                 "YDST"        : 'int32',
                 "RDVA"        : 'int16',
                 "DRVA"        : 'int32',
                 "EWCT"        : 'int16',
                 "NSCT"        : 'int16',
                 "MAXV"        : 'int16',
                 "MINV"        : 'int16',
                 "ESPC"        : 'int16',
                 "ETMP"        : 'int16',
                 "ERSC"        : 'int16',
                 "ERTC"        : 'int16',
                 "SPRC"        : 'int16',
                 "NARX"        : 'int8',
                 "NATX"        : 'int8',
                 "SLTR"        : 'int32',
                 "SLNR"        : 'int32',
                 "SLTT"        : 'int16',
                 "SLNT"        : 'int16',
                 "TIME_QC"     : 'int8',
                 "POSITION_QC" : 'int8',
                 "DEPH_QC"     : 'int8',
                 "QCflag"      : 'int8',
                 "OWTR_QC"     : 'int8',
                 "MDFL_QC"     : 'int8',
                 "VART_QC"     : 'int8',
                 "CSPD_QC"     : 'int8',
                 "AVRB_QC"     : 'int8',
                 "RDCT_QC"     : 'int8'}

scale_factors = {"XDST"        : 0.001,
                 "YDST"        : 0.001,
                 "RDVA"        : 0.001,
                 "DRVA"        : 0.001,
                 "EWCT"        : 0.001,
                 "NSCT"        : 0.001,
                 "ESPC"        : 0.001,
                 "ETMP"        : 0.001,
                 "MAXV"        : 0.001,
                 "MINV"        : 0.001,
                 "ERSC"        : 1,
                 "ERTC"        : 1,
                 "XDST"        : 0.001,
                 "YDST"        : 0.001,
                 "SPRC"        : 1,
                 "NARX"        : 1,
                 "NATX"        : 1,
                 "SLTR"        : 0.001,
                 "SLNR"        : 0.001,
                 "SLTT"        : 0.001,
                 "SLNT"        : 0.001,
                 "TIME_QC"     : 1,
                 "POSITION_QC" : 1,
                 "DEPH_QC"     : 1,
                 "QCflag"      : 1,
                 "OWTR_QC"     : 1,
                 "MDFL_QC"     : 1,
                 "VART_QC"     : 1,
                 "CSPD_QC"     : 1,
                 "AVRB_QC"     : 1,
                 "RDCT_QC"     : 1}

add_offsets = {}

for key, value in scale_factors.items():
    if isinstance(value, float):

        scale_factors[key] = np.float32(scale_factors[key])
        add_offsets[key]   = np.float32(0)

    else:

        # Generamos un conversor de tipo a partir del tipo de la variable:
        conversor           = np.dtype(dtypes[key])
    
        # Utilizamos el conversor para recodificar un tipo nativo de python a un escalar tipo numpy:
        scale_factors[key] = np.int_(scale_factors[key]).astype(conversor)
        add_offsets[key]   = np.int_(0).astype(conversor)


_FillValues   = {}

for key, value in dtypes.items():
    if 'float' in value:
        _FillValues[key] = np.finfo(dtypes[key]).min+1
    else:
        _FillValues[key] = np.iinfo(dtypes[key]).min+1

class Radial():

    """
    Clase de abstracción para la lectura y procesamiento de ficheros radiales (.ruv)

    Atributos
    ---------

    Metodos
    -------
    """

    def __init__(self, fichero):

        """
        Constructor

        Parametros
        ----------
        fichero: Fichero .ruv con las velocidades radiales
        """

        # El archivo tiene que ser abierto como binary:
        contenido = [linea.decode('utf-8').replace('%','').replace('\n','') for linea in open(fichero, 'rb').readlines()
                     if '%%' not in  str(linea)]

        metadatos = [linea for linea in contenido if 'Table' not in linea]
        metadatos = dict([(linea.split(':')[0],linea.split(':')[1]) for linea in metadatos if ':' in str(linea) ])

        # Parseamos algunos metadatos que necesitaremos:
        self.Origin                   = np.array(metadatos['Origin'].split(),dtype=float)
        self.RangeEnd                 = int(metadatos['RangeEnd'])
        self.RangeResolutionKMeters   = float(metadatos['RangeResolutionKMeters'])
        self.AntennaBearing           = float(metadatos['AntennaBearing'].replace('True',''))
        self.AngularResolution        = float(metadatos['AngularResolution'].replace('Deg',''))
        self.TimeStamp                = datetime.strptime(metadatos['TimeStamp'],' %Y %m %d %H %M %S')

        # Líneas inicial y final de las tablas:
        starts  = np.arange(len(contenido))[['TableStart' in linea for linea in contenido]]
        ends    = np.arange(len(contenido))[['TableEnd' in linea for linea in contenido]]
        lengths = ends - starts - 1

        # Linea que contiene el header:
        columns = np.arange(len(contenido))[['TableColumnTypes' in linea for linea in contenido]]

        tablas = []

        # Aquí podemos aplicar los cambios en los nombres de las variables:

        headers    = [contenido[indice].split(':')[1].split() for indice in columns]
        headers[0] = ['LOND', 'LATD', 'EWCT', 'NSCT', 'OWTR_QC', 'ESPC', 'ETMP', 'MAXV', 'MINV', 'ERSC', 'ERTC', 'XDST', 'YDST', 'RNGE', 'BEAR', 'RDVA', 'DRVA', 'SPRC']
        ## Originales: LOND    LATD    VELU    VELV    VFLG       ESPC    ETMP    MAXV    MINV    ERSC    ERTC    XDST    YDST    RNGE    BEAR    VELO    HEAD    SPRC 

        for i in range(3):

            if lengths[i] != 0:

                start = starts[i] + 1
                end   = ends[i]

                tablas.append(pd.DataFrame(np.array([linea.split() for linea in contenido[start:end]],dtype=float), columns=headers[i]))

        # Eventualmente pueden aparecer datos erroneos en estas variables:
        tablas[0].ESPC[tablas[0].ESPC==999.00] =  np.nan
        tablas[0].ETMP[tablas[0].ETMP==999.00] =  np.nan

        # Aquí aplicamos los factores de conversión necesarios:
        tablas[0].EWCT /=  100.
        tablas[0].NSCT /=  100.
        tablas[0].RDVA /= -100.
        tablas[0].MINV /= -100.
        tablas[0].MAXV /= -100.
        tablas[0].ESPC /=  100.
        tablas[0].ETMP /=  100.
        tablas[0].ERSC /=  1.
        tablas[0].ERTC /=  1.
        tablas[0].SPRC /=  1.
       
        self.metadatos = metadatos
        self.tablas    = tablas

    def to_grid(self, grid):

        # Busqueda cKDTree:
        nearest = cKDTree(np.column_stack([grid.longitud.values.flatten(), grid.latitud.values.flatten()]))
        puntos  = np.column_stack([self.tablas[0].LOND.values, self.tablas[0].LATD.values])
        distancias, vecinos = nearest.query(puntos)

        variables = ['EWCT', 'NSCT', 'OWTR_QC', 'MINV', 'MAXV', 'RDVA', 'DRVA', 'ESPC', 'ETMP', 'ERSC', 'ERTC', 'SPRC']

        self.variables = OrderedDict()

        # Necesitamos completar la lista de coordenadas:
        delta = self.TimeStamp - datetime(1950,1,1)
        self.variables['TIME'] = xr.DataArray([delta.days + delta.seconds/86400], dims   = {'TIME' : 1})
        self.variables['DEPH'] = xr.DataArray([0.], dims   = {'DEPTH' : 1})

        for variable in variables:
            
            # Creamos la matriz que se llenará con los datos:
            tmp          = np.ones_like(grid.longitud.values.flatten())*np.nan

            # Asignamos los vecinos más próximos:
            tmp[vecinos] = self.tablas[0][variable]

            # Volvemos a la forma original:
            tmp          = tmp.reshape(grid.longitud.shape)

            # Creamos el DataArray:
            if variable in ['EWCT', 'NSCT', 'OWTR_QC', 'MINV', 'MAXV', 'RDVA', 'DRVA', 'ESPC', 'ETMP', 'ERSC', 'ERTC', 'SPRC']:

                # Crecemos en DEPTH:
                tmp = np.expand_dims(tmp,axis=0)

                # Crecemos en TIME:
                tmp = np.expand_dims(tmp,axis=0)

                self.variables[variable] = xr.DataArray(tmp, 
                                           dims   = {'TIME' : 1, 'DEPTH' : 1, 'BEAR' : grid.nBEAR, 'RNGE' : grid.nRNGE},
                                           coords = {'TIME' : self.variables['TIME'], 'DEPH' : self.variables['DEPH'], 'BEAR' : grid.BEAR , 
                                                     'RNGE' : grid.RNGE, 'LONGITUDE' : grid.longitud, 'LATITUDE' : grid.latitud,
                                                     'XDST' : grid.X, 'YDST' : grid.Y})

                # Encoding de las variables en el fichero:
                self.variables[variable].encoding["scale_factor"] = scale_factors[variable]
                self.variables[variable].encoding["add_offset"  ] = add_offsets[variable]
                self.variables[variable].encoding["dtype"       ] = dtypes[variable]
                self.variables[variable].encoding["_FillValue"  ] = _FillValues[variable]

    def QC_control(self):

        """
        Método para el control de calidad de los datos

        Parametros
        ----------
        Utiliza la lista de variables del objeto, que deben estar ya ongrid.
        """

        # Construimos alias de las variables para no tener que reescribir mucho código:

        bear   = self.variables['RDVA'].BEAR*deg2rad
        lond   = self.variables['RDVA'].LONGITUDE
        latd   = self.variables['RDVA'].LATITUDE
        X      = self.variables['RDVA'].XDST
        Y      = self.variables['RDVA'].YDST

        # owtr   = self.variables['']

        etmp   = self.variables['ETMP']
        head   = self.variables['DRVA']
        radVel = self.variables['RDVA']
        owtr   = self.variables['OWTR_QC']

        # Time quality flag:
        sdnTime_QCflag = 1
        self.variables['TIME_QC'] = xr.DataArray(sdnTime_QCflag, dims   = {'TIME' : 1}, coords = {'TIME' : self.variables['TIME']})
        
        # Position quality flag:
        sdnPosition_QCflag = radVel.copy()*0 + 1
        self.variables['POSITION_QC'] = sdnPosition_QCflag

        #sdnPosition_QCflag = netcdf.getConstant('NC_FILL_BYTE').*int8(ones(size(velu,1),size(velu,2),1));
        #sdnPosition_QCflag(velu~=netcdf.getConstant('NC_FILL_SHORT')) = 1;
        
        # Depth quality flag
        sdnDepth_QCflag = 1
        self.variables['DEPH_QC'] = xr.DataArray(sdnTime_QCflag, dims   = {'TIME' : 1}, coords = {'TIME' : self.variables['TIME']})

        # Variance Threshold QC test
        varVec = etmp**2 #varVec = etmp.^2;

        # Average Radial Bearing QC test:
        # avgBear_HEAD = head.values[~np.isnan(head).values].mean() # avgBear_HEAD = mean(head(~isnan(head)));
        avgBear      = bear.values[~np.isnan(bear).values].mean() # avgBear = mean(bear(~isnan(bear)));

        # Radial Count QC test:
        radVectors = int((~np.isnan(radVel)).sum()) # radVectors = sum(sum(~isnan(radVel)));

        # Over Water quality flags (cambiamos la codificación de la variable):
        self.variables['OWTR_QC'].values = np.select([owtr.values==0, owtr.values==128, np.isnan(owtr)], [1,4,np.nan])

        '''
        % Over Water quality flags
        overWater(owtr==0) = 1;
        overWater(owtr==128) = 4;
        '''
        
        # Velocity Threshold quality flags
        ## Velocity Threshold QC test
        maxspd_R = Radial_QC_params.VelThr # Dato basado en el netCDF # maxspd_R = Radial_QC_params.VelThr;
        velThr   = radVel.copy()
        velThr.values = np.select([np.abs(radVel) <= maxspd_R, np.abs(radVel) > maxspd_R, np.isnan(radVel)], [1,4,np.nan])

        self.variables['CSPD_QC'] = velThr

        '''
        velThr((abs(radVel) <= maxspd_R)) = 1
        velThr((abs(radVel) > maxspd_R)) = 4
        '''
        
        # Variance Threshold quality flags:
        varThr = varVec.copy()
        varThr.values = np.select([varVec <= Radial_QC_params.VarThr, varVec > Radial_QC_params.VarThr, np.isnan(varVec)], [1,4,np.nan])

        '''
        varThr((varVec > Radial_QC_params.VarThr)) = 4;
        varThr((varVec <= Radial_QC_params.VarThr)) = 1;
        '''

        # Set the QC flag for the current hour to 0 (no QC performed)
        tempDer  = radVel.copy()
        tempDer *= 0

        self.variables['VART_QC'] = tempDer

        # Median Filter QC test
        radVelMedianFiltered = radVel.copy()
        medFilt              = radVel.copy()

        condicion = ~np.isnan(radVel[0,0])

        nt, nd, nBear, nRange = radVel.shape

        tmp, bear = np.meshgrid( radVel.RNGE, radVel.BEAR)

        for i in range(nRange):
            for j in range(nBear):

                if condicion[j,i]:

                    refBear = bear[j,i]

                    # Condición de puntos que están a menos de una distancia dada:
                    ventana  = np.sqrt((X-X[j,i])**2 + (Y-Y[j,i])**2) <= Radial_QC_params.MedFilt[0]

                    # Condición de distancia angular (con codigo para controlar el círculo):
                    dif = bear - refBear

                    dif[dif>=np.pi] -= 2*np.pi
                    dif[dif<-np.pi] += 2*np.pi

                    ventana &= np.abs(dif) <= Radial_QC_params.MedFilt[1]

                    # Los datos no pueden ser nan:
                    ventana &= condicion

                    radVelMedianFiltered[0,0,j,i] = np.median(radVel[0,0].values[ventana])
                    #print('Número de puntos: %i' % ventana.sum())
                    #print('Velocidad filtrada: %f' % radVelMedianFiltered[0,0,j,i])


        medFilt[:] = np.select([np.abs(radVelMedianFiltered - radVel) <= Radial_QC_params.MedFilt[2], 
                                np.abs(radVelMedianFiltered - radVel) >  Radial_QC_params.MedFilt[2], np.isnan(radVelMedianFiltered)], [1,4,np.nan])

        self.variables['MDFL_QC'] = medFilt

        # Average Radial Bearing quality flag:
        if ((avgBear >= Radial_QC_params.AvgRadBear[0]) & (avgBear <= Radial_QC_params.AvgRadBear[1])):
            avgRadBear = 1
        else:
            avgRadBear = 4

        self.variables['AVRB_QC'] = xr.DataArray(avgRadBear, dims   = {'TIME' : 1}, coords = {'TIME' : self.variables['TIME']})
        
        # Radial Count quality flag:
        if (radVectors > Radial_QC_params.RadCnt):
            radCount = 1
        else:
            radCount = 4

        self.variables['RDCT_QC'] = xr.DataArray(avgRadBear, dims   = {'TIME' : 1}, coords = {'TIME' : self.variables['TIME']})

        # Populate the overall quality variable:
        condicion = (self.variables['CSPD_QC'] == 1) & \
                    (self.variables['OWTR_QC'] == 1) & \
                    (self.variables['MDFL_QC'] == 1) & \
                    (self.variables['AVRB_QC'] == 1) & \
                    (self.variables['RDCT_QC'] == 1)

        isNan = np.isnan(self.variables['CSPD_QC'])

        self.variables['QCflag'] = self.variables['CSPD_QC'].copy()
        self.variables['QCflag'].values = np.select([condicion & ~isNan, ~condicion & ~isNan, isNan], [1,4,np.nan])

        '''
        if(velThr(ii,jj) ~= netcdf.getConstant('NC_FILL_BYTE'))
            if((velThr(ii,jj) == 1) && (overWater(ii,jj) == 1) && (medFilt(ii,jj) == 1) && (avgRadBear == 1) && (radCount == 1))
                overall(ii,jj) = 1;
            else
                overall(ii,jj) = 4;
        '''

        # Terminamos ajustando algunos parámetros de las variables:        
        for variable in ['TIME_QC', 'POSITION_QC', 'DEPH_QC', 'QCflag', 'OWTR_QC', 'MDFL_QC', 'VART_QC', 'CSPD_QC', 'AVRB_QC', 'RDCT_QC']:

                # Encoding de las variables en el fichero:
                self.variables[variable].encoding["scale_factor"] = scale_factors[variable]
                self.variables[variable].encoding["add_offset"  ] = add_offsets[variable]
                self.variables[variable].encoding["dtype"       ] = dtypes[variable]
                self.variables[variable].encoding["_FillValue"  ] = _FillValues[variable]


    def to_netcdf(self, fichero):

        radar = re.findall("[A-Z]{4}", fichero.split('/')[-1])[0]
        fecha = datetime.strptime('%s%s%s%s' % tuple(re.findall("\d+", fichero.split('/')[-1])),'%Y%m%d%H%M' ) 

        logging.info('Fichero: %s Radar: %s' % (fichero, radar))

        # Info de la proyección:
        self.variables['crs'] = xr.DataArray(np.int16(0),)

        # Datos SDN:
        SDN_EDMO_CODEs = {'PRIO' : 4841, 'SILL' : 2751, 'VILA' : 4841}

        self.variables['SDN_EDMO_CODE'] = xr.DataArray(np.int16([[SDN_EDMO_CODEs[radar]]]), dims   = {'TIME' : 1, 'MAXINST' : 1})

        cadena = b'HFR-Galicia'
        n      = len(cadena)
        self.variables['SDN_CRUISE'] = xr.DataArray(np.array([cadena]), dims   = {'TIME' : 1})

        cadena = ('HFR-Galicia-%s' % radar).encode()
        n      = len(cadena)
        self.variables['SDN_STATION'] = xr.DataArray(np.array([cadena]), dims   = {'TIME' : 1})

        cadena = ('HFR-Galicia-%s_%sZ' % (radar, self.TimeStamp.isoformat())).encode()
        n      = len(cadena)
        self.variables['SDN_LOCAL_CDI_ID'] = xr.DataArray(np.array([cadena]), dims   = {'TIME' : 1})

        cadena = b'http://opendap.intecmar.gal/thredds/catalog/data/nc/RADAR_HF/Galicia/catalog.html'
        n      = len(cadena)
        self.variables['SDN_REFERENCES'] = xr.DataArray(np.array([cadena]), dims   = {'TIME' : 1})

        cadena = b"<sdn_reference xlink:href=\"http://opendap.intecmar.gal/thredds/catalog/data/nc/RADAR_HF/Galicia/catalog.html\" xlink:role=\"\" xlink:type=\"URL\"/>"
        n      = len(cadena)
        self.variables['SDN_XLINK'] = xr.DataArray(np.array([[cadena]]), dims   = {'TIME' : 1, 'REFMAX' : 1})

        # Otras:
        siteLat, siteLon = self.Origin
        self.variables['SLTR'] = xr.DataArray([[siteLat]], dims   = {'TIME' : 1, 'MAXSITE' : 1})
        self.variables['SLNR'] = xr.DataArray([[siteLon]], dims   = {'TIME' : 1, 'MAXSITE' : 1})
        self.variables['SLTT'] = xr.DataArray([[siteLat]], dims   = {'TIME' : 1, 'MAXSITE' : 1})
        self.variables['SLNT'] = xr.DataArray([[siteLon]], dims   = {'TIME' : 1, 'MAXSITE' : 1})

        cadena = ('%s' % radar).encode()
        n      = len(cadena)
        self.variables['SCDR'] = xr.DataArray(np.array([[cadena]]), dims   = {'TIME' : 1, 'MAXSITE' : 1})
        self.variables['SCDT'] = xr.DataArray(np.array([[cadena]]), dims   = {'TIME' : 1, 'MAXSITE' : 1})

        numSites = 1
        self.variables['NARX'] = xr.DataArray([numSites], dims   = {'TIME' : 1})
        self.variables['NATX'] = xr.DataArray([numSites], dims   = {'TIME' : 1})

        for variable in ['SLTT','SLNT','SLTR','SLNR','NARX','NATX']:
                # Encoding de las variables en el fichero:
                self.variables[variable].encoding["scale_factor"] = scale_factors[variable]
                self.variables[variable].encoding["add_offset"  ] = add_offsets[variable]
                self.variables[variable].encoding["dtype"       ] = dtypes[variable]
                self.variables[variable].encoding["_FillValue"  ] = _FillValues[variable]

        # Generamos el xarra.Dataset. radial.variables contienen los xr.DataArray necesarios:
        dataset = xr.Dataset(self.variables)

        # Atributos globales:
        ## Leemos los atributos específicos de cada radar:
        f = open('%s.json' % radar)
        atributos = json.loads(f.read())
        f.close()

        ## Atributos del fichero radial que serán sobreescritos con los datos del fichero radial:
        atributos_fichero  = ['AngularResolution', 'AntennaBearing', 'BraggHasSecondOrder', 'BraggSmoothingPoints', 
                              'DopplerResolutionHzPerBin', 'FirstOrderCalc', 'FirstOrderMethod', 'MergeMethod', 'MergedCount', 
                              'PatternAmplitudeCalculations', 'PatternAmplitudeCorrections', 'PatternMethod', 'PatternPhaseCalculations', 
                              'PatternPhaseCorrections', 'PatternResolution', 'RadialBraggNoiseThreshold', 'RadialBraggPeakDropOff', 'RadialBraggPeakNull', 
                              'RadialMinimumMergePoints', 'RadialMusicParameters', 'RangeEnd', 'RangeResolutionKMeters', 'RangeStart', 
                              'ReferenceBearing', 'SpatialResolution', 'SpectraDopplerCells', 'SpectraRangeCells', 'TransmitBandwidthKHz', 
                              'TransmitCenterFreqMHz', 'TransmitSweepRateHz', 'UUID']

        ## Creamos algunos atributos:
        atributos['id'] = 'HFR-Galicia-%s_%sZ' % (radar, self.TimeStamp.isoformat())

        atributos['time_coverage_start'] = '%sZ' % (self.TimeStamp-timedelta(minutes=30)).isoformat()
        atributos['time_coverage_end']   = '%sZ' % (self.TimeStamp+timedelta(minutes=30)).isoformat()

        ahora = datetime(*datetime.now().timetuple()[0:6]).isoformat()
        atributos['date_created']        = '%sZ' % ahora
        atributos['metadata_date_stamp'] = '%sZ' % ahora
        atributos['date_modified']       = '%sZ' % ahora
        atributos['date_issued']         = '%sZ' % ahora

        atributos['history']             = '%s data collected. %s netCDF file created and sent to European HFR Node' % (self.TimeStamp.isoformat(), ahora)

        for atributo_fichero in atributos_fichero:
            try:
                atributos[atributo_fichero] = self.metadatos[atributo_fichero]
            except:
                logging.info('No puedo cargar el atributo --> %s del fichero radial' % atributo_fichero)

        ## ... y los insertamos
        dataset.attrs = atributos

        # Atributos de las variables:
        f = open('variables.json')
        atributos = json.loads(f.read())
        f.close()

        # Los tipos de los atributos valid_min/max son deserializados incorrectamente:
        for var in dataset:
            for key, value in atributos[var].items():
                if isinstance(atributos[var][key], int):

                    # Generamos un conversor de tipo a partir del tipo de la variable:
                    conversor           = np.dtype(dtypes[var])
                
                    # Utilizamos el conversor para recodificar un tipo nativo de python a un escalar tipo numpy:
                    atributos[var][key] = np.int_(atributos[var][key]).astype(conversor)

                elif isinstance(atributos[var][key],list):
    
                    # Generamos un conversor de tipo a partir del tipo de la variable:
                    conversor           = np.dtype(dtypes[var])
                
                    # Utilizamos el conversor para recodificar un tipo nativo de python a un escalar tipo numpy:
                    atributos[var][key] = np.array(atributos[var][key]).astype(conversor)

        for var in dataset:
            dataset[var].attrs = atributos[var]
        
        # Completamos coordenadas y dimensiones que xArray procesa de forma automática una vez creado el xr.Dataset a partir del diccionario de variables:
        for var in ['TIME', 'DEPH', 'BEAR', 'RNGE']:

            dataset[var].encoding["dtype"       ] = dtypes[var]
            dataset[var].encoding["_FillValue"  ] = None

            dataset[var].attrs = atributos[var]

        for var in ['LONGITUDE','LATITUDE']:

            dataset[var].encoding["dtype"       ] = dtypes[var]
            dataset[var].encoding["_FillValue"  ] = _FillValues[var]

            dataset[var].attrs = atributos[var]

        for var in ['XDST','YDST']:
            dataset[var].encoding["scale_factor"] = scale_factors[var]
            dataset[var].encoding["add_offset"  ] = add_offsets[var]
            dataset[var].encoding["dtype"       ] = dtypes[var]
            dataset[var].encoding["_FillValue"  ] = _FillValues[var]

            dataset[var].attrs = atributos[var]

        for var in ['DEPH', 'BEAR', 'RNGE','LONGITUDE','LATITUDE','XDST','YDST']:

            # Generamos un conversor de tipo a partir del tipo de la variable:
            conversor           = np.dtype(dtypes[var])
        
            # Utilizamos el conversor para recodificar un tipo nativo de python a un escalar tipo numpy:
            dataset[var].attrs['valid_min'] = np.float_(atributos[var]['valid_min']).astype(conversor)
            dataset[var].attrs['valid_max'] = np.float_(atributos[var]['valid_max']).astype(conversor)

        # Escribimos el netCDF:
        dataset.reset_coords(drop=False).to_netcdf('HFR-Galicia-%s_%s.nc' % (radar, fecha.strftime('%Y_%m_%d_%H%M'))) 

    def __repr__(self):
    
        return '<Radial class>'

class Radial_QC_params():

    """
    Clase estática para contener los umbrales.
    """

    VelThr      = 1.2 # (m/s)
    VarThr      = 1.  # (m2/s2?)
    tempDer_Thr = 0
    AvgRadBear  = [0., 70.]
    RadCnt      = 100
    MedFilt     = [5000,30*np.pi/180,1]  # 5km, 30 grados y 1m/s


class Grid():

    """
    Clase para la generación de la malla para albergar los datos

    Atributos
    ---------
    longitud, latitud: Matriz con las longitudes y latitudes reconstruidas

    Metodos
    -------
    """

    def __init__(self, radial, nBEAR=72, nRNGE=40):

        """
        Constructor

        Parametros
        ----------
        radial: Objeto de la clase Radial

        Parametros por defecto
        ----------------------
        nBEAR:  Número de direcciones
        nRNGE:  Número de distancias
        """

        self.nBEAR, self.nRNGE = nBEAR, nRNGE

        origen_lat, origen_lon = radial.Origin

        # Escogemos una proyección. Tmercator está bien. La idea es trabajar en un plano:
        #m = Basemap(llcrnrlon=-11.0, llcrnrlat=41.8, urcrnrlon=-8, urcrnrlat=44.5, resolution='h', projection='tmerc', lon_0=-8, lat_0=45)
        m = Basemap(llcrnrlon=-11.0, llcrnrlat=41.8, urcrnrlon=-8, urcrnrlat=44.5, resolution='h', projection='tmerc', lon_0=origen_lon, lat_0=origen_lat)

        # Necesito las coordenadas del origen y su proyección:
        origen_x, origen_y     = m(origen_lon, origen_lat)

        # Coordenadas polares de los puntos:
        RangeResolutionKMeters = radial.RangeResolutionKMeters

        AntennaBearing    = radial.AntennaBearing
        AngularResolution = radial.AngularResolution

        # Radios:
        RNGE = np.arange(nRNGE)*RangeResolutionKMeters*1000

        # Ángulos:
        BEAR = np.arange(nBEAR)*AngularResolution + AntennaBearing 
        # BEAR = np.sort(BEAR%360)*deg2rad
        BEAR = np.sort(BEAR%360)

        # Generamos la lista de vectores unitarios en las direcciones:
        X, Y = m.rotate_vector(np.sin(BEAR*deg2rad), np.cos(BEAR*deg2rad), np.repeat(origen_lon, len(BEAR)), np.repeat(origen_lat,len(BEAR)))

        X = np.array([RNGE*x + origen_x for x in X])
        Y = np.array([RNGE*y + origen_y for y in Y])

        # ... y las coordenadas esféricas reconstruidas:
        longitud, latitud = m(X,Y,inverse=True)

        # Preparamos las variables para guardar (las queremos en km no en m y referidas al origen de coordenadas):
        X    -= origen_x
        Y    -= origen_y
        X    /= 1000
        Y    /= 1000
        RNGE /= 1000

        # Guardamos las coordenadas proyectadas para trabajar en el plano:
        self.X = xr.DataArray(X, dims={'BEAR' : nBEAR, 'RNGE' : nRNGE}, coords={'BEAR' : BEAR, 'RNGE' : RNGE})
        self.Y = xr.DataArray(Y, dims={'BEAR' : nBEAR, 'RNGE' : nRNGE}, coords={'BEAR' : BEAR, 'RNGE' : RNGE})

        # ... que se guardan como xr.DataArray para su uso futuro en la definición de las variables:
        self.longitud = xr.DataArray(longitud, dims={'BEAR' : nBEAR, 'RNGE' : nRNGE}, coords={'BEAR' : BEAR, 'RNGE' : RNGE})
        self.latitud  = xr.DataArray(latitud , dims={'BEAR' : nBEAR, 'RNGE' : nRNGE}, coords={'BEAR' : BEAR, 'RNGE' : RNGE})

        # ... y las coordenadas polares de los puntos:
        self.RNGE,self.BEAR          = RNGE,BEAR

    def __repr__(self):
    
        return '<Grid class -> nBEAR: %i, nRNGE: %i>' % (self.nBEAR, self.nRNGE)


def VART_QC(ficheros):

    datasets = [xr.open_dataset(fichero) for fichero in ficheros] 
    radiales = [dataset.RDVA[0,0].values for dataset in datasets]

    radVel2h, radVel1h, radVel = radiales

    tempDer1h = np.full_like(radVel1h,4)

    tempDer_Thr = 1

    condicion  = np.abs(radVel   - radVel1h) < tempDer_Thr
    condicion &= np.abs(radVel2h - radVel1h) < tempDer_Thr

    tempDer1h[condicion] = 1

    condicion = np.isnan(radVel) | np.isnan(radVel2h)

    tempDer1h[condicion] = 0

    tempDer1h[np.isnan(radVel1h)] = np.nan

    datasets[1].VART_QC.values[0,0,:] = tempDer1h[:]

    # Redefinimos overall quality variable para incluir los cambios recientes en VART_QC:
    condicion = (datasets[1].variables['VART_QC'] == 1) & \
                (datasets[1].variables['CSPD_QC'] == 1) & \
                (datasets[1].variables['OWTR_QC'] == 1) & \
                (datasets[1].variables['MDFL_QC'] == 1) & \
                (datasets[1].variables['AVRB_QC'] == 1) & \
                (datasets[1].variables['RDCT_QC'] == 1)

    isNan = np.isnan(datasets[1].variables['CSPD_QC'])

    datasets[1].variables['QCflag'].values = np.select([condicion & ~isNan, ~condicion & ~isNan, isNan], [1,4,np.nan])

    datasets[1].to_netcdf('%s_new.nc' % ficheros[1].split('.')[0])


if __name__ == '__main__':

    # Parseamos la linea de comandos:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', help='File to transform', type=str)
    args = parser.parse_args()
    fichero = args.fichero
    radar = re.findall("[A-Z]{4}", fichero.split('/')[-1])[0]
    fecha = datetime.strptime('%s%s%s%s' % tuple(re.findall("\d+", fichero.split('/')[-1])), '%Y%m%d%H%M' )

    # Creamos el objeto radial para leer el fichero:
    radial = Radial(fichero)

    # Creamos la malla donde queremos inscribir la tabla:
    grd = Grid(radial)

    # Metemos la tabla en la malla:
    radial.to_grid(grd)

    # Generamos las variables de control de calidad:
    radial.QC_control()

    # Generamos el fichero NetCDF:
    radial.to_netcdf(fichero)

    ficheros = ['HFR-Galicia-%s_%s.nc' % (radar, (fecha + timedelta(hours=-i)).strftime('%Y_%m_%d_%H%M')) for i in range(3)]
    print(ficheros)
    condiciones = [os.path.isfile(fichero) for fichero in ficheros]
    print(condiciones)

    if np.all(condiciones):
        logging.info('Procesando VART_QC en %s' % ficheros[1])
        VART_QC(ficheros)
    else:
        logging.info('No VART_QC')


'''
    plt.pcolormesh(grd.longitud,grd.latitud,radial.variables['RNGE'])
    plt.grid()
    plt.colorbar()
    plt.show()

    plt.pcolormesh(grd.longitud,grd.latitud,radial.variables['BEAR'])
    plt.grid()
    plt.colorbar()
    plt.show()

    plt.plot(grd.longitud,grd.latitud,'k.')
    plt.plot(radial.tablas[0].LOND, radial.tablas[0].LATD,'r.')
    plt.grid()
    plt.show()

    ax = plt.subplot(111, projection='polar')
    ax.pcolormesh(grd.theta, grd.r,radial.variables['VELO'])
    ax.grid()
    plt.show()

    test = './datos/netcdf/HFR-Galicia-PRIO_2020_08_01_1200.nc'
    datos_test = xr.open_dataset(test)


    # Representación con CartoPy:
    land_50m = cfeature.NaturalEarthFeature('physical', 'land', '50m',
                                            edgecolor='k',
                                            facecolor=cfeature.COLORS['land'])

    land = cfeature.GSHHSFeature(scale='auto')


    proyeccion = ccrs.LambertConformal(central_longitude=-15.0, central_latitude=40)

    fig, ax = plt.subplots(1, 1, figsize=(9,6), subplot_kw=dict(projection=proyeccion))

    ax.set_extent([-11, -5.8, 41.7, 45.4], crs=ccrs.PlateCarree())
    #ax.stock_img()

    ax.add_feature(land)
    ax.add_feature(cfeature.BORDERS, edgecolor='gray')

    # Solo soportado para PlateCarree y Mercator:
    # gl = ax.gridlines(draw_labels=True, linewidth=1, color='black', alpha=0.5, linestyle='--')
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.5, linestyle='--')

    gl.xlines = True
    gl.ylines = True
    gl.xlocator = mticker.FixedLocator(np.arange(-11,-4,1))
    gl.ylocator = mticker.FixedLocator(np.arange(41,46,1))
       
    cset = ax.pcolormesh(grd.longitud, grd.latitud, radial.variables['VELO'][0,0,:], transform=ccrs.PlateCarree())

    #plt.show()


    fig, ax = plt.subplots(1, 1, figsize=(9,6), subplot_kw=dict(projection=proyeccion))

    ax.set_extent([-11, -5.8, 41.7, 45.4], crs=ccrs.PlateCarree())
    #ax.stock_img()

    ax.add_feature(land)
    ax.add_feature(cfeature.BORDERS, edgecolor='gray')

    # Solo soportado para PlateCarree y Mercator:
    # gl = ax.gridlines(draw_labels=True, linewidth=1, color='black', alpha=0.5, linestyle='--')
    gl = ax.gridlines(linewidth=1, color='black', alpha=0.5, linestyle='--')

    gl.xlines = True
    gl.ylines = True
    gl.xlocator = mticker.FixedLocator(np.arange(-11,-4,1))
    gl.ylocator = mticker.FixedLocator(np.arange(41,46,1))
       
    cset = ax.pcolormesh(datos_test.LONGITUDE, datos_test.LATITUDE, -100*datos_test.RDVA[0,0,:], transform=ccrs.PlateCarree())

    plt.show()
'''
