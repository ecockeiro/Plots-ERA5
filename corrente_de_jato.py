#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  4 12:40:52 2022

@author: coqueiro
"""

#importando bibliotecas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cmocean

#dataset

file_1 = xr.open_dataset(
    '/home/coqueiro/Downloads/GFS_Global_0p25deg_20220903_1200.grib2.nc4'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

#extent
lon_slice = slice(-120., 10.)
lat_slice = slice(10., -70.)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

#seta as variaveis
level = 250 * units('hPa')

for i in range(len(file_1.variables['time1'])):
    
    
    u = file_1['u-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time1[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    v = file_1['v-component_of_wind_isobaric'].metpy.sel(
        time = file_1.time1[i], 
        vertical=level, 
        latitude=lat_slice, 
        longitude=lon_slice
        ).metpy.unit_array.squeeze()
    
    mag = np.sqrt(u**2+v**2)
    
    #data
    vtime = file_1.time1.data[i].astype('datetime64[ms]').astype('O')
    
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    
    shapefile = list(
        shpreader.Reader(
        '/home/coqueiro/ufrj/Estagio_supervisionado/shapefile/unidades_federativas/Brasil/UFEBRASIL.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, 
        ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=0.5
        )
    
    # adiciona continente e bordas
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5,
                      xlocs=np.arange(-180, 180, 10), 
                      ylocs=np.arange(-90, 90, 10), 
                      draw_labels=True
                      )
    gl.top_labels = False
    gl.right_labels = False
    gl.xlabel_style = {'size': 29, 'color': 'black'}
    gl.ylabel_style = {'size': 29, 'color': 'black'}
    
    # intevalos da corrente de jato
    intervalo_min3 = 30
    intervalo_max3 = 100
    interval_3 = 10            # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # corrente de jato
    img = plt.contourf(lons,
                       lats, 
                       mag, 
                       cmap=cmocean.cm.amp, 
                       levels = levels_3, 
                       extend='both')
    
    img2 = ax.contour(lons, 
                      lats, 
                      mag, 
                      colors='white', 
                      linewidths=0.3, 
                      levels=levels_3, 
                      transform=ccrs.PlateCarree())
    
    ax.streamplot(lons, 
                  lats, 
                  u, 
                  v, 
                  density=[6,6], 
                  linewidth=1.5, 
                  color='black', 
                  transform=ccrs.PlateCarree())
    
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(img, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # Add a title
    plt.title('Vorticidade relativa (1/s) - 1000 hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    #previsao
    #plt.title('Valid Time: {}'.format(vtime), fontsize=35, loc='right')
    #analise
    plt.title('Análise: {}'.format(vtime), fontsize=35, loc='right')
    
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'/home/coqueiro/ufrj/Estagio_supervisionado/imagens/corrente_de_jato/corrente_de_jato_{format(vtime)}.png', bbox_inches='tight')

