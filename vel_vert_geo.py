#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 14:06:13 2022

@author: ladsin
"""
#importando bibliotecas
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo

#dataset

file_0 = xr.open_dataset(
    '/home/ladsin/Área de Trabalho/Coqueiro/dados/fis_nuvens/dados_fis_nuvens/24_02_19_continental_multi.nc',
    decode_times = False
    )

file_1 = xr.open_dataset(
    '/home/ladsin/Área de Trabalho/Coqueiro/dados/fis_nuvens/dados_fis_nuvens/24_02_19_continental_multi.nc'
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

file_2 = xr.open_dataset(
    '/home/ladsin/Área de Trabalho/Coqueiro/dados/fis_nuvens/dados_fis_nuvens/24_02_19_continental_single.nc'
    ).metpy.parse_cf()


file_2 = file_2.assign_coords(dict(
    longitude = (((file_2.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')


#
#extent
lon_slice = slice(-75., -20.)
lat_slice = slice(10., -60.)

#pega as lat/lon
lats = file_1.latitude.sel(latitude=lat_slice).values
lons = file_1.longitude.sel(longitude=lon_slice).values

#seta as variaveis
level = 500

for i in range(len(file_1.variables['time'])):
    
    geopotencial = file_1.z.metpy.sel(
        time = file_1.time[i] , vertical=level, latitude=lat_slice, longitude=lon_slice).metpy.unit_array.squeeze()/10
    
    100
    w = file_1.w.metpy.sel(
        time = file_1.time[i] , vertical=level, latitude=lat_slice, longitude=lon_slice).metpy.unit_array.squeeze()
    
    pnmm = file_2.msl.metpy.sel(
        time = file_1.time[i] , latitude=lat_slice, longitude=lon_slice).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa
    
    
    dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
    
    
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    
    shapefile = list(
        shpreader.Reader(
        '/home/ladsin/br_unidades_da_federacao/BR_UF_2019.shp'
        ).geometries()
        )
    
    ax.add_geometries(
        shapefile, ccrs.PlateCarree(), 
        edgecolor = 'black', 
        facecolor='none', 
        linewidth=1.0
        )
    
    # adiciona continente e bordas
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(),
                      color='gray',
                      alpha=1.0, 
                      linestyle='--', 
                      linewidth=0.5, 
                      xlocs=np.arange(-180, 180, 5), 
                      ylocs=np.arange(-90, 90, 5), 
                      draw_labels=True
                      )
    gl.top_labels = False
    gl.right_labels = False
    
    # intevalos da pnmm
    intervalo_min2 = np.amin(np.array(geopotencial))
    intervalo_max2 = np.amax(np.array(geopotencial))
    interval_2 = 50            # de quanto em quanto voce quer que varie
    levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)
    
    # intevalos da divergencia - umidade
    intervalo_min3 = np.amin(np.array(w))
    intervalo_max3 = np.amax(np.array(w))
    interval_3 = 1             # de quanto em quanto voce quer que varie
    levels_3 = np.arange(intervalo_min3, intervalo_max3, interval_3)
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # plota a imagem geopotencial
    sombreado = ax.contourf(lons, 
                            lats, 
                            w, 
                            cmap = 'Blues', 
                            levels = levels_3, 
                            extend = 'both'
                            )
    
    # plota a imagem pressao
    contorno = ax.contour(lons,
                          lats, 
                          geopotencial, 
                          colors='black', 
                          linewidths=1, 
                          levels=levels_2
                          )
    
    ax.clabel(contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
              fmt = '%5.0f', 
              colors= 'black'
              )
    
    # adiciona legenda 
    barra_de_cores = plt.colorbar(sombreado, 
                                  orientation = 'horizontal', 
                                  pad=0.04, 
                                  fraction=0.04
                                  )
    
    font_size = 20 # Adjust as appropriate.
    barra_de_cores.ax.tick_params(labelsize=font_size)
    
    # Getting the file time and date
    add_seconds = int(file_0.variables['time'][i])
    
    date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
    date_formatted = date.strftime('%Y-%m-%d %H')
    	
    # Add a title
    plt.title('Velocidade vertical e geopotencial - '+ str(level) +' hPa',
              fontweight='bold', 
              fontsize=35, 
              loc='left'
              )
    
    plt.title(f'{date_formatted}', 
              fontsize=20, 
              loc='right'
              )
    
    
    #--------------------------------------------------------------------------
    # Salva imagem
    plt.savefig(f'vel_vert_geo {date_formatted}.png', bbox_inches='tight')
    