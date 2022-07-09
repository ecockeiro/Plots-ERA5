#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  3 15:34:56 2022

@author: coqueiro
"""

import xarray as xr # ler netCDF e muito mais <3
import matplotlib.pyplot as plt  # biblioteca de plots
import matplotlib.colors                  # Matplotlib colors  
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cartopy, cartopy.crs as ccrs  # plot mapas
import cartopy.feature as cfeature  # operacao de desenho em imagens
import numpy as np  # importa o pacote numpy
import cmocean 
from cartopy.util import add_cyclic_point
import mygrads as mg #pip install mygrads
import metpy.calc as mpcalc
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader # Import shapefiles
from metpy.units import units
from metpy.interpolate import cross_section


file_1 = xr.open_dataset('/home/coqueiro/ufrj/micro/multi_bomba.nc', decode_times = False)
file_2 = xr.open_dataset('/home/coqueiro/ufrj/micro/pnmm_bomba.nc', decode_times = False)

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))).sortby('longitude')
file_2 = file_2.assign_coords(dict(
    longitude = (((file_2.longitude.values + 180) % 360) - 180))).sortby('longitude')

# limites de lat e lon
lat_min = -70.00
lat_max = 10.00
lon_min = -105.0
lon_max = -20.00

# seleciona uma extensao minima para o recorte de dados
extent = [lon_min, lon_max, lat_min, lat_max]

#geopotencial
lats= np.array(file_1.variables['latitude'][400:641])
lons= np.array(file_1.variables['longitude'][400:641])

#pressao
lats = file_2.variables['latitude'][400:641]
lons = file_2.variables['longitude'][400:641]

# define a extensao da imagem
img_extent = [lon_min, lon_max, lat_min, lat_max]
#------------------------------------------------------------------------------
# escolha o tamanho do plot em polegadas (largura x altura)
plt.figure(figsize=(20,20))

# usando a projeção da coordenada cilindrica equidistante 
ax = plt.axes(projection=ccrs.PlateCarree())
#ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())


# adiciona continente e bordas
shapefile = list(shpreader.Reader('/home/coqueiro/Downloads/br_unidades_da_federacao/BR_UF_2019.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.5)

ax.coastlines(resolution='50m', color='black', linewidth=1.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-105, -20, 5), ylocs=np.arange(-70, 15, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

i = 22
t = file_1.variables['t'][i,3,400:641,400:641]
pressao = file_2.variables['msl'][i,400:641,400:641]/100
c = 0.286
theta = t*(1000/800)**c

# temp
vmin= 270
vmax= 324
data_min = vmin
data_max = vmax
interval = 2
levels = np.arange(data_min,data_max,interval)

# pnmm
vmin2 = 930
vmax2 = 1028
 
data_min2 = vmin2
data_max2 = vmax2
interval2 = 2
levels2 = np.arange(data_min2,data_max2,interval2)

# vort relativa
vmin4= -22
vmax4= -2
data_min4 = vmin4
data_max4 = vmax4
interval4 = 2
levels4 = np.arange(data_min4,data_max4,interval4)

img = plt.contourf(lons,lats, theta, vmin=vmin, vmax=vmax,   cmap='jet', levels=levels, extend='both')
img2 = ax.contour(lons, lats, theta, vmin=vmin, vmax=vmax, colors='white', linewidths=0.3,  transform=ccrs.PlateCarree(), levels=levels)
img3 = ax.contour(lons, lats, pressao, vmin=vmin2, vmax=vmax2, colors='black', linewidths=0.8, levels=levels2)
ax.clabel(img3, inline=1, inline_spacing=3, fontsize=14, fmt = '%3.0f', colors= 'black')

# adiciona legenda 
cb = plt.colorbar(img, extend ='both', orientation = 'horizontal', pad=0.04, fraction=0.04)
font_size = 20 # Adjust as appropriate.
cb.ax.tick_params(labelsize=font_size)

# Getting the file time and date
add_seconds = int(file_1.variables['time'][i])
date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
date_formatted = date.strftime('%Y-%m-%d %H')
	
# Add a title
plt.title('Temperatura potencial (K) - 800 hPa', fontweight='bold', fontsize=35, loc='left')
plt.title(f'{date_formatted}', fontsize=15, loc='right')
#----------------------------------------------------------------------------------------------------------- 
# Salva imagem
plt.savefig(f'adv_temp {date_formatted}.png', bbox_inches='tight')