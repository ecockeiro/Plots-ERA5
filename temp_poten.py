#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 18:18:16 2022

@author: coqueiro
"""

from netCDF4 import Dataset  # ler e escreve netcdf ## uma merda essa biblioteca usa xarray mesmo ##
import xarray as xr # ler netCDF e muito mais <3
import matplotlib.pyplot as plt  # biblioteca de plots
import matplotlib.colors                  # Matplotlib colors  
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cartopy, cartopy.crs as ccrs  # plot mapas
import cartopy.feature as cfeature  # operacao de desenho em imagens
import numpy as np  # importa o pacote numpy
import cmocean 
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader # Import shapefiles
#------------------------------------------------------------------------------
# abre os arquivos netCDF4

file_1 = xr.open_dataset('/home/coqueiro/ufrj/fis_nuvens/dados/fdn_multi.nc', decode_times = False)
file_2 = xr.open_dataset('/home/coqueiro/ufrj/fis_nuvens/dados/fdn_pnmm.nc', decode_times = False)

# transformando as coordenadas de longitude dos arquivos
file_1 = file_1.assign_coords(dict(longitude = (((file_1.longitude.values + 180) % 360) - 180))).sortby('longitude')
file_2 = file_2.assign_coords(dict(longitude = (((file_2.longitude.values + 180) % 360
) - 180))).sortby('longitude')


# limites de lat e lon
lat_min = -70.00
lat_max = 10.00
lon_min = -105.0
lon_max = -20.00

# seleciona uma extensao minima para o recorte de dados
extent = [lon_min, lon_max, lat_min, lat_max]

#geopotencial
lats= file_1.variables['latitude'][400:641]
lons= file_1.variables['longitude'][400:641]

#pressao
lats= file_2.variables['latitude'][400:641]
lons= file_2.variables['longitude'][400:641]

# define a extensao da imagem
img_extent = [lon_min, lon_max, lat_min, lat_max]
#------------------------------------------------------------------------------
# escolha o tamanho do plot em polegadas (largura x altura)
plt.figure(figsize=(20,20))

# usando a projeção da coordenada cilindrica equidistante 
*ax = plt.axes(projection=ccrs.PlateCarree())
#ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())



# adiciona continente e bordas
shapefile = list(shpreader.Reader('/home/coqueiro/Downloads/br_unidades_da_federacao/BR_UF_2019.shp').geometries())
ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.5)

ax.coastlines(resolution='50m', color='black', linewidth=1.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-105, -20, 5), ylocs=np.arange(-70, 15, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False


#------------------------------------------------------------------------------
# extrai a variavel
# >>> como setar os valores corretamente no dado nc <<<
# >>> ex: file_1.variables['var'][time[pos],level[pos],lat[pos_inicial]:lat[pos_final],lon[pos_inicial]:lon[pos_inicial]]
# var = variavel que voce quer extrair dê um print(file_1.variables) para ver todas as variaveis
# pos = posicao que tal latitude ocupa no dado .nc (180W = -180 = file_1['longitude'][0])
# variação de 10 graus = 40 posições

i=0
T = file_1.variables['t'][i,2,400:641,400:641]

Po = 1000
P = file_2.variables['msl'][i,400:641,400:641]/100
Rcp = 0.286

theta = (T*(Po/P)**Rcp)-273.15

vmin= -50
vmax= 70
data_min = vmin
data_max = vmax
interval = 10
levels = np.arange(data_min,data_max,interval)

#plt.imshow(theta, cmap='jet', extent=img_extent )
img = plt.contourf(lons,lats, theta,  cmap='seismic', levels = levels , extend='both')
img2 = ax.contour(lons, lats, theta, colors='white', linewidths=0.3,  transform=ccrs.PlateCarree())