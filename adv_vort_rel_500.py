#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 10:00:35 2022

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
from cartopy.util import add_cyclic_point
import mygrads as mg #pip install mygrads
import metpy.calc as mpcalc
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader # Import shapefiles
 
#------------------------------------------------------------------------------
# abre os arquivos netCDF4

file_1 = xr.open_dataset('/home/coqueiro/ufrj/micro/multi_bomba.nc', decode_times = False)
file_2 = xr.open_dataset('/home/coqueiro/ufrj/micro/pnmm_bomba.nc', decode_times = False)

# transformando as coordenadas de longitude dos arquivos
file_1 = file_1.assign_coords(dict(longitude = (((file_1.longitude.values + 180) % 360) - 180))).sortby('longitude')
file_2 = file_2.assign_coords(dict(longitude = (((file_2.longitude.values + 180) % 360) - 180))).sortby('longitude')


# limites de lat e lon
lat_min = -70.00
lat_max = 10.00
lon_min = -105.0
lon_max = -20.00

# seleciona uma extensao minima para o recorte de dados
extent = [lon_min, lon_max, lat_min, lat_max]

# ler as lats e lons
# velocidades
lats = file_1.variables['latitude'][400:641]
lons = file_1.variables['longitude'][400:641]

# define a extensao da imagem
img_extent = [lon_min, lon_max, lat_min, lat_max]

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
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 30), ylocs=np.arange(-90, 90, 10), draw_labels=True)
gl.top_labels = False
gl.right_labels = False

#adiciona mascara de terra
ax.add_feature(cfeature.LAND)

#------------------------------------------------------------------------------
# extrai a variavel
# >>> como setar os valores corretamente no dado nc <<<
# >>> ex: file_1.variables['var'][time[pos],level[pos],lat[pos_inicial]:lat[pos_final],lon[pos_inicial]:lon[pos_inicial]]
# var = variavel que voce quer extrair dê um print(file_1.variables) para ver todas as variaveis
# pos = posicao que tal latitude ocupa no dado .nc (180W = -180 = file_1['longitude'][0])
# variação de 10 graus = 40 posições

# np.array pra colocar a densidade linha senão dá error: AttributeError: 'Variable' object has no attribute 'ravel'
#3 a 6
i = 22
# pressao
pressao = file_2.variables['msl'][i,400:641,400:641]/100

# componente do geopotencial
geopotencial = file_1.variables['z'][i,1,400:641,400:641]/10

# 500 hPa
vmin1 = 4950
vmax1 = 5850

data_min1 = vmin1
data_max1 = vmax1
interval1 = 25
levels1 = np.arange(data_min1,data_max1,interval1)


# componente do vento
u_comp = np.array(file_1.variables['u'][i,1,400:641,400:641])
v_comp = np.array(file_1.variables['v'][i,1,400:641,400:641])
rel_vort = mg.hcurl(u_comp, v_comp, lats, lons)
adv_rel_vort = mg.hadv(u_comp, v_comp, rel_vort, lats, lons)*10000000
rel_vort = 100000*rel_vort

# geopotencial 500
vmin1 = 4900
vmax1 = 5900

data_min1 = vmin1
data_max1 = vmax1
interval1 = 100
levels1 = np.arange(data_min1,data_max1,interval1)


# adv vort relativa
vmin2= -1.5
vmax2= -0.1
data_min2 = vmin2
data_max2 = vmax2
interval2 = 0.10
levels2 = np.arange(data_min2,data_max2,interval2)

# pnmm
vmin3 = 930
vmax3 = 1028
data_min3 = vmin3
data_max3 = vmax3
interval3 = 4
levels3 = np.arange(data_min3,data_max3,interval3)

# vort relativa
vmin4= -22
vmax4= -2
data_min4 = vmin4
data_max4 = vmax4
interval4 = 2
levels4 = np.arange(data_min4,data_max4,interval4)

# plot adv. vorticidade relativa 500 hPa
img = plt.contourf(lons,lats, adv_rel_vort, vmin=vmin2, vmax=vmax2, cmap='gist_heat', levels = levels2, extend='both')
img2 = ax.contour(lons, lats, geopotencial, colors='green', linewidths=1, levels=levels1)
ax.clabel(img2, inline=1, inline_spacing=3, fontsize='14',fmt = '%1.0f', colors= 'black')

img4 = ax.contour(lons, lats, pressao, colors='black', linewidths=1, linestyles='-.', levels=levels3)
ax.clabel(img4, inline=1, inline_spacing=3, fontsize='14',fmt = '%1.0f', colors= 'black')

# adiciona legenda 
cb = plt.colorbar(img, extend ='both', orientation = 'horizontal', pad=0.04, fraction=0.04)
font_size = 20 # Adjust as appropriate.
cb.ax.tick_params(labelsize=font_size)

# Getting the file time and date
add_seconds = int(file_1.variables['time'][i])
date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
date_formatted = date.strftime('%Y-%m-%d %H')
	
# Add a title
plt.title('Adv. vort rel (10⁷.s⁻¹)/PNMM/Geopotencial - 500 hPa', fontweight='bold', fontsize=30, loc='left')
plt.title(f'{date_formatted}', fontsize=15, loc='right')


#----------------------------------------------------------------------------------------------------------- 
# Salva imagem
plt.savefig(f'adv_vort_rel {date_formatted}.png', bbox_inches='tight')