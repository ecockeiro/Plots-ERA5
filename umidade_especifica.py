#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 23:14:33 2022

@author: coqueiro
"""

from netCDF4 import Dataset  # ler e escreve netcdf ## uma merda essa biblioteca usa xarray mesmo ##
import xarray as xr # ler netCDF e muito mais <3
import matplotlib.pyplot as plt  # biblioteca de plots
import matplotlib.colors as mcolors                  # Matplotlib colors  
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import cartopy, cartopy.crs as ccrs  # plot mapas
import cartopy.feature as cfeature  # operacao de desenho em imagens
import numpy as np  # importa o pacote numpy
import cmocean 
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


#------------------------------------------------------------------------------
# extrai a variavel
# >>> como setar os valores corretamente no dado nc <<<
# >>> ex: file_1.variables['var'][time[pos],level[pos],lat[pos_inicial]:lat[pos_final],lon[pos_inicial]:lon[pos_inicial]]
# var = variavel que voce quer extrair dê um print(file_1.variables) para ver todas as variaveis
# pos = posicao que tal latitude ocupa no dado .nc (180W = -180 = file_1['longitude'][0])
# variação de 10 graus = 40 posições

i=22
# np.array pra colocar a densidade linha senão dá error: AttributeError: 'Variable' object has no attribute 'ravel'
u_comp = np.array(file_1.variables['u'][i,4,400:641,400:641])
v_comp = np.array(file_1.variables['v'][i,4,400:641,400:641])
q = file_1.variables['q'][i,4,400:641,400:641]*1000

# magnitude do vento
mag = np.sqrt(u_comp**2+v_comp**2)


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



#------------------------------------------------------------------------------
# define os intervalos da legenda
clevs = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18])

# lista de cores, em ordem crescete. RGBA
colors = np.array([ # (R, G, B, A)
    [242, 98, 0, 255],
    [249, 155, 77, 255],
    [254, 217, 118, 255],
    [255, 247, 188, 255],
    [190, 220, 230, 255],
    [156, 194, 255, 255],
    [59, 118, 255, 255],
    [0, 77, 182, 255]
]) / 255 # divide por 255

# cria um novo cmap a partir do pre-existente
cmap = mcolors.LinearSegmentedColormap.from_list(
    'Custom cmap', colors, clevs.shape[0] - 1)
cmap.set_over(np.array([0, 37, 89, 255])/255)
cmap.set_under('white')

# nromaliza com base nos intervalos
norm = mcolors.BoundaryNorm(clevs, cmap.N)

# corrente de jato
img = plt.contourf(lons,lats, q, cmap = cmap, levels = clevs, extend='both')
img2 = ax.contour(lons, lats, q, colors='white', linewidths=0.3,  transform=ccrs.PlateCarree())
ax.streamplot(lons, lats, u_comp, v_comp, density=[6,6], linewidth=1.5, color='black', transform=ccrs.PlateCarree())

# adiciona legenda 
cb = plt.colorbar(img, extend ='both', orientation = 'horizontal', pad=0.04, fraction=0.04)
font_size = 20 # Adjust as appropriate.
cb.ax.tick_params(labelsize=font_size)

# Getting the file time and date
add_seconds = int(file_1.variables['time'][i])
date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
date_formatted = date.strftime('%Y-%m-%d %H')
	
# Add a title
plt.title('Umidade específica (g/kg) - 850 hPa', fontweight='bold', fontsize=35, loc='left')
plt.title(f'{date_formatted}', fontsize=15, loc='right')
#----------------------------------------------------------------------------------------------------------- 

# Salva imagem
plt.savefig(f'umidade {date_formatted}.png', bbox_inches='tight')