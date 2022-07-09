#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 23:14:33 2022

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

file_1 = xr.open_dataset('/home/coqueiro/ufrj/micro/multi_bomba.nc', decode_times = False)
file_2 = xr.open_dataset('/home/coqueiro/ufrj/micro/pnmm_bomba.nc', decode_times = False)


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

for i in range(len(file_1.variables['time'])):
    # np.array pra colocar a densidade linha senão dá error: AttributeError: 'Variable' object has no attribute 'ravel'
    u_comp = np.array(file_1.variables['u'][i,0,400:641,400:641])
    v_comp = np.array(file_1.variables['v'][i,0,400:641,400:641])
    
    # magnitude do vento
    mag = np.sqrt(u_comp**2+v_comp**2)
    
    
    #------------------------------------------------------------------------------
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(20,20))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())
    
    shapefile = list(shpreader.Reader('/home/coqueiro/Downloads/br_unidades_da_federacao/BR_UF_2019.shp').geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.5)
    
    
    # adiciona continente e bordas
    ax.coastlines(resolution='50m', color='black', linewidth=1.5)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 30), ylocs=np.arange(-90, 90, 10), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    
    
    #------------------------------------------------------------------------------
    vmin= 30
    vmax= 80
    data_min = vmin
    data_max = vmax
    interval = 10
    levels = np.arange(data_min,data_max,interval)
    
    # corrente de jato
    img = plt.contourf(lons,lats, mag, vmin=vmin, vmax=vmax, cmap=cmocean.cm.amp, levels = levels, extend='both')
    img2 = ax.contour(lons, lats, mag, colors='white', linewidths=0.3, levels=levels, transform=ccrs.PlateCarree())
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
    plt.title('Corrente de jato (m/s) - 250 hPa', fontweight='bold', fontsize=35, loc='left')
    plt.title(f'{date_formatted}', fontsize=15, loc='right')
    #----------------------------------------------------------------------------------------------------------- 
    # Salva imagem
    plt.savefig(f'jato_250 {date_formatted}.png', bbox_inches='tight')
'''
i=0
# np.array pra colocar a densidade linha senão dá error: AttributeError: 'Variable' object has no attribute 'ravel'
u_comp = np.array(file_1.variables['u'][i,2,300:641,300:641])
v_comp = np.array(file_1.variables['v'][i,2,300:641,300:641])

# magnitude do vento
mag = np.sqrt(u_comp**2+v_comp**2)


#------------------------------------------------------------------------------
# escolha o tamanho do plot em polegadas (largura x altura)
plt.figure(figsize=(20,20))

# usando a projeção da coordenada cilindrica equidistante 
ax = plt.axes(projection=ccrs.PlateCarree())
#ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())



# adiciona continente e bordas
ax.coastlines(resolution='50m', color='black', linewidth=1.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 30), ylocs=np.arange(-90, 90, 10), draw_labels=True)
gl.top_labels = False
gl.right_labels = False



#------------------------------------------------------------------------------
vmin= 30
vmax= 80
data_min = vmin
data_max = vmax
interval = 10
levels = np.arange(data_min,data_max,interval)

# corrente de jato
img = plt.contourf(lons,lats, mag, vmin=vmin, vmax=vmax, cmap=cmocean.cm.amp, levels = levels, extend='both')
img2 = ax.contour(lons, lats, mag, colors='white', linewidths=0.3, levels=levels, transform=ccrs.PlateCarree())
ax.streamplot(lons, lats, u_comp, v_comp, density=[6,6], linewidth=1.5, color='black', transform=ccrs.PlateCarree())

# adiciona legenda 
plt.colorbar(img, label = 'vento (m/s)', extend ='both', orientation = 'vertical', pad=0.04, fraction=0.04)

# Getting the file time and date
add_seconds = int(file_1.variables['time'][i])
date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
date_formatted = date.strftime('%Y-%m-%d %H')
	
# Add a title
plt.title('Corrente de jato (250 hPa)', fontweight='bold', fontsize=20, loc='left')
plt.title(f'{date_formatted}', fontsize=15, loc='right')
#----------------------------------------------------------------------------------------------------------- 
# Salva imagem
plt.savefig(f'jato_250 {date_formatted}.png')
'''