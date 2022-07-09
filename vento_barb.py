#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 21:52:36 2022

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
#------------------------------------------------------------------------------
# abre os arquivos netCDF4

file_1 = xr.open_dataset('/home/coqueiro/ufrj/micro/multi_bomba.nc', decode_times = False)
file_2 = xr.open_dataset('/home/coqueiro/ufrj/micro/pnmm_bomba.nc', decode_times = False)


# transformando as coordenadas de longitude dos arquivos
file_1 = file_1.assign_coords(dict(longitude = (((file_1.longitude.values + 180) % 360) - 180))).sortby('longitude')
file_2 = file_2.assign_coords(dict(longitude = (((file_2.longitude.values + 180) % 360) - 180))).sortby('longitude')


# limites de lat e lon
lat_min = -70.00
lat_max = 15.00
lon_min = -105.0
lon_max = -20.00

# seleciona uma extensao minima para o recorte de dados
extent = [lon_min, lon_max, lat_min, lat_max]

# ler as lats e lons
# velocidades
lats = file_1.variables['latitude'][300:641:9]
lons = file_1.variables['longitude'][300:641:9]

# define a extensao da imagem
img_extent = [lon_min, lon_max, lat_min, lat_max]


#------------------------------------------------------------------------------
# extrai a variavel
# >>> como setar os valores corretamente no dado nc <<<
# >>> ex: file_1.variables['var'][time[pos],level[pos],lat[pos_inicial]:lat[pos_final],lon[pos_inicial]:lon[pos_inicial]]
# var = variavel que voce quer extrair dê um print(file_1.variables) para ver todas as variaveis
# pos = posicao que tal latitude ocupa no dado .nc (180W = -180 = file_1['longitude'][0])
# variação de 10 graus = 40 posições

#i=3
for i in range(len(file_1.variables['time'])):
    u_comp = file_1.variables['u'][i,2,300:641:9,300:641:9]*1.94384
    v_comp = file_1.variables['v'][i,2,300:641:9,300:641:9]*1.94384
    #mag = (u_comp**2+v_comp**2)**0.5
    pressao = file_2.variables['msl'][i,300:641:9,300:641:9]/100
    #------------------------------------------------------------------------------
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(20,20))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())
    
    
    
    # adiciona continente e bordas
    ax.coastlines(resolution='50m', color='black', linewidth=2.5)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1.5)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='white', alpha=1, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 30), ylocs=np.arange(-90, 90, 10), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    
    #------------------------------------------------------------------------------
    # vento
    vmin= 40
    vmax= 120
    data_min = vmin
    data_max = vmax
    interval = 10
    levels = np.arange(data_min,data_max,interval)
    
    # pnmm
    vmin2 = 930
    vmax2 = 1028
    
    data_min2 = vmin2
    data_max2 = vmax2
    interval2 = 2
    levels2 = np.arange(data_min2,data_max2,interval2)
    
    img1 = ax.contourf(lons, lats, pressao, vmin= vmin2, vmax=vmax2, cmap=cmocean.cm.balance, levels=levels2, extend='both')
    img3 = ax.contour(lons, lats, pressao, vmin= vmin2, vmax=vmax2, colors='black', linewidths=0.8, levels=levels2)
    ax.clabel(img1, inline=1, inline_spacing=3, fontsize='14',fmt = '%1.0f', colors= 'white')
    
    
    
    # plot barbela
    #img = plt.contour(lons,lats,mag, vmin=vmin, vmax=vmax, levels = levels, extend='both')
    plt.barbs(
        lons,
        lats,
        u_comp,
        v_comp,
        fill_empty=True,
        length=7,
        sizes=dict(emptybarb=0, height=0.6),
        barbcolor="black",
        barb_increments=dict(flag=50),)
    
    
    
    # Add a title
    #valid = str(geopotencial.validDate)     # Valid date / time 
    plt.title('Vento' , fontweight='bold', fontsize=10, loc='left')
    #plt.title('Valid: ' + valid, fontsize=10, loc='right')
    #----------------------------------------------------------------------------------------------------------- 
    # Getting the file time and date
    add_seconds = int(file_1.variables['time'][i])
    date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
    date_formatted = date.strftime('%Y-%m-%d %H')
    	
    # Add a title
    plt.title('Pnmm e vento (850 hPa)', fontweight='bold', fontsize=20, loc='left')
    plt.title(f'{date_formatted}', fontsize=15, loc='right')
    
    # adiciona legenda 
    plt.colorbar(img1, label = 'pressão', extend ='both', orientation = 'vertical', pad=0.04, fraction=0.04)
    
    #----------------------------------------------------------------------------------------------------------- 
    # Salva imagem
    plt.savefig(f'pnmm_barb_850 {date_formatted}.png')
