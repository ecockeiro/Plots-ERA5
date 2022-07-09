#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 18:52:45 2022

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
#geopotencial
lats = file_1.variables['latitude'][400:641]
lons = file_1.variables['longitude'][400:641]

#pressao
lats = file_2.variables['latitude'][400:641]
lons = file_2.variables['longitude'][400:641]

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
    geopotencial = file_1.variables['z'][i,1,400:641,400:641]/10
    pressao = file_2.variables['msl'][i,400:641,400:641]/100
    
    #------------------------------------------------------------------------------ 
    
    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    #ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())
    
    shapefile = list(shpreader.Reader('/home/coqueiro/Downloads/br_unidades_da_federacao/BR_UF_2019.shp').geometries())
    ax.add_geometries(shapefile, ccrs.PlateCarree(), edgecolor='black', facecolor='none', linewidth=0.5)
    
    # adiciona continente e bordas
    ax.coastlines(resolution='50m', color='black', linewidth=1)
    ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1)
    gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.5, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
    gl.top_labels = False
    gl.right_labels = False
    
    # cria uma escala de cores:
    colors = ["#2d001c", "#5b0351", "#780777", "#480a5e", "#1e1552", 
              "#1f337d", "#214c9f", "#2776c6", "#2fa5f1", "#1bad1d", 
              "#8ad900", "#ffec00", "#ffab00", "#f46300", "#de3b00", 
              "#ab1900", "#6b0200", '#3c0000']
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    cmap.set_over('#3c0000')
    cmap.set_under('#28000a')
    
    # intevalos entre isobaras
    # geopotencial 500
    vmin1 = 4900
    vmax1 = 6000
    
    data_min1 = vmin1
    data_max1 = vmax1
    interval1 = 50
    levels1 = np.arange(data_min1,data_max1,interval1)
    
    # pnmm
    vmin2 = 930
    vmax2 = 1028
    
    data_min2 = vmin2
    data_max2 = vmax2
    interval2 = 2
    levels2 = np.arange(data_min2,data_max2,interval2)
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    #------------------------------------------------------------------------------

    # plota a imagem geopotencial
    img1 = ax.contourf(lons, lats, geopotencial,  cmap=cmap, levels=levels1, extend='both')
    img2 = ax.contour(lons, lats, geopotencial, colors='white', linewidths=0.3, levels=levels1)
    
    # plota a imagem pressao
    img3 = ax.contour(lons, lats, pressao, colors='black', linewidths=0.8, levels=levels2)
    ax.clabel(img3, inline=1, inline_spacing=3, fontsize=14, fmt = '%1.0f', colors= 'white')
    
    
    # Getting the file time and date
    add_seconds = int(file_1.variables['time'][i])
    date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
    date_formatted = date.strftime('%Y-%m-%d %H')
    	
    # Add a title
    plt.title('Pnmm/Geopotencial (m) - 500 hPa', fontweight='bold', fontsize=35, loc='left')
    plt.title(f'{date_formatted}', fontsize=15, loc='right')
    
    # adiciona legenda 
    cb = plt.colorbar(img1, orientation = 'horizontal', pad=0.04, fraction=0.04)
    font_size = 20 # Adjust as appropriate.
    cb.ax.tick_params(labelsize=font_size)
    
    #----------------------------------------------------------------------------------------------------------- 
    # Salva imagem
    plt.savefig(f'pnmm_geo_500 {date_formatted}.png', bbox_inches='tight')


