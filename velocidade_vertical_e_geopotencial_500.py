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
import matplotlib.colors as mcolors
import metpy.calc as mpcalc
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy.io.shapereader as shpreader # Import shapefiles
from datetime import datetime, timedelta  # basicas datas e tipos de tempo
import os

#dataset
input_dir = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Dados'
saida_dir = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Figuras'
continental_multi = os.path.join(input_dir, '24_02_19_continental_multi.nc')
continental_single = os.path.join(input_dir, '24_02_19_continental_single.nc')
oceanico_multi = os.path.join(input_dir, '01_07_19_oceanico_multi.nc')
oceanico_single = os.path.join(input_dir, '01_07_19_oceanico_single.nc')


file_0 = xr.open_dataset(
    oceanico_multi,
    decode_times = False
    )

file_1 = xr.open_dataset(
    oceanico_multi
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

file_2 = xr.open_dataset(
    oceanico_single
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
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)
g = 9.8 # m/s^2

# colorbar da velocidade vertical
# intevalos da velocidade vertical
w_min = -14
w_max = 14
n_levs = 8             # de quanto em quanto voce quer que varie
wlevs = np.array([-12, -10, -8, -6, -4, -2, 2, 4, 6, 8, 10, 12])

# lista de cores, em ordem crescete. RGBA
colors = np.array([ # [R, G, B, A]
    [88, 0, 36, 255],
    [189, 96, 78, 255],
    [255, 255 ,255, 255],
    [90, 142, 191, 255],
    [24, 43, 90, 255]
]) / 255

# cria um novo cmap a partir do pre-existente
cmap = mcolors.LinearSegmentedColormap.from_list(
    'Custom cmap', colors, wlevs.shape[0] - 1)
cmap.set_over(colors[-1])
cmap.set_under(colors[0])

# nromaliza com base nos intervalos
norm = mcolors.BoundaryNorm(wlevs, cmap.N) # usa no PColormesh, nao no Contourf

for i in range(len(file_1.variables['time'])):
    args = dict(
        time = file_1.time[i],
        vertical=level,
        latitude=lat_slice,
        longitude=lon_slice
    )
    geopotencial = file_1.z.metpy.sel(**args).metpy.unit_array.squeeze()/g
    w = file_1.w.metpy.sel(**args).metpy.unit_array.squeeze()*10

    # escolha o tamanho do plot em polegadas (largura x altura)
    plt.figure(figsize=(25,25))
    
    # usando a projeção da coordenada cilindrica equidistante 
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    
    shapefile = list(
        shpreader.Reader(
        os.path.join(input_dir, 'BR_UF_2019.shp')
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
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # plota a imagem geopotencial
    sombreado = ax.contourf(lons, 
                            lats, 
                            w, 
                            cmap = cmap, 
                            levels = wlevs, 
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
    barra_de_cores.ax.set_xticks(wlevs)
    
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
    fname = f'vel_vert_geo {date_formatted}.png'
    plt.savefig(os.path.join(saida_dir, fname), bbox_inches='tight')
    
