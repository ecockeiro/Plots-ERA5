#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 24 18:47:59 2022

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
import os

#dataset
input_dir = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Dados'
saida_dir = r'F:\Lucas\Conteudo\Fisica das nuvens e precipitacao\Figuras'
continental_multi = os.path.join(input_dir, '24_02_19_continental_multi.nc')
continental_single = os.path.join(input_dir, '24_02_19_continental_single.nc')
oceanico_multi = os.path.join(input_dir, '01_07_19_oceanico_multi.nc')
oceanico_single = os.path.join(input_dir, '01_07_19_oceanico_single.nc')

file_0 = xr.open_dataset(
    continental_multi,
    decode_times = False
    )

file_1 = xr.open_dataset(
    continental_multi
    ).metpy.parse_cf()

file_1 = file_1.assign_coords(dict(
    longitude = (((file_1.longitude.values + 180) % 360) - 180))
    ).sortby('longitude')

file_2 = xr.open_dataset(
    continental_single
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
level = 1000

# intevalos da umidade especifica
# define os intervalos da legenda
q_levs = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18])

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
    'Custom cmap', colors, q_levs.shape[0] - 1)
cmap.set_over(np.array([0, 37, 89, 255])/255)
cmap.set_under('white')

# nromaliza com base nos intervalos
norm = mcolors.BoundaryNorm(q_levs, cmap.N) # usa no PColormesh, nao no Contourf

# intevalos da temperatura
temp_min = -40
temp_max = 60
interval_3 = 2             # de quanto em quanto voce quer que varie
temp_levs = np.arange(temp_min, temp_max, interval_3)

for i in range(len(file_1.variables['time'])):
    args = dict(
        time = file_1.time[i],  
        vertical=level,
        latitude=lat_slice,
        longitude=lon_slice
    )
    
    t = file_1.t.metpy.sel(**args).metpy.unit_array.squeeze()
    t = t.to('degC')
    
    q = file_1.q.metpy.sel(**args).metpy.unit_array.squeeze()*1000
    
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
                      xlocs=np.arange(-180, 180, 5), 
                      ylocs=np.arange(-90, 90, 5), 
                      draw_labels=True
                      )
    gl.top_labels = False
    gl.right_labels = False
    
    # adiciona mascara de terra
    ax.add_feature(cfeature.LAND)
    
    # plota a imagem geopotencial
    sombreado = ax.contourf(lons, 
                            lats, 
                            q, 
                            cmap=cmap, 
                            levels = q_levs, 
                            extend = 'both'
                            )
    
    # plota a imagem pressao
    contorno = ax.contour(lons,
                          lats, 
                          t, 
                          colors='black', 
                          linewidths=0.8, 
                          levels=temp_levs
                          )
    
    ax.clabel(contorno, 
              inline = 1, 
              inline_spacing = 1, 
              fontsize=20, 
              fmt = '%3.0f', 
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
    plt.title('Temperatura (°C) e umidade (g/kg) - '+ str(level) +' hPa',
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
    fname = f'temperatura_umidade_{date_formatted}.png'
    plt.savefig(os.path.join(saida_dir, fname), bbox_inches='tight', dpi = 200)
    plt.close()