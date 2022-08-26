#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 19:04:43 2022

@author: coqueiro
"""
#importando bibliotecas
from doctest import OutputChecker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors # para a colorbar
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
level = 1000


# CRIANDO O CMAP E NORM PARA A COLORBAR

# intevalos da pnmm
intervalo_min2 = 900
intervalo_max2 = 1040
interval_2 = 2              # de quanto em quanto voce quer que varie
levels_2 = np.arange(intervalo_min2, intervalo_max2, interval_2)

# intevalos da divergencia - umidade
divq_min = -2.4
divq_max = 2.4
n_levs = 16 # numero de intervalos
divlevs = np.round(np.linspace(divq_min, divq_max, n_levs), 1)

# lista de cores, em ordem crescete. RGBA
colors = ['mediumseagreen', 'mediumaquamarine', 'palegreen', 'white', 'wheat', 'gold', 'goldenrod']

# cria um novo cmap a partir do pre-existente
cmap = mcolors.LinearSegmentedColormap.from_list(
    'Custom cmap', colors, divlevs.shape[0] - 1)
cmap.set_over('darkgoldenrod')
cmap.set_under("seagreen")

# nromaliza com base nos intervalos
norm = mcolors.BoundaryNorm(divlevs, cmap.N) # util para o PCOLORMESH, CONTOURF nao usa

# variaveis repetidas em cada loop
dx, dy = mpcalc.lat_lon_grid_deltas(lons, lats)


for i in range(len(file_1.variables['time'])):
    args = dict(
        time = file_1.time[i] ,
        vertical=level,
        latitude=lat_slice,
        longitude=lon_slice
    )
    geopotencial = file_1.z.metpy.sel(**args).metpy.unit_array.squeeze()
    u = file_1.u.metpy.sel(**args).metpy.unit_array.squeeze()
    v = file_1.v.metpy.sel(**args).metpy.unit_array.squeeze()
    q = file_1.q.metpy.sel(**args).metpy.unit_array.squeeze()

    del args["vertical"]
    pnmm = file_2.msl.metpy.sel(**args).metpy.unit_array.squeeze()* 0.01 * units.hPa/units.Pa
    
    divergencia = mpcalc.divergence(u, v, dx=dx, dy=dy, x_dim=- 1, y_dim=- 2)
    divergencia_umidade = divergencia * q * 1e6
    
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
                            divergencia_umidade, 
                            cmap = cmap, 
                            levels = divlevs, 
                            extend = 'both'
                            )
    
    # plota a imagem pressao
    contorno = ax.contour(lons,
                          lats, 
                          pnmm, 
                          colors='black', 
                          linewidths=0.8, 
                          levels=levels_2
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
    barra_de_cores.ax.set_xticks(divlevs)
    
    # Getting the file time and date
    add_seconds = int(file_0.variables['time'][i])
    
    date = datetime(1900,1,1,0) + timedelta(hours=add_seconds)
    date_formatted = date.strftime('%Y-%m-%d %H')
    	
    # Add a title
    plt.title('Divergencia de umidade (1/s)*10^6 - '+ str(level) +' hPa',
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
    fname = f'div_umi {date_formatted}.png'
    plt.savefig(os.path.join(saida_dir, fname), bbox_inches='tight')
    plt.close()
