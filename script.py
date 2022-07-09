#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 17:55:27 2022

@author: coqueiro
"""
#imports 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from metpy.units import units
import numpy as np
import xarray as xr
import cartopy, cartopy.crs as ccrs  # plot mapas
import cartopy.feature as cfeature  # operacao de desenho em imagens
from metpy.calc import temperature_from_potential_temperature
ds = xr.open_dataset('/home/coqueiro/ufrj/sinotica_2/dados_nc/yakecan_multi_16_21.nc').metpy.parse_cf()

# Set subset slice for the geographic extent of data to limit download
lon_slice = slice(300, 360)
lat_slice = slice(-15, -50)

# Grab lat/lon values (GFS will be 1D)
lats = ds.latitude.sel(latitude=lat_slice).values
lons = ds.longitude.sel(longitude=lon_slice).values

# Set level to plot/analyze
level = 850 

# Grad individual data arrays with units from our file, selecting for level and lat/lon slice
hght_850 = (ds['z'].metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).squeeze().metpy.unit_array)*units.hPa
tmpk_850 = ds['t'].metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).squeeze().metpy.unit_array
uwnd_850 = ds['u'].metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).squeeze().metpy.unit_array
vwnd_850 = ds['v'].metpy.sel(
    vertical=level, latitude=lat_slice, longitude=lon_slice).squeeze().metpy.unit_array

# Convert temperatures to degree Celsius for plotting purposes
tmpc_850 = tmpk_850.to('degC')

# Get a sensible datetime format
vtime = ds.time.data[0].astype('datetime64[ms]').astype('O')

# escolha o tamanho do plot em polegadas (largura x altura)
plt.figure(figsize=(20,20))
 
# usando a projeção da coordenada cilindrica equidistante 
ax = plt.axes(projection=ccrs.PlateCarree())
#ax.set_extent([extent[0], extent[2], extent[1], extent[3]], ccrs.PlateCarree())
 
datacrs = ccrs.PlateCarree() 
 
# adiciona continente e bordas
ax.coastlines(resolution='50m', color='black', linewidth=1.5)
ax.add_feature(cartopy.feature.BORDERS, edgecolor='black', linewidth=1)
gl = ax.gridlines(crs=ccrs.PlateCarree(), color='gray', alpha=1.0, linestyle='--', linewidth=0.25, xlocs=np.arange(-180, 180, 5), ylocs=np.arange(-90, 90, 5), draw_labels=True)
gl.top_labels = False
gl.right_labels = False
 

# Plot colorfill and dashed contours of 850-hPa temperatures in Celsius
clevs_850_tmpc = np.arange(-40, 41, 2)
cf = ax.contourf(lons, lats, tmpc_850, clevs_850_tmpc, cmap='jet', transform=datacrs)
cf = ax.contourf(lons, lats, tmpc_850, cmap='jet', transform=datacrs)
cb = plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50)
cb.set_label('Temperature (C)')
csf = ax.contour(lons, lats, tmpc_850, clevs_850_tmpc, colors='grey', linestyles='dashed', transform=datacrs)
plt.clabel(csf, fmt='%d')

# Plot contours of 850-hPa geopotential heights in meters
# clevs_850_hght = np.arange(0, 8000, 30)
# cs = ax.contour(lons, lats, hght_850, clevs_850_hght, colors='black', transform=datacrs)
# plt.clabel(cs, fmt='%d')

# Plot wind barbs every fifth element
wind_slice = (slice(None, None, 5), slice(None, None, 5))
ax.barbs(lons[wind_slice[0]], lats[wind_slice[1]],
         uwnd_850[wind_slice[0], wind_slice[1]].to('kt').m,
         vwnd_850[wind_slice[0], wind_slice[1]].to('kt').m,
         pivot='middle', color='black', transform=datacrs)

# Add some titles
plt.title('850-hPa GFS Geopotential Heights (m), Temperature (C), '
          'and Wind Barbs (kt)', loc='left')
plt.title('Valid Time: {}'.format(vtime), loc='right')

plt.show()

# Add coastline and state map features
ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
ax.add_feature(cfeature.STATES.with_scale('50m'))

# Plot colorfill and dashed contours of 850-hPa temperatures in Celsius
clevs_850_tmpc = np.arange(-40, 41, 2)
cf = ax.contourf(lons, lats, tmpc_850, clevs_850_tmpc, cmap='jet', transform=datacrs)
cb = plt.colorbar(cf, orientation='horizontal', pad=0, aspect=50)
cb.set_label('Temperature (C)')
csf = ax.contour(lons, lats, tmpc_850, clevs_850_tmpc, colors='grey',
                 linestyles='dashed', transform=datacrs)
plt.clabel(csf, fmt='%d')

# Plot contours of 850-hPa geopotential heights in meters
clevs_850_hght = np.arange(0, 8000, 30)
cs = ax.contour(lons, lats, hght_850, clevs_850_hght, colors='black', transform=datacrs)
plt.clabel(cs, fmt='%d')

# Plot wind barbs every fifth element
wind_slice = (slice(None, None, 5), slice(None, None, 5))
ax.barbs(lons[wind_slice[0]], lats[wind_slice[1]],
         uwnd_850[wind_slice[0], wind_slice[1]].to('kt').m,
         vwnd_850[wind_slice[0], wind_slice[1]].to('kt').m,
         pivot='middle', color='black', transform=datacrs)

# Add some titles
plt.title('850-hPa GFS Geopotential Heights (m), Temperature (C), '
          'and Wind Barbs (kt)', loc='left')
plt.title('Valid Time: {}'.format(vtime), loc='right')

plt.show()