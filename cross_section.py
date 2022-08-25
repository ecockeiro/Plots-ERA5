# Copyright (c) 2018 MetPy Developers.
# Distributed under the terms of the BSD 3-Clause License.
# SPDX-License-Identifier: BSD-3-Clause
"""
======================
Cross Section Analysis
======================

The MetPy function `metpy.interpolate.cross_section` can obtain a cross-sectional slice through
gridded data.
"""

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.interpolate import cross_section

##############################
# **Getting the data**
#
# This example uses [NARR reanalysis data](
# https://www.ncei.noaa.gov/products/weather-climate-models/north-american-regional)
# for 18 UTC 04 April 1987 from NCEI.
#
# We use MetPy's CF parsing to get the data ready for use, and squeeze down the size-one time
# dimension.

data = xr.open_dataset('/home/coqueiro/ufrj/fis_nuvens/dados/fdn_multi_hora_18.nc')
data = data.metpy.parse_cf().squeeze()
print(data)

##############################
# Define start and end points:

start = (-27., -56.)
end = (-38., -37.)

##############################
# Get the cross section, and convert lat/lon to supplementary coordinates:

cross = cross_section(data, start, end).set_coords(('latitude', 'longitude'))
print(cross)

##############################
# For this example, we will be plotting potential temperature, relative humidity, and
# tangential/normal winds. And so, we need to calculate those, and add them to the dataset:

cross['Potential_temperature'] = mpcalc.potential_temperature(
    cross['level'],
    cross['t']
)
cross['Relative_humidity'] = mpcalc.relative_humidity_from_specific_humidity(
    cross['level'],
    cross['t'],
    cross['q']
)
cross['u_wind'] = cross['u'].metpy.convert_units('knots')
cross['v_wind'] = cross['v'].metpy.convert_units('knots')
cross['t_wind'], cross['n_wind'] = mpcalc.cross_section_components(
    cross['u_wind'],
    cross['v_wind']
)

print(cross)

##############################
# Now, we can make the plot.

# Define the figure object and primary axes
fig = plt.figure(1, figsize=(16., 9.))
ax = plt.axes()

# Plot RH using contourf
rh_contour = ax.contourf(cross['longitude'], cross['level'], cross['Relative_humidity'],
                         levels=np.arange(0, 1.05, .05), cmap='YlGnBu')

rh_colorbar = fig.colorbar(rh_contour)

# Plot potential temperature using contour, with some custom labeling
theta_contour = ax.contour(cross['longitude'], cross['level'], cross['Potential_temperature'],
                           levels=np.arange(250, 450, 5), colors='k', linewidths=2)
theta_contour.clabel(theta_contour.levels[1::2], fontsize=8, colors='k', inline=1,
                     inline_spacing=8, fmt='%i', rightside_up=True, use_clabeltext=True)

# Plot winds using the axes interface directly, with some custom indexing to make the barbs
# less crowded
wind_slc_vert = list(range(0, 4, 1))
wind_slc_horz = slice(5, 100, 5)
ax.barbs(cross['longitude'][wind_slc_horz], cross['level'][wind_slc_vert],
         cross['t_wind'][wind_slc_vert, wind_slc_horz],
         cross['n_wind'][wind_slc_vert, wind_slc_horz], color='k')

# Adjust the y-axis to be logarithmic
ax.set_yscale('symlog')
ax.set_yticklabels(np.arange(1000, 50, -100))
ax.set_ylim(cross['level'].max(), cross['level'].min())
ax.set_yticks(np.arange(1000, 50, -100))

# Define the CRS and inset axes
data_crs = data['z'].metpy.cartopy_crs
ax_inset = fig.add_axes([0.097, 0.63, 0.25, 0.25], projection=data_crs)

# Plot geopotential height at 500 hPa using xarray's contour wrapper
ax_inset.contour(data['longitude'], data['latitude'], data['z'].sel(level=500.)/10,
                 levels=np.arange(4500, 6000, 100), cmap='inferno')

# Plot the path of the cross section
endpoints = data_crs.transform_points(ccrs.PlateCarree(),
                                      *np.vstack([start, end]).transpose()[::-1])
ax_inset.scatter(endpoints[:, 0], endpoints[:, 1], c='k', zorder=2)
ax_inset.plot(cross['longitude'], cross['latitude'], c='k', zorder=2)

# Add geographic features
ax_inset.coastlines()
ax_inset.add_feature(cfeature.STATES.with_scale('50m'), edgecolor='k', alpha=0.2, zorder=0)

# Set the titles and axes labels
ax_inset.set_title('')
ax.set_title(f'ERA5 Cross-Section \u2013 {start} to {end} \u2013 '
             f'Valid: {cross["time"].dt.strftime("%Y-%m-%d %H:%MZ").item()}\n'
             'Potential Temperature (K), Tangential/Normal Winds (knots), Relative Humidity '
             '(dimensionless)\nInset: Cross-Section Path and 500 hPa Geopotential Height')
ax.set_ylabel('Pressure (hPa)')
ax.set_xlabel('Longitude (degrees east)')
rh_colorbar.set_label('Relative Humidity (dimensionless)')

plt.show()
#plt.savefig(f'div_umi {date_formatted}.png', bbox_inches='tight')