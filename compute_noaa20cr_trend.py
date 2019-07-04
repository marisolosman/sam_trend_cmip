#!/usr/bin/python
#import sys
import numpy as np
import xarray as xr
import glob
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.basemap as bm
#vector de variables: 0 for zg200 1 for slp 2 for t2m
INDEX = 2
VAR = ['zg200', 'slp', 't2m']
var = ['hgt', 'prmsl', 'air']
#contourf levels for each variable
NIVELES = [np.linspace(-60, 60, 21), np.linspace(-2, 2, 21), np.linspace(-2, 2, 21)]
UNITS = ['m/30year', 'hPa/30year', 'K/30year']
PATH = '/datos3/noaa20cr/'
LATMIN = -20
LATMAX = -90
TIMEMIN = '1957-01-01'
TIMEMAX = '2014-12-31'
ds = xr.open_dataset(PATH + VAR[INDEX] + '_noaa20cr.nc')
if i == 0:
    ds_sub = ds.sel(lat=slice(LATMIN, LATMAX), time=slice(TIMEMIN, TIMEMAX), level=200)
else:
    ds_sub = ds.sel(lat=slice(LATMIN, LATMAX), time=slice(TIMEMIN, TIMEMAX))

ds_exp = ds_exp.resample(time='6M').mean(dim='time')
#select DJFMAM
ds_exp = ds_exp.sel(time=ds_exp['time.month']==12)
#compute trend and plot
ds1_exp =ds_exp.stack(points=['lat', 'lon'])
trends = np.zeros_like(ds1_exp[var[INDEX]][0, :])
pval = np.zeros_like(ds1_exp[var[INDEX]][0, :])
for i in range(ds1_exp[var[INDEX]].shape[1]):
    y = ds1_exp[var[INDEX]][:, i]
    [trends[i], interc, r_va, pval[i], z] = stats.linregress(np.arange(len(ds1_exp.time.values)), y)
trends = np.reshape(trends, (len(ds_exp['lat']), len(ds_exp['lon'])))
pval = np.reshape(pval, (len(ds_exp['lat']), len(ds_exp['lon'])))
slopes = xr.Dataset({'slope': (('lat', 'lon'), trends), 'p_value': (('lat', 'lon'),
                                                                    pval)},
                      coords={'lat': ds_exp['lat'], 'lon': ds_exp['lon']})
RUTA = '/datos/osman/sam_trend_cmip_figures/'
[dx, dy] = np.meshgrid(np.append(slopes.lon.values,360), slopes.lat.values)
fig = plt.figure(2,(7, 7),300)
mapproj = bm.Basemap(projection='spstere',lon_0=120,lat_0=-90,
                     boundinglat=-30, resolution='l', round=True)
mapproj.drawcoastlines()
lonproj, latproj = mapproj(dx, dy)      #poject grid
# set desired contour levels.
clevs = NIVELES[INDEX]
barra = plt.cm.PuOr_r #colorbar
y = np.concatenate((slopes.slope.values[:, :], slopes.slope.values[:, 0][:, np.newaxis]), axis=1)
CS1 = mapproj.contourf(lonproj, latproj, y * 30, clevs, cmap=barra,
                        extend='both', vmin=clevs[0], vmax=clevs[-1])
barra.set_under(barra(0))
barra.set_over(barra(barra.N-1))
#titulo general
plt.title('DJFMAM ' + VAR[INDEX] + ' Trend ' + UNITS[INDEX], fontsize=12)
fig.colorbar(CS1, orientation='horizontal')
plt.savefig(RUTA + 'djfmam_' + VAR[INDEX] + '_noaa20cr_trends.png', dpi=300,
            bbox_inches='tight', orientation='landscape', papertype='A4')


