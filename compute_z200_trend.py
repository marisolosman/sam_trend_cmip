#!/usr/bin/python
#import sys
import numpy as np
import xarray as xr
import glob
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import LG_xarray_tools
import mpl_toolkits.basemap as bm
#from cartopy.util import add_cyclic_point
VAR = 'zg200'
MODEL = 'IPSL-CM6A-LR'
EXPERIMENT = ['Historical', 'HistoricalO3','HistoricalGHG', 'HistoricalNat']#,'HistoricalSol']
FREQ = 'mon'
LATMIN = -90
LATMAX = -20
TIMEMIN = '1951-01-01'
TIMEMAX = '2014-12-31'
#path to model data
ds_exp = []
sam = []
for i in range(len(EXPERIMENT)):
    PATH = '/datos3/CMIP6/' + EXPERIMENT[i] + '/' + FREQ + '/' + VAR + '/'
    FILE = glob.glob(PATH + '*IPSL*.nc*')
    #VAR + '_A' + FREQ + '_' + MODEL + '_' + 'historical_185001-201412.nc'
    #Open data with xarray
    ds = xr.open_dataset(FILE[0], decode_coords=False)
    #Select latitudinal, longitudinal and temporal range
    ds_sub = ds.sel(lat=slice(LATMIN, LATMAX), time=slice(TIMEMIN, TIMEMAX)).mean(dim='ensemble')
    ds_exp.append(ds_sub)
ds_exp = xr.concat(ds_exp, dim='experiment')
ds_exp['experiment'] = EXPERIMENT
ds_exp = ds_exp.sel(time=slice('1979-12-01', '2014-03-31'))
ds_exp = ds_exp.resample(time='6M').mean(dim='time')
ds_exp = ds_exp.sel(time=ds_exp['time.month']==12)
#compute trend and plot
ds1_exp =ds_exp.stack(points=['experiment', 'lat', 'lon'])

trends = np.zeros_like(ds1_exp.zg[0, 0, :])
pval = np.zeros_like(ds1_exp.zg[0, 0, :])

for i in range(ds1_exp.zg.shape[2]):
    y = ds1_exp.zg[:, 0, i]
    [trends[i], interc, r_va, pval[i], z] = stats.linregress(np.arange(len(ds1_exp.time.values)), y)
trends = np.reshape(trends, (4, len(ds_exp['lat']), len(ds_exp['lon'])))
pval = np.reshape(pval, (4, len(ds_exp['lat']), len(ds_exp['lon'])))

slopes = xr.Dataset({'slope': (('experiment', 'lat', 'lon'), trends), 'p_value': (('experiment',
                                                                                      'lat', 'lon'),
                                                                                     pval)},
                      coords={'experiment': ds_exp['experiment'],
                              'lat': ds_exp['lat'], 'lon': ds_exp['lon']})
RUTA = '/datos/osman/sam_trend_cmip_figures/'
[dx, dy] = np.meshgrid(np.append(slopes.lon.values,360), slopes.lat.values)
fig = plt.figure(2,(18,7),300)
for i in range(4):
 #   ax = fig.add_subplot(1,5,i+1)
     plt.subplot(1,4,i+1)
     mapproj = bm.Basemap(projection='spstere',lon_0=120,lat_0=-90,
                           boundinglat=-30, resolution='l',round=True)
     mapproj.drawcoastlines()
     lonproj, latproj = mapproj(dx, dy)      #poject grid
     # set desired contour levels.
     clevs = np.linspace(-40, 40, 21)
     barra = plt.cm.PuOr_r #colorbar
     y = np.concatenate((slopes.slope.values[i, :, :], slopes.slope.values[i, :, 0][:, np.newaxis]), axis=1)
     CS1 = mapproj.contourf(lonproj, latproj, y * 30, clevs, cmap=barra,
                            extend='both', vmin=-40, vmax=40)
     barra.set_under(barra(0))
     barra.set_over(barra(barra.N-1))
     plt.title(slopes.experiment.values[i], fontsize=10)
     #titulo general
     plt.suptitle('DJFMAM Z 200hPa Trend (m/30year)', fontsize=12, x=0.52, y=0.85)
     fig.subplots_adjust(top=0.9, bottom=0.08)
cbar_ax = fig.add_axes([0.29, 0.15, 0.45, 0.05])
fig.colorbar(CS1, cax=cbar_ax,orientation='horizontal')
#cbar_ax.set_xticklabels([-0.9, -0.6,-0.3,0,0.3,0.6,0.9])#,size = 9)
plt.savefig(RUTA + 'djfmam_1979_present_zg200_trends.png', dpi=300, bbox_inches='tight', orientation='landscape', 
            papertype='A4')

