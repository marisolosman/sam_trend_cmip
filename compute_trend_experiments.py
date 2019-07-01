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
II = 0
VAR = 'psl'
MODEL = ['IPSL-CM6A-LR', 'CanESM5']
NEXPER=[4, 5]
ensemble = np.concatenate([np.repeat('H', 31), np.repeat('O3', 10), np.repeat('GHG', 10),
                           np.repeat('Nat', 10)])
#ensemble = np.concatenate([np.repeat('H', 25), np.repeat('O3', 10), np.repeat('GHG', 10),
#                           np.repeat('Nat', 10), np.repeat('Sol', 10)])

EXPERIMENT = ['Historical', 'HistoricalO3','HistoricalGHG', 'HistoricalNat', 'HistoricalSol']
PERIOD = ['185001-201412', '185001-202012',  '185001-202012',  '185001-202012', '185001-202012']
FREQ = 'mon'
LATMIN = -90
LATMAX = -20
TIMEMIN = '1951-01-01'
TIMEMAX = '2014-12-31'
#path to model data
ds_exp = []
sam = []
for i in range(NEXPER[II]):
    PATH = '/datos3/CMIP6/' + EXPERIMENT[i] + '/' + FREQ + '/' + VAR + '/'
    ARCHIVO = VAR + '_A' + FREQ + '_' + MODEL[II] + '_' + EXPERIMENT[i] +'_' + PERIOD[i] + '.nc4'
    FILE = glob.glob(PATH + ARCHIVO)
    #Open data with xarray
    ds = xr.open_dataset(FILE[0], decode_coords=False)
    #Select latitudinal, longitudinal and temporal range
    ds_sub = ds.sel(lat=slice(LATMIN, LATMAX), time=slice(TIMEMIN, TIMEMAX))
    #compute sam index: zonal mean normalized pressure between 40S and 65S
    ds1 = ds_sub.sel(lat=-40, method='nearest').mean(dim='lon')
    ds2 = ds_sub.sel(lat=-65, method='nearest').mean(dim='lon')
    sam_index = (ds1 - ds1.mean(dim='time')) / ds1.std(dim='time') - (ds2-ds2.mean(dim='time')) / ds2.std(dim='time')
    sam.append(sam_index)

RUTA = '/datos/osman/sam_trend_cmip_figures/'
sam = xr.concat(sam, dim='ensemble')
sam ['ensemble'] = ensemble
trends = np.zeros([5, len(ensemble)])
pval = np.zeros([5, len(ensemble)])
sam_seasonal = sam.resample(time='QS-MAR').mean(dim='time')
sam_annual = sam.groupby('time.year').mean(dim='time')

for j in range(len(ensemble)):
    for i in range(4):
        y = sam_seasonal.isel(ensemble=j).sel(time=sam_seasonal['time.month']== 3 + i * 3)
        [trends[i, j], interc, r_va, pval[i, j], z] = stats.linregress(np.arange(len(y.psl.values)),
                                                                       y.psl.values)
    y = sam_annual.isel(ensemble=j)
    [trends[4, j], interc, r_va, pval[4, j], z] = stats.linregress(np.arange(len(y.psl.values)),
                                                                   y.psl.values)
slopes = xr.Dataset({'slope': (('season', 'ensemble'), trends), 'p_value': (('season', 'ensemble'),
                                                                            pval)},
                    coords={'season': np.array(['MAM', 'JJA', 'SON', 'DJF', 'Annual']),
                            'ensemble': ensemble})
fig1 = plt.figure(1, (10, 4), dpi=500)  #fig size in inches
color = np.array( ['k', 'r', 'b', 'g', 'violet'])
simbolos = [(5,1), (3, 0), (4, 0), (5, 0), (4, 1)]
ens = np.array(['H', 'O3', 'GHG', 'Nat', 'Sol'])
fig, ax = plt.subplots()
for i in range(NEXPER[II]):
    colores = color[i]
    simb = simbolos[i]
    x = np.repeat(np.arange(5) + (1 + 0.10 * i), np.sum(slopes.ensemble.values==ens[i])).reshape([5, np.sum(slopes.ensemble.values==ens[i])])
    ax.scatter(x, slopes.slope.values[:, slopes.ensemble.values==ens[i]]*30,
              color=colores, alpha=0.5, label=ens[i], edgecolors='none')
ax.legend()
for i in range(NEXPER[II]):
    colores = color[i]#[slopes.ensemble.values[i] == ens]
    simb = simbolos[i]#[np.int(np.where(slopes.ensemble.values[i] == ens)[0])]
    x = np.arange(5) + (1 + 0.10 * i)
    ax.scatter(x, np.mean(slopes.slope.values[:, slopes.ensemble.values==ens[i]]*30, axis=1),
               s=np.repeat(180,5),  color=colores, marker='_')
    print(np.mean(slopes.slope.values[:, slopes.ensemble.values==ens[i]]*30, axis=1))
ax.axhline(y=0, xmin=0, xmax=6, linestyle='-',linewidth=0.5, color='k')
ax.set_ylabel('Trend (std units /30years)')
ax.set_xticks(np.arange(5) + 1.10)
ax.set_xticklabels(slopes.season.values)
ax.set_title(MODEL[II] + ' Southern Annular Mode 1957-2014 Trend', fontsize=12)
fig.savefig(RUTA + MODEL[II] + '_simulated_sam_trend_experiments.png', dpi=500, bbox_inches='tight',
             papertype='A4', orientation='landscape')

