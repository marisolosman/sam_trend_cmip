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

II = 1
RUTA = '/datos/osman/sam_trend_cmip_figures/'
VAR = 'psl'
MODEL = ['IPSL-CM6A-LR', 'CanESM5']
NEXPER=[4, 5]
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
    ARCHIVO = VAR + '_A' + FREQ + '_' + MODEL[1] + '_' + EXPERIMENT[i] +'_' + PERIOD[i] + '.nc4'
    FILE = glob.glob(PATH + ARCHIVO)
    #Open data with xarray
    ds = xr.open_dataset(FILE[0], decode_coords=False)
    #Select latitudinal, longitudinal and temporal range
    ds_sub = ds.sel(lat=slice(LATMIN, LATMAX), time=slice(TIMEMIN, TIMEMAX))
    #compute sam index: zonal mean normalized pressure between 40S and 65S
    ds1 = ds_sub.sel(lat=-40, method='nearest').mean(dim='lon')
    ds2 = ds_sub.sel(lat=-65, method='nearest').mean(dim='lon')
    sam_index = (ds1 - ds1.mean(dim='time')) / ds1.std(dim='time') - (ds2-ds2.mean(dim='time')) / ds2.std(dim='time')
    sam.append(sam_index.mean(dim='ensemble', skipna=True))

#====================================================================
sam = xr.concat(sam, dim='experiment')
sam['experiment'] = EXPERIMENT[0: NEXPER[II]]

seas = ['MAM', 'JJA', 'SON', 'DJF']
sam_seasonal = sam.resample(time='QS-MAR').mean(skipna=True)
sam_annual = sam.groupby('time.year').mean(dim='time')
for j in range(NEXPER[II]):
    [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(sam_annual.year.values)),
                                                       sam_annual.psl.isel(experiment=j).values)
    print("Annual")
    print(slope * 30, p_val)
    for i in range(4):
        print(seas[i])
        [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0,
                                                                     len(sam_seasonal.time.values[sam_seasonal['time.month'] == 3 +  i * 3])),
                                                           sam_seasonal.psl.sel(time=sam_seasonal['time.month']==
                                                                            3 + i *
                                                                            3).isel(experiment=j))
        print(slope*30, p_val)

#plot and compute trends
#plot sam time serie
#fig1 = plt.figure(1, (10, 4), dpi=500)  #fig size in inches
#color = iter(cm.rainbow(np.linspace(0, 1, 5)))
#for i in range(len(EXPERIMENT)):
#    c = next(color)
#    plt.plot(sam.time.values, sam.psl.isel(experiment=i).mean(dim='ensemble', skipna=True).values,
#
#    plt.plot(sam.time.values, sam.psl.isel(experiment=i).values,
#             color=c, linewidth=1, label=sam.experiment.values[i])
#plt.axhline(y=0,xmin=0,linestyle='--',linewidth=0.8,
#           color='k',alpha=0.5)
#plt.ylim((-3.2, 3.2))
#plt.ylabel('SAM')
#plt.legend(bbox_to_anchor=(1.05, 0), loc='lower left', borderaxespad=0.)
#plt.title(MODEL + ' Southern Annular Mode 1957-2014',fontsize=12)
#fig1.savefig(RUTA + MODEL + '_simulated_sam_time_series.png', dpi=500, bbox_inches='tight',
#             papertype='A4', orientation='landscape')
