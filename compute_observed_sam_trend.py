#!/usr/bin/python
import numpy as np
from scipy import stats
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
FILE = '/datos/osman/newsam.1957.2019.txt'
#f = open(FILE, 'r')
#header1 = f.readline()
data = []
year = []
with open(FILE) as f:
        f.readline()
        f.readline()
        for line in f:
            columns = line.split()
            year.append(np.float(columns[0]))
            data.append(np.asarray([float(i) for i in columns[1:]]))

data = np.concatenate(data, axis=0)

time = pd.date_range('1957-01-15', freq='M', periods=np.shape(data)[0])

sam = xr.DataArray(data, coords=[time], dims=['time'])

#selecciono periodo comun con modelos
TIMEMIN = '1957-01-01'
TIMEMAX = '2014-12-31'

sam = sam.sel(time=slice(TIMEMIN, TIMEMAX))
RUTA = '/datos/osman/sam_trend_cmip_figures/'
#plot sam time serie
#fig1 = plt.figure(1, (10, 4), dpi=500)  #fig size in inches
#plt.plot(sam.time.values,sam.data, 'r', linewidth=1.5)
#plt.axhline(y=0,xmin=0,linestyle='--',linewidth=0.8,
#           color='k',alpha=0.5)
##ax.plot(2015,proyeccion_hgt_2015[i], 'b.', linewidth=2)
##ax.set_xlim((ti-1,2016))
#plt.ylim((-8, 8))
#plt.ylabel('SAM')
#plt.title('Southern Annular Mode 1957-2014',fontsize=12)
#fig1.savefig(RUTA + 'sam_time_series.png', dpi=500, bbox_inches='tight', papertype='A4', 
#             orientation='landscape')
#
#compute seasonal time series
seas = ['MAM', 'JJA', 'SON', 'DJF']
sam_seasonal = sam.resample(dim='time', freq='QS-MAR', how='mean')
#fig2 = plt.figure(1, (10, 14), dpi=500)
#x = np.arange(1957, 2015)
#
#for i in range(4):
#    ax = plt.subplot(5, 1, i+1)
#    line1 = ax.plot(sam_seasonal.time.values[sam_seasonal['time.month']==3 + i*3], sam_seasonal.sel(time=sam_seasonal['time.month']==3 + i*3), 'r', linewidth=1.5)
#    ax.axhline(y=0, xmin=0, linestyle='--', linewidth=0.8,
#               color='k', alpha=0.5)
#    # tidy up the figure
#    ax.set_ylim((-4.5, 4.5))
#    plt.title(seas[i], fontsize=10)
#plt.tight_layout()#pad=0.4, w_pad=0.5, h_pad=1.0)
#plt.subplots_adjust(top=0.93)
#plt.suptitle('Seasonal Southern Annular Mode', fontsize=12)
#fig2.savefig(RUTA + 'sam_seasonal.png',dpi=500, bbox_inches='tight', papertype='A4', 
#             orientation='landscape')
#
#compute trends
#annual
sam_annual = sam.groupby('time.year').mean(dim='time')
[slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(sam_annual.year.values)),
                                                   sam_annual.data)
print("Annual")
print(slope, p_val)
for i in range(4):
    print(seas[i])
    [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(sam_seasonal.time.values[sam_seasonal['time.month']==3 + i*3])),
                                                       sam_seasonal.sel(time=sam_seasonal['time.month']==3 + i*3))
    print(slope, p_val)


#half period
sam_early = sam_annual.sel(year=slice(1957, 1984))
sam_late = sam_annual.sel(year=slice(1985, 2014))
sam_seasonal_early = sam_seasonal.sel(time=slice('1957-01-01', '1984-12-31'))
sam_seasonal_late = sam_seasonal.sel(time=slice('1985-01-01', '2014-11-30'))

print('Early')
[slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(sam_early.data)),
                                                   sam_early.data)
print("Annual")
print(slope, p_val)
for i in range(4):
    print(seas[i])
    [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0,
                                                                 len(sam_seasonal_early.time.values[sam_seasonal_early['time.month']==3 + i*3])),
                                                       sam_seasonal_early.sel(time=sam_seasonal_early['time.month']==3 + i*3))
    print(slope, p_val)

print('Late')
[slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(sam_late.data)),
                                                   sam_late.data)
print("Annual")
print(slope, p_val)
for i in range(4):
    print(seas[i])
    [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0,
                                                                 len(sam_seasonal_late.time.values[sam_seasonal_late['time.month']==3 + i*3])),
                                                       sam_seasonal_late.sel(time=sam_seasonal_late['time.month']==3 + i*3))
    print(slope, p_val)

#10-year period
datei = ['1957-01-01', '1972-01-01', '1987-01-01', '2002-12-31']
yeari = [1957, 1972, 1987, 2002]
datef = ['1971-12-31', '1986-12-31', '2001-01-01', '2014-11-30']
yearf = [1971, 1986, 2001, 2014]
for i in np.arange(4):
    sam_decadal = sam_annual.sel(year=slice(yeari[i], yearf[i]))
    sam_seasonal_decadal = sam_seasonal.sel(time=slice(datei[i], datef[i]))
    print(datei[i], '-', datef[i])
    [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0, len(sam_decadal.data)),
                                                       sam_decadal.data)
    print("Annual")
    print(slope, p_val)
    for i in range(4):
        print(seas[i])
        [slope, interc, r_va, p_val, z] = stats.linregress(np.arange(0,
                                                                     len(sam_seasonal_decadal.time.values[sam_seasonal_decadal['time.month']==3 + i*3])),
                                                           sam_seasonal_decadal.sel(time=sam_seasonal_decadal['time.month']==3 + i*3))
        print(slope, p_val)


