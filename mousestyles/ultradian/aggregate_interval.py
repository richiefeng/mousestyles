import numpy as np
from matplotlib import pyplot as plt
from matplotlib.externals import six


import pandas as pd

import numpy.random

import datetime
from scipy.signal import lombscargle
import os
import sys

os.chdir('../../')
print(os.getcwd())

import mousestyles.data as data

def aggegate_interval(strain, mouse, feature, bin_width):
    """
    data loaded from data.load_intervals(feature)

    Parameters
    ---------------
    feature: {"AS", "Food", "IS", "M_AS", "M_IS", "Water", "Distance", "AS_Intensity", "AS_prob"}
    bin_width: number of minutes of time interval for data aggregation

    Returns
    ----------
    ts: pandas.tseries
        a pandas time series of length 12(day)*24(hour)*60(minute)/n
    """
    
    path = 'mousestyles/data/txy_coords/recordingStartTimeEndTime'
    files = []
    le=24*60/bin_width
    tim=list()
    a = datetime.datetime(100,1,1,0,0,0)
    for i in range(int(le)):
        b = a + datetime.timedelta(0,i*bin_width*60)
    tim.append(b.time())
    intervals=data.load_intervals(feature)
    inn=intervals.loc[intervals['strain'] == strain]
    innn=inn.loc[intervals['mouse'] == mouse]
    a='recordingStartTimeEndTime_strain'
    a=a+str(strain)
    a=a+'_mouse'
    a=a+str(mouse)
    for i in os.listdir(path):
        if os.path.isfile(os.path.join(path,i)) and a in i:
            files.append(i)
    tt=list()
    for i in range(len(files)):
        pa=path+'/'+str(files[i])
        b=np.load(pa)
        start=b[0]
        end=b[0]+3600*24
        ini=innn.loc[innn['day'] == i]
        for j in range(int(le)):
            aa=ini.loc[((ini['start'] >start+j*bin_width*60) & (ini['start'] <start+(j+1)*bin_width*60))|((ini['stop'] >start+j*bin_width*60) & (ini['stop'] <start+(j+1)*bin_width*60))]
            g=0
            gg=list(ini.iloc[:]['stop']-ini.iloc[:]['start'])
            if len(aa)==0:
                g=0
            if len(aa)!=0:
                if aa.iloc[0]['start']>=start+j*bin_width*60 and aa.iloc[len(aa)-1]['stop']<=start+(j+1)*bin_width*60:
                    for k in range(len(aa)):
                        g=g+gg[k]
                if aa.iloc[0]['start']<start+j*bin_width*60 and aa.iloc[len(aa)-1]['stop']<=start+(j+1)*bin_width*60:
                    for k in range(len(aa)-1):
                        g=g+gg[k+1]
                    g=g+aa.iloc[0]['stop']-start-j*bin_width*60
                if aa.iloc[0]['start']>=start+j*bin_width*60 and aa.iloc[len(aa)-1]['stop']>start+(j+1)*bin_width*60:
                    for k in range(len(aa)-1):
                        g=g+gg[k]
                    g=g+start+(j+1)*bin_width*60-aa.iloc[len(aa)-1]['start']
                if aa.iloc[0]['start']<start+j*bin_width*60 and aa.iloc[len(aa)-1]['stop']>start+(j+1)*bin_width*60:
                    for k in range(len(aa)-2):
                        g=g+gg[k+1]
                    g=g+aa.iloc[0]['stop']-start-j*bin_width*60
                    g=g+start+(j+1)*bin_width*60-aa.iloc[len(aa)-1]['start']
            tt.append(g)
    tre=str(bin_width)+'min'
    ts=pd.DataFrame(tt, index=pd.date_range('2014-01-01',periods=len(tt),freq=tre))
    ts.columns=[feature]
    return(ts)


#tsg=aggegate_interval(0, 0, 'F', 20)
#print (tsg)