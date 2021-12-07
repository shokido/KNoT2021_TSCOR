import numpy as np
import netCDF4 as ncdf
import Ngl
import datetime as dt
from default_ncl import *
from stat_ncl import *
from ncdf_read import select_region_TLLL_files
import matplotlib.pyplot as plt
from scipy.stats import t
from scipy import signal

dtrend=True
dtrend=False
tsmooth=True

tsmooth=False

lat1=-70.0; lat2=70.0
lon1=0.0; lon2=379.0
lon1=0.0; lon2=359.0
#lat1=-50.0; lat2=50.0
lev1=0.0; lev2=10.0
dt1=dt.datetime(2004,1,1,0,0,0)
dt2=dt.datetime(2019,12,31,0,0,0)
#dt2=dt.datetime(2013,12,31,0,0,0)
out_form="X11"
#out_form="eps"
data_flag="RG";
#data_flag="MOAA"
lonname="lon"; latname="lat"
#out_name="figure_variance_sscale_"+data_flag+"_moving15x15"
fflags=["","_moving9_res","_moving9"]
out_name="figure2_std_sscale_"+data_flag+"_moving9x9_"+data_flag
#fflags=["","_move15"]
nflags=len(fflags)
def dtrend_n(var,miss):
    var_mean=np.nanmean(var,axis=0)
    var_dt=signal.detrend(var,axis=0)+var_mean
    var_dt[np.abs(var)==miss]=np.nan
    return(var_dt)

def get_cor_array(x1,x2):
    ndim=np.shape(x1)
    x1_tmp=np.reshape(x1,[ndim[0],ndim[1]*ndim[2]])
    x2_tmp=np.reshape(x2,[ndim[0],ndim[1]*ndim[2]])
    x1_tmp=x1_tmp-np.mean(x1_tmp,axis=0)
    x2_tmp=x2_tmp-np.mean(x2_tmp,axis=0)
    x1x1=np.mean(x1_tmp*x1_tmp,axis=0)
    x2x2=np.mean(x2_tmp*x2_tmp,axis=0)
    x1x2=np.mean(x1_tmp*x2_tmp,axis=0)
    cor=x1x2/(np.sqrt(x1x1)*np.sqrt(x2x2))
    cor_2=np.reshape(cor,[ndim[1],ndim[2]])
    return(cor_2)

wgt=np.ones(3)
wgt=wgt/np.sum(wgt)
i=0
fname_temp="MOAA_temp_merge_2004_2019"+fflags[i]+".nc"
fname_salt="MOAA_salt_merge_2004_2019"+fflags[i]+".nc"
varname_temp="temp"
varname_salt="salt"
lon,lat,lev,time,temp,miss_temp=select_region_TLLL_files([fname_temp],varname_temp,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)

# Get correlation
ntime=len(time)
nlat=len(lat); nlon=len(lon)
temp_std=np.zeros((nflags,nlat,nlon))
salt_std=np.zeros((nflags,nlat,nlon))
ts_cor=np.zeros((nflags,nlat,nlon))
ts_dof=np.zeros((nflags,nlat,nlon))
for i in range(0,nflags):
    fname_temp="MOAA_temp_merge_2004_2019"+fflags[i]+".nc"
    fname_salt="MOAA_salt_merge_2004_2019"+fflags[i]+".nc"
    temp,miss_temp=select_region_TLLL_files([fname_temp],varname_temp,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
    salt,miss_salt=select_region_TLLL_files([fname_salt],varname_salt,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
    temp[np.abs(temp-miss_temp)<1.0e-3]=np.nan
    temp=np.nanmean(temp,axis=1)
    salt,miss_salt=select_region_TLLL_files([fname_salt],varname_salt,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
    salt[np.abs(salt-miss_salt)<1.0e-3]=np.nan
    salt=np.nanmean(salt,axis=1)

    if (dtrend == True):
        temp=dtrend_n(temp,miss_temp)
        salt=dtrend_n(salt,miss_salt)

    temp_clm=clmmonTLL(temp); temp_anm=calcmonanomTLL(temp,temp_clm)
    salt_clm=clmmonTLL(salt); salt_anm=calcmonanomTLL(salt,salt_clm)

    if (tsmooth == True):
        temp_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_anm)
        salt_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,salt_anm)
    temp_std[i,:,:]=np.var(temp_anm,axis=0)
    salt_std[i,:,:]=np.var(salt_anm,axis=0)
    temp_std[i,:,:]=np.std(temp_anm,axis=0)
    salt_std[i,:,:]=np.std(salt_anm,axis=0)
    ts_cor[i,:,:]=get_cor_array(temp_anm,salt_anm)
    r1=get_cor_array(temp_anm[1:],temp_anm[:-1])
    r2=get_cor_array(salt_anm[1:],salt_anm[:-1])
    dof=ntime*(1-r1*r2)/(1+r1*r2)
    tval=t.ppf(q=0.975,df= dof.astype(int))
    q=(tval**2)/(dof-2)
    ts_dof[i,:,:]=np.sqrt(q/(1+q))
    #ntime*(1-r1*r2)/(1+r1*r2)
#    dof=
#    tval=dof
    
temp_std=np.ma.masked_invalid(temp_std)
salt_std=np.ma.masked_invalid(salt_std)

ts_cor=np.ma.masked_invalid(ts_cor)
ts_dof=np.ma.masked_invalid(ts_dof)
ts_cor[np.abs(ts_cor)<ts_dof]=np.nan
ts_cor=np.ma.masked_invalid(ts_cor)



wks=Ngl.open_wks(out_form,out_name)
res=Ngl.Resources()
res.cnFillOn=True
res.nglDraw=False
res.nglFrame=False
res.cnLineLabelsOn=False
res.mpFillOn=True
res.lbBoxEndCapStyle="TriangleBothEnds"
res.lbOrientation="Horizontal"
res.lbOrientation="Vertical"
res.lbLabelFontHeightF=0.02
res.mpGridLatSpacingF= 1000#lat_spacing
res.mpGridLonSpacingF= 1000#lat_spacing
res.tmXBLabelFontHeightF=0.025
res.tmYLLabelFontHeightF=0.025
res.pmLabelBarHeightF=0.6
res.pmLabelBarWidthF=0.1
res.sfXArray=lon
res.sfYArray=lat
res=setmap(res,lon1,lon2,lat1,lat2)
res.cnLinesOn=False
plot=[]#*9
titles=[" (Total)"," (Small scale)"," (Large scale)"]
for i in range(0,nflags):
    res.tiMainString="SST variance"+titles[i]
    res=setcint(res,0.0,1.5,0.1)
    res.cnFillPalette="MPL_Oranges"
    res.cnFillPalette="MPL_jet"
    plot.append(Ngl.contour_map(wks,temp_std[i,:,:],res))
    res.tiMainString="SSS variance"+titles[i]
 #   res=setcint(res,0.0,0.1,0.01)
    res=setcint(res,0.0,0.3,0.03)
    plot.append(Ngl.contour_map(wks,salt_std[i,:,:],res))
    res.cnFillPalette="MPL_bwr"
    #res.cnFillPalette="testcmap"
    res.tiMainString="Cor(SST,SSS) "+titles[i]
    res=setcint(res,-1.0,1.0,0.1)
#    res=setcint(res,0.0,ntime,10)
    plot.append(Ngl.contour_map(wks,ts_cor[i,:,:],res))
#    plot.append(Ngl.contour_map(wks,ts_dof[i,:,:],res))

# res.cnFillPalette="MPL_jet"    
# res=setcint(res,10.0,100.0,10)
# res.tiMainString="Percentage of small scale (SST)"
# plot.append(Ngl.contour_map(wks,100*temp_std[nflags-1,:,:]/temp_std[0,:,:],res))
# res.tiMainString="Percentage of small scale (SSS)"
# plot.append(Ngl.contour_map(wks,100*salt_std[nflags-1,:,:]/salt_std[0,:,:],res))

arand=list([0,3,6,1,4,7,2,5,8])
plot2=[]
plot2.append(plot[0])
plot2.append(plot[3])
plot2.append(plot[6])
plot2.append(plot[1])
plot2.append(plot[4])
plot2.append(plot[7])
plot2.append(plot[2])
plot2.append(plot[5])
plot2.append(plot[8])
pres=Ngl.Resources()
Ngl.panel(wks,plot2[6:9],[1,3],pres)
Ngl.frame(wks)
Ngl.end()

