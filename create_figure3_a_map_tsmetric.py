import numpy as np
import netCDF4 as ncdf
import Ngl
import datetime as dt
from default_ncl import *
from stat_ncl import *
from ncdf_read import select_region_TLL_files,select_region_TLLL_files
import matplotlib.pyplot as plt
from scipy import signal
from scipy.stats import t
import glob
dtrend=True

#dtrend=False
lon1=0.0; lon2=359.0; lat1=-60.0; lat2=60.0
lon1=0.0; lon2=360.0; lat1=-60.0; lat2=60.0
#lon1=120.0; lon2=240.0; lat1=-20.0; lat2=20.0
lev1=0.0; lev2=10.0
dt1=dt.datetime(1800,1,1,0,0,0)
dt2=dt.datetime(2019,12,31,0,0,0)
lonname="lon"; latname="lat"


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

fname_temp="MOAA_temp_merge_2004_2019.nc"
fname_salt="MOAA_salt_merge_2004_2019.nc"
fname_temp_grad_x="MOAA_grad_temp_x_merge_2004_2019.nc"
fname_temp_grad_y="MOAA_grad_temp_y_merge_2004_2019.nc"
fname_salt_grad_x="MOAA_grad_salt_x_merge_2004_2019.nc"
fname_salt_grad_y="MOAA_grad_salt_y_merge_2004_2019.nc"
fflag_sm="_moving9_res"
fname_temp_sm="MOAA_temp_merge_2004_2019"+fflag_sm+".nc"
fname_salt_sm="MOAA_salt_merge_2004_2019"+fflag_sm+".nc"
out_name="figure3_a_map_gradient_MOAA_2004_2019"
out_form="X11"
#out_form="eps"
varname_temp="temp";varname_salt="salt"

varname_temp_grad_x="temp_x"; varname_temp_grad_y="temp_y"
varname_salt_grad_x="salt_x"; varname_salt_grad_y="salt_y"

temp,miss_tt=select_region_TLLL_files([fname_temp],varname_temp,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
temp=np.nanmean(temp,axis=1)
salt,miss_tt=select_region_TLLL_files([fname_salt],varname_salt,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
salt=np.nanmean(salt,axis=1)

temp_sm,miss_tt=select_region_TLLL_files([fname_temp_sm],varname_temp,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
temp_sm=np.nanmean(temp_sm,axis=1)
salt_sm,miss_tt=select_region_TLLL_files([fname_salt_sm],varname_salt,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
salt_sm=np.nanmean(salt_sm,axis=1)
temp_sm_clm=clmmonTLL(temp_sm); temp_sm_anm=calcmonanomTLL(temp_sm,temp_sm_clm)
salt_sm_clm=clmmonTLL(salt_sm); salt_sm_anm=calcmonanomTLL(salt_sm,salt_sm_clm)


lon,lat,time,lev,temp_grad_x,miss_tt=select_region_TLLL_files([fname_temp_grad_x],varname_temp_grad_x,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)
temp_grad_x=np.nanmean(temp_grad_x,axis=1)

temp_grad_y=select_region_TLLL_files([fname_temp_grad_y],varname_temp_grad_y,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2)
temp_grad_y=np.nanmean(temp_grad_y,axis=1)


salt_grad_x=select_region_TLLL_files([fname_salt_grad_x],varname_salt_grad_x,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2)
salt_grad_x=np.nanmean(salt_grad_x,axis=1)
salt_grad_y=select_region_TLLL_files([fname_salt_grad_y],varname_salt_grad_y,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2)
salt_grad_y=np.nanmean(salt_grad_y,axis=1)

temp_grad_x=np.ma.masked_where(np.abs(temp_grad_x-miss_tt)<1.0e-3,temp_grad_x)
temp_grad_y=np.ma.masked_where(np.abs(temp_grad_y-miss_tt)<1.0e-3,temp_grad_y)
salt_grad_x=np.ma.masked_where(np.abs(salt_grad_x-miss_tt)<1.0e-3,salt_grad_x)
salt_grad_y=np.ma.masked_where(np.abs(salt_grad_y-miss_tt)<1.0e-3,salt_grad_y)


dtrend=True
dtrend=False
tsmooth=True
tsmooth=False
lonname="lon"; latname="lat"
fflag_sm=""
def dtrend_n(var,miss):
    var_mean=np.nanmean(var,axis=0)
    var_dt=signal.detrend(var,axis=0)+var_mean
    var_dt[np.abs(var)==miss]=np.nan
    return(var_dt)
wgt=np.ones(3)
wgt=wgt/np.sum(wgt)

if (tsmooth == True):
    temp_sm_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,temp_sm_anm)
    salt_sm_anm=np.apply_along_axis(lambda m:np.convolve(m,wgt,"same"),0,salt_sm_anm)


i=0
temp_mean=np.nanmean(temp,axis=0)
salt_mean=np.nanmean(salt,axis=0)
temp_xgrad_mean=np.nanmean(temp_grad_x,axis=0)
temp_xgrad_mean=np.ma.masked_invalid(temp_xgrad_mean)
temp_ygrad_mean=np.nanmean(temp_grad_y,axis=0)
temp_ygrad_mean=np.ma.masked_invalid(temp_ygrad_mean)
salt_xgrad_mean=np.nanmean(salt_grad_x,axis=0)
salt_xgrad_mean=np.ma.masked_invalid(salt_xgrad_mean)
salt_ygrad_mean=np.nanmean(salt_grad_y,axis=0)
salt_ygrad_mean=np.ma.masked_invalid(salt_ygrad_mean)

temp_hgrad_mean=np.nanmean(np.sqrt(temp_grad_x*temp_grad_x+temp_grad_y*temp_grad_y),axis=0)
temp_hgrad_mean=np.ma.masked_invalid(temp_hgrad_mean)
salt_hgrad_mean=np.nanmean(np.sqrt(salt_grad_x*salt_grad_x+salt_grad_y*salt_grad_y),axis=0)
salt_hgrad_mean=np.ma.masked_invalid(salt_hgrad_mean)


denom1=np.sqrt((temp_xgrad_mean*temp_xgrad_mean)+(temp_ygrad_mean*temp_ygrad_mean))
denom2=np.sqrt((salt_xgrad_mean*salt_xgrad_mean)+(salt_ygrad_mean*salt_ygrad_mean))
numer1=temp_xgrad_mean*salt_xgrad_mean+temp_ygrad_mean*salt_ygrad_mean



numer1=np.ma.masked_invalid(numer1)
denom1=np.ma.masked_invalid(denom1)
denom2=np.ma.masked_invalid(denom2)

angle_iso=numer1/(denom1*denom2)
angle_iso=angle_iso
angle_iso=np.ma.masked_invalid(angle_iso)

#print(np.max(angle_iso))
#print(np.min(angle_iso))

temp[np.abs(temp)>=50.0]=np.nan
temp_mean=np.nanmean(temp,axis=0)
temp_mean=np.ma.masked_invalid(temp_mean)
salt[np.abs(salt)>=50.0]=np.nan
salt_mean=np.nanmean(salt,axis=0)
salt_mean=np.ma.masked_invalid(salt_mean)
ndim=np.shape(temp)
ntime=ndim[0]

ts_cor_sm=get_cor_array(temp_sm_anm,salt_sm_anm)
r1=get_cor_array(temp_sm_anm[1:],temp_sm_anm[:-1])
r2=get_cor_array(salt_sm_anm[1:],salt_sm_anm[:-1])
dof=ntime*(1-r1*r2)/(1+r1*r2)
tval=t.ppf(q=0.975,df= dof.astype(int))
q=(tval**2)/(dof-2)
ts_dof=np.sqrt(q/(1+q))
#ts_cor_sm[np.abs(ts_cor_sm)<ts_dof]=np.nan
ts_cor_sm=np.ma.masked_invalid(ts_cor_sm)

mask_pos=np.copy(ts_dof)
mask_pos[:,:]=np.nan
mask_pos[ts_cor_sm>ts_dof]=1.0
mask_pos[ts_cor_sm<(-1*ts_dof)]=-1.0
mask_pos=np.ma.masked_invalid(mask_pos)



# angle_mean=np.nanmean(angle,axis=0)
# angle_mean=np.ma.masked_invalid(angle_mean)
# angle_x_mean=np.nanmean(angle_x,axis=0)
# angle_x_mean=np.ma.masked_invalid(angle_x_mean)
# angle_y_mean=np.nanmean(angle_y,axis=0)
# angle_y_mean=np.ma.masked_invalid(angle_y_mean)

wks=Ngl.open_wks(out_form,out_name)
res=Ngl.Resources()
res.cnFillOn=True
res.nglDraw=False
res.nglFrame=False
res.cnLineLabelsOn=False
res.mpFillOn=True
res.lbLabelBarOn=True
res.lbBoxEndCapStyle="TriangleBothEnds"
res.lbOrientation="Vertical"
res.lbOrientation="Horizontal"
res.mpGridLatSpacingF= 1000#lat_spacing
res.mpGridLonSpacingF= 1000#lat_spacing
res.tmXBLabelFontHeightF=0.025
res.tmYLLabelFontHeightF=0.025
res.lbLabelFontHeightF=0.02
res.tiMainFontHeightF=0.02
res.sfXArray=lon
res.sfYArray=lat
res.cnFillPalette="MPL_jet"
res.cnFillPalette="MPL_bwr"
res=setmap(res,lon1,lon2,lat1,lat2)
plot=[]

res.cnLinesOn=False
factor=1.0e5
#levels=[-1.5e-5,-1.0e-5,-0.5e-5,0.5e-5,1.0e-5,1.5e-5]
#res.cnLevels             = levels #     ; note 90, and not 100!
#plot[0]=Ngl.contour_map(wks,temp_hgrad_mean,res)
#plot.append(Ngl.contour(wks,temp_hgrad_mean,res))
res=setcint(res,-1.5e-5,1.5e-5,0.5e-5)
res=setcint(res,-1.5,1.5,0.5)
plot.append(Ngl.contour_map(wks,temp_ygrad_mean*factor,res))
levels=[-5.0e-6,-3.0e-6,-1.0e-6,1.0e-6,3.0e-6,5.0e-6]
levels=[-1.5e-6,-1.0e-6,-0.5e-6,0.5e-6,1.0e-6,1.5e-6]
#res.cnLineDashPatterns=dashes
#res.cnLevels = levels #     ; note 90, and not 100!
#plot[1]=Ngl.contour_map(wks,salt_hgrad_mean,res)
#plot.append(Ngl.contour(wks,salt_hgrad_mean,res))
res=setcint(res,-5.0e-6,5.0e-6,1.0e-6)
factor=1.0e6
res=setcint(res,-5.0,5.0,1.0)
plot.append(Ngl.contour_map(wks,salt_ygrad_mean*factor,res))

levels=[-1.0,-0.99,-0.9,-0.6,-0.3,0,0.3,0.6,0.9,0.99,1.0]
levels=[-1.0,-0.99,-0.9,-0.6,0.6,0.9,0.99,1.0]
#dashes=[0.0]*len(levels)
#dashes[0:5]=[5.0]*5
res.cnLevelSelectionMode = "ExplicitLevels"
res.cnFillOn=True
res.cnLinesOn=False
res.cnLevels             = levels #     ; note 90, and not 100!
res.cnInfoLabelOn=False
res.cnFillMode="RasterFill"
res.cnFillPalette="MPL_bwr"
res.tiMainString="SST"
plot.append(Ngl.contour_map(wks,angle_iso,res))

#plot.append(Ngl.contour_map(wks,mask_pos,res))

res.cnFillOn=False
res.cnLinesOn=True
res.cnLineColor="Grey20"
#res.cnLineColor="black"
res.cnLineThicknessF=1.0
levels=np.arange(0.0,30.0,2.0)
res.cnLineDashPatterns=[0]*len(levels)
res.cnLevels             = levels #     ; note 90, and not 100!
plot.append(Ngl.contour(wks,temp_mean,res))
levels=np.arange(31.0,38.0,0.4)
res.cnLineDashPatterns=[0]*len(levels)
res.cnLevels             = levels #     ; note 90, and not 100!
res.tiMainString="SSS"
plot.append(Ngl.contour(wks,salt_mean,res))
Ngl.overlay(plot[0],plot[3])
Ngl.overlay(plot[1],plot[4])

pattern=[6,0,17]
levels=[-0.3,0.3]
res.cnFillOpacityF=0.0
#res.cnFillOpacityF=0.1
res.cnMonoFillScale=True
res.cnFillScaleF=0.3
#res.cnFillDotSizeF=0.0015
#res.cnFillScaleF=0.5
res.cnFillDotSizeF=0.004
#res.cnFillDotSizeF=0.008
res.cnFillDotSizeF=0.02
res.cnFillScaleF=0.4
res.cnFillDotSizeF=0.002
res.cnFillColors=["Black","White","Black"]
res.cnFillColors=["White","White","White"]
res.cnFillColors=["Grey","White","Grey"]
res.cnFillColors=["cyan","White","cyan"]
#res.cnFillColors=["yellow2","White","lightcyan3"]
#res.cnFillColors=["brown","White","brown"]
#res.cnFillColors=["Black","White","Black"]
#res.cnFillColors=["Orange","White","Green"]
res.cnMonoFillPattern=False
res.cnFillPatterns=pattern
res.cnFillOn=True
res.cnLevelSelectionMode = "ExplicitLevels"
res.cnLevels             = levels
res.lbLabelBarOn=False
res.cnFillMode="AreaFill"
res.tiMainString="alpha"
plot.append(Ngl.contour(wks,mask_pos,res))
#Ngl.overlay(plot[2],plot[3])
Ngl.overlay(plot[2],plot[5])



pres=Ngl.Resources()
Ngl.panel(wks,plot[0:3],[1,3],pres)
#Ngl.panel(wks,plot[2:3],[1,1],pres)
#Ngl.panel(wks,plot[3:4],[1,1],pres)
Ngl.frame(wks)
Ngl.end()

