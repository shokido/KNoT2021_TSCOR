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
dtrend=True
lonname="lon"; latname="lat"
fflag_sm="_moving9"
t_thres=0.975
fname_nhffwf="out_nhffwf_lag_JOFURO3_V1.1_2004_2017"+fflag_sm+".nc"
fname_tttend="out_ttendnhf_lag_MOAA_JOFURO3_V1.1_2004_2017"+fflag_sm+".nc"
fname_sstend="out_stendfwf_lag_MOAA_JOFURO3_V1.1_2004_2017"+fflag_sm+".nc"

fname_temp_sm="MOAA_temp_merge_2004_2019"+fflag_sm+".nc"
fname_salt_sm="MOAA_salt_merge_2004_2019"+fflag_sm+".nc"
varname_temp="temp";varname_salt="salt"

out_name="figure3_b_cor_nhffwf_JOFURO3_2004_2017"+fflag_sm
out_form="X11"
#out_form="eps"
#dtrend=False
lon1=0.0; lon2=360.0; lat1=-60.0; lat2=60.0; area_name="KE"
lev1=0.0; lev2=10.0

varname_nhffwf="cor_nhffwf"
varname_r1_nhffwf="cor_nhfnhf"
varname_r2_nhffwf="cor_fwffwf"

varname_tttend="cor_ttendnhf"
varname_sstend="cor_stendfwf"
varname_r1_tttend="cor_nhfnhf"; varname_r2_tttend="cor_ttendttend"
varname_r1_sstend="cor_fwffwf"; varname_r2_sstend="cor_stendstend"


dt1=dt.datetime(1800,1,1,0,0,0)
dt2=dt.datetime(2019,12,31,0,0,0)

nc_nhffwf=ncdf.Dataset(fname_nhffwf,"r")
lags=nc_nhffwf.variables["lags"][:]
nc_nhffwf.close()

lon,lat,time,cor_nhffwf,miss_nhffwf=select_region_TLL_files([fname_nhffwf],varname_nhffwf,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)
r1_nhffwf=select_region_TLL_files([fname_nhffwf],varname_r1_nhffwf,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
r2_nhffwf=select_region_TLL_files([fname_nhffwf],varname_r2_nhffwf,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
lags_out=[0]
nlag=len(lags_out)
index_out=np.where(lags==1)[0][0]
r1_nhffwf=r1_nhffwf[index_out,:,:]
r2_nhffwf=r2_nhffwf[index_out,:,:]
ntime=12*10
dof=ntime*(1-r1_nhffwf*r2_nhffwf)/(1+r1_nhffwf*r2_nhffwf)
tval=t.ppf(q=t_thres,df= dof.astype(int))
q=(tval**2)/(dof-2)
r_thres=np.sqrt(q/(1+q))
#ntime*(1-r1_nhffwf*r2_nhffwf)/(1+r1_nhffwf*r2_nhffwf)

cor_nhffwf=cor_nhffwf*(-1) ## Positive NHF means loss of oceanic heat, so sign needs to be flipped
cor_nhffwf[np.abs(cor_nhffwf)<=r_thres]=np.nan
cor_nhffwf=np.ma.masked_invalid(cor_nhffwf)


# Ttend
lon,lat,time,cor_tttend,miss_tttend=select_region_TLL_files([fname_tttend],varname_tttend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)
r1_tttend=select_region_TLL_files([fname_tttend],varname_r1_tttend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
r2_tttend=select_region_TLL_files([fname_tttend],varname_r2_tttend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
lags_out=[0]
nlag=len(lags_out)
index_out=np.where(lags==1)[0][0]
r1_tttend=r1_tttend[index_out,:,:]
r2_tttend=r2_tttend[index_out,:,:]
ntime=12*10
dof=ntime*(1-r1_tttend*r2_tttend)/(1+r1_tttend*r2_tttend)
tval=t.ppf(q=t_thres,df= dof.astype(int))
q=(tval**2)/(dof-2)
r_thres=np.sqrt(q/(1+q))

cor_tttend=cor_tttend*(-1) ## Positive NHF means loss of oceanic heat, so sign needs to be flipped
cor_tttend[np.abs(cor_tttend)<=r_thres]=np.nan
cor_tttend=np.ma.masked_invalid(cor_tttend)

lon,lat,time,cor_sstend,miss_sstend=select_region_TLL_files([fname_sstend],varname_sstend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)
r1_sstend=select_region_TLL_files([fname_sstend],varname_r1_sstend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
r2_sstend=select_region_TLL_files([fname_sstend],varname_r2_sstend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
lags_out=[0]
cor_sstend=cor_sstend*(-1.0)
nlag=len(lags_out)
index_out=np.where(lags==1)[0][0]
r1_sstend=r1_sstend[index_out,:,:]
r2_sstend=r2_sstend[index_out,:,:]
ntime=12*10
dof=ntime*(1-r1_sstend*r2_sstend)/(1+r1_sstend*r2_sstend)
tval=t.ppf(q=t_thres,df= dof.astype(int))
q=(tval**2)/(dof-2)
r_thres=np.sqrt(q/(1+q))

cor_sstend=cor_sstend*(-1) ## Positive NHF means loss of oceanic heat, so sign needs to be flipped
cor_sstend[np.abs(cor_sstend)<=r_thres]=np.nan
cor_sstend=np.ma.masked_invalid(cor_sstend)

# T-S correlation
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
temp_sm,miss_tt=select_region_TLLL_files([fname_temp_sm],varname_temp,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
temp_sm=np.nanmean(temp_sm,axis=1)
salt_sm,miss_tt=select_region_TLLL_files([fname_salt_sm],varname_salt,dt1,dt2,lev1,lev2,lat1,lat2,lon1,lon2,return_missing=True)
salt_sm=np.nanmean(salt_sm,axis=1)
temp_sm_clm=clmmonTLL(temp_sm); temp_sm_anm=calcmonanomTLL(temp_sm,temp_sm_clm)
salt_sm_clm=clmmonTLL(salt_sm); salt_sm_anm=calcmonanomTLL(salt_sm,salt_sm_clm)

ndim=np.shape(temp_sm_anm)
ntime=ndim[0]
ts_cor_sm=get_cor_array(temp_sm_anm,salt_sm_anm)
r1=get_cor_array(temp_sm_anm[1:],temp_sm_anm[:-1])
r2=get_cor_array(salt_sm_anm[1:],salt_sm_anm[:-1])
dof=ntime*(1-r1*r2)/(1+r1*r2)
tval=t.ppf(q=t_thres,df= dof.astype(int))
q=(tval**2)/(dof-2)
ts_dof=np.sqrt(q/(1+q))
ts_cor_sm=np.ma.masked_invalid(ts_cor_sm)
mask_pos=np.copy(ts_dof)
mask_pos[:,:]=np.nan
mask_pos[ts_cor_sm>ts_dof]=1.0
mask_pos[ts_cor_sm<(-1*ts_dof)]=-1.0
mask_pos=np.ma.masked_invalid(mask_pos)

wks=Ngl.open_wks(out_form,out_name)
res=Ngl.Resources()
res.cnFillOn=True
res.nglDraw=False
res.nglFrame=False
res.cnLineLabelsOn=False
res.cnLinesOn=False
res.mpFillOn=True
res.lbLabelBarOn=True
res.lbBoxEndCapStyle="TriangleBothEnds"
res.lbOrientation="Horizontal"
res.lbLabelFontHeightF=0.02
res.mpGridLatSpacingF= 1000#lat_spacing
res.mpGridLonSpacingF= 1000#lat_spacing
res.tmXBLabelFontHeightF=0.025
res.tmYLLabelFontHeightF=0.025
res.sfXArray=lon
res.sfYArray=lat
res.cnFillPalette="MPL_bwr"
res=setmap(res,lon1,lon2,lat1,lat2)
plot=[]

lags_out=[0]
nlag=len(lags_out)
res=setcint(res,-1.0,1.0,0.1)

index_out=np.where(lags==lags_out[0])[0][0]
res.tiMainString="Cor(NHF,dSST/dt)"
plot.append(Ngl.contour_map(wks,cor_tttend[index_out,:,:],res))
res.tiMainString="Cor(FWF,dSSS/dt)"
plot.append(Ngl.contour_map(wks,cor_sstend[index_out,:,:],res))
res.tiMainString="Cor(NHF,FWF)"
plot.append(Ngl.contour_map(wks,cor_nhffwf[index_out,:,:],res))
#plot.append(Ngl.contour_map(wks,ts_cor_sm,res))
#plot.append(Ngl.contour_map(wks,ts_dof,res))

# Stripping
pattern=[6,0,17]
#pattern=[17,0,17]
#pattern=[17,0,0]
#pattern=[0,0,0]
levels=[-0.3,0.3]
res.cnFillOpacityF=0.0
#res.cnFillOpacityF=0.1
res.cnMonoFillScale=True
res.cnMissingValFillColor=-1
res.cnFillScaleF=0.3
res.cnFillDotSizeF=0.002
res.cnFillColors=["Black","White","Black"]
res.cnFillColors=["White","White","White"]
res.cnFillColors=["Grey80","White","Grey80"]
#res.cnFillColors=["purple","White","purple"]
res.cnFillColors=["green","White","w"]
res.cnFillColors=["cyan","White","cyan"]
res.cnMonoFillPattern=False
res.cnFillPatterns=pattern
res.cnFillOn=True
res.cnLevelSelectionMode = "ExplicitLevels"
res.cnLevels             = levels
res.lbLabelBarOn=False
res.cnFillMode="AreaFill"
#res.cnLinesOn=True
plot.append(Ngl.contour(wks,mask_pos,res))
plot.append(Ngl.contour(wks,mask_pos,res))
plot.append(Ngl.contour(wks,mask_pos,res))
Ngl.overlay(plot[0],plot[3])
Ngl.overlay(plot[1],plot[4])
Ngl.overlay(plot[2],plot[5])



pres=Ngl.Resources()
Ngl.panel(wks,plot[0:3],[3,1],pres)
Ngl.frame(wks)
Ngl.end()

