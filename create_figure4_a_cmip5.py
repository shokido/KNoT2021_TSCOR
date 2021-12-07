import numpy as np
import netCDF4 as ncdf
import Ngl
import datetime as dt
from default_ncl import *
from stat_ncl import *
from ncdf_read import select_region_TLL_files
import matplotlib.pyplot as plt
from scipy import signal

out_form="X11"
#out_form="eps"
out_name="figure4_a_cmip5"
#dtrend=False
lon1=0.0; lon2=359.0; lat1=-60.0; lat2=60.0; area_name="KE"
lon1=0.0; lon2=379.0
#lon1=140.0; lon2=160.0; lat1=30.0; lat2=50.0; area_name="KE"
#lon1=290.0; lon2=310.0; lat1=30.0; lat2=50.0; area_name="KE"

dt1=dt.datetime(1800,1,1,0,0,0)
dt2=dt.datetime(2019,12,31,0,0,0)
lonname="lon"; latname="lat"
varname_tt="cor_tt";varname_ts="cor_ts";varname_ss="cor_ss"
fflag_sm=""; fflag_sm2=""
fflag_sm="_move9_res";fflag_sm2="_moving9_res"
fflags=["","_moving9_res","_moving9"]
fname_model_list="namelist_model.txt"
model_names=[]
f_list=open(fname_model_list,"r")
l=f_list.readlines()
f_list.close()
nmodel_names=len(l)
print(nmodel_names)
for i in range(0,nmodel_names):
    model_name=l[i].strip()
    model_names.append(model_name)

i=0
iflag=0
nflags=len(fflags)
fname_in_model="out_ts_cor_historical_"+model_names[i]+"_360x180"+fflags[iflag]+".nc"    
lon_model,lat_model,time_model,cor_ts_model,miss_ts_model=select_region_TLL_files([fname_in_model],varname_ts,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)
nc_in_model=ncdf.Dataset(fname_in_model,"r")
lags=nc_in_model.variables["lags"][:]
nc_in_model.close()
cor_ts_models=np.zeros((nmodel_names,nflags,len(time_model),len(lat_model),len(lon_model)))
for i in range(0,nmodel_names):
#    print(i)
    for iflag in range(0,nflags):
        fname_in_model="out_ts_cor_historical_"+model_names[i]+"_360x180"+fflags[iflag]+".nc"
        cor_ts_model,miss_ts_model=select_region_TLL_files([fname_in_model],varname_ts,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True)
        cor_ts_model=np.nan_to_num(cor_ts_model,nan=miss_ts_model)
        cor_ts_model[cor_ts_model==miss_ts_model]=np.nan
    
        #cor_ts_model=np.ma.masked_invalid(cor_ts_model)
        cor_ts_models[i,iflag,:,:,:]=cor_ts_model[:,:,:]

cor_ts_models_mean=np.nanmean(cor_ts_models,axis=0)
cor_ts_models_mean=np.ma.masked_invalid(cor_ts_models_mean)
cor_ts_models=np.ma.masked_invalid(cor_ts_models)

# ogcm_names=["OFES1","OFES2"]
# dir_ogcm="/latent/kido/Salinity/TS_COR/CALC/"
# for i in range(0,nogcm_names):
# #    print(i)
#     fname_in_ogcm=dir_ogcm+"out_ts_cor_"+ogcm_names[i]+"_360x180"+fflag_sm2+".nc"



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
res.lbOrientation="Vertical"
res.pmLabelBarHeightF=0.6
res.pmLabelBarWidthF=0.1
res.lbLabelFontHeightF=0.02
res.mpGridLatSpacingF= 1000#lat_spacing
res.mpGridLonSpacingF= 1000#lat_spacing
res.tmXBLabelFontHeightF=0.025
res.tmYLLabelFontHeightF=0.025
res.sfXArray=lon_model
res.sfYArray=lat_model
res.cnFillPalette="MPL_bwr"
res=setmap(res,lon1,lon2,lat1,lat2)
res=setcint(res,-1.0,1.0,0.1)
plot=[]

lags_out=[0]
nlag=len(lags_out)
#for i in range(0,nlag):
index_out=np.where(lags==0)[0][0]

res.tiMainString="CMIP5 (Lags="+str(int(lags[index_out]))+")"
plot.append(Ngl.contour_map(wks,cor_ts_models_mean[0,index_out,:,:],res))
plot.append(Ngl.contour_map(wks,cor_ts_models_mean[1,index_out,:,:],res))
plot.append(Ngl.contour_map(wks,cor_ts_models_mean[2,index_out,:,:],res))

pres=Ngl.Resources()
Ngl.panel(wks,plot[0:3],[1,3],pres)
#Ngl.frame(wks)
Ngl.end()

