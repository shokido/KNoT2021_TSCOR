import numpy as np
import netCDF4 as ncdf
import Ngl
import datetime as dt
from default_ncl import *
from stat_ncl import *
from ncdf_read import select_region_TLL_files
import matplotlib.pyplot as plt
from scipy import signal
from scipy.stats import t
out_form="X11"
#out_form="eps"
out_name="figure4_b_scatter_multimodel.eps"

#dtrend=False
lon1=0.0; lon2=359.0; lat1=-60.0; lat2=60.0; area_name="KE"
#lon1=40.0; lon2=260.0; lat1=-60.0; lat2=50.0; area_name="KE"
#lon1=300.0; lon2=320.0; lat1=35.0; lat2=42.0; area_name="KE"

dt1=dt.datetime(1800,1,1,0,0,0)
dt2=dt.datetime(2019,12,31,0,0,0)
lonname="lon"; latname="lat"
varname_tt="cor_tt";varname_ts="cor_ts";varname_ss="cor_ss"
fflags_1=["","_moving9_res","_moving9"]
fflags_2=["","_moving9_res","_moving9"]
captions=["Total","Small scale","Large scale"]
nflags=len(fflags_1)
fflag_sm=""; fflag_sm2=""
fnames_in_obs=["out_ts_lag_MOAA_2004_2019"+fflags_1[0]+".nc","out_ts_lag_RG_2004_2019"+fflags_1[0]+".nc"]
obs_names=["MOAA","RG"]

lags_out=[0]
nlag=len(lags_out)

nobs=len(obs_names)
lon_obs,lat_obs,time_obs,cor_ts_obs,miss_ts_obs=select_region_TLL_files([fnames_in_obs[1]],varname_ts,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)
cor_ts_obss=np.zeros((nobs,nflags,1,len(lat_obs),len(lon_obs)))
for i in range(0,nobs):
    for iflag in range(0,nflags):
        fname_in_obs="out_ts_lag_"+obs_names[i]+"_2004_2019"+fflags_1[iflag]+".nc"
        nc_in=ncdf.Dataset(fname_in_obs,"r")
        lags=nc_in.variables["lags"][:]
        index_out=np.where(lags==0)[0][0]
        nc_in.close()

        cor_ts_obs,miss_ts_obs=select_region_TLL_files([fname_in_obs],varname_ts,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True)
        cor_ts_obs=np.nan_to_num(cor_ts_obs,nan=miss_ts_obs)
        cor_ts_obs[cor_ts_obs==miss_ts_obs]=np.nan
        cor_ts_obss[i,iflag,0,:,:]=cor_ts_obs[index_out,:,:]
    
#cor_ts_obss=np.nanmean(cor_ts_obss,axis=3)
cor_ts_obss=np.ma.masked_invalid(cor_ts_obss)
cor_ts_obss_mean=np.nanmean(cor_ts_obss,axis=0)

fname_model_list="namelist_model.txt"
model_names=[]
f_list=open(fname_model_list,"r")
l=f_list.readlines()
f_list.close()
nmodel_names=len(l)

for i in range(0,nmodel_names):
    model_name=l[i].strip()
    model_names.append(model_name)
cor_ts_models=np.zeros((nmodel_names,nflags,1,len(lat_obs),len(lon_obs)))
for i in range(0,nmodel_names):
    for iflag in range(0,nflags):
        fname_in_model="out_ts_cor_historical_"+model_names[i]+"_360x180"+fflags_2[iflag]+".nc"
        nc_in=ncdf.Dataset(fname_in_model,"r")
        lags=nc_in.variables["lags"][:]
        index_out=np.where(lags==0)[0][0]
        nc_in.close()
        cor_ts_model,miss_ts_model=select_region_TLL_files([fname_in_model],varname_ts,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True)
        cor_ts_model=np.nan_to_num(cor_ts_model,nan=miss_ts_model)
        cor_ts_model[cor_ts_model==miss_ts_model]=np.nan

        cor_ts_model=np.ma.masked_invalid(cor_ts_model)
        cor_ts_models[i,iflag,0,:,:]=cor_ts_model[index_out,:,:]

#cor_ts_models=np.nanmean(cor_ts_models,axis=3)
cor_ts_models_std=np.nanstd(cor_ts_models,axis=0)
cor_ts_models_mean=np.nanmean(cor_ts_models,axis=0)
cor_ts_models_mean=np.ma.masked_invalid(cor_ts_models_mean)
cor_ts_models=np.ma.masked_invalid(cor_ts_models)

ogcm_names=["OFES1","OFES2"]
dir_ogcm=""
nogcm_names=len(ogcm_names)
cor_ts_ogcms=np.zeros((nogcm_names,nflags,1,len(lat_obs),len(lon_obs)))
for i in range(0,nogcm_names):
    for iflag in range(0,nflags):
        fname_in_ogcm=dir_ogcm+"out_ts_cor_"+ogcm_names[i]+"_360x180"+fflags_2[iflag]+".nc"
        nc_in=ncdf.Dataset(fname_in_ogcm,"r")
        lags=nc_in.variables["lags"][:]
        index_out=np.where(lags==0)[0][0]
        nc_in.close()
        cor_ts_ogcm,miss_ts_ogcm=select_region_TLL_files([fname_in_ogcm],varname_ts,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True)
        cor_ts_ogcm=np.nan_to_num(cor_ts_ogcm,nan=miss_ts_ogcm)
        cor_ts_ogcm[cor_ts_ogcm==miss_ts_ogcm]=np.nan

        #cor_ts_ogcm=np.ma.masked_invalid(cor_ts_ogcm)
        cor_ts_ogcms[i,iflag,0,:,:]=cor_ts_ogcm[index_out,:,:]

#cor_ts_ogcms=np.nanmean(cor_ts_ogcms,axis=3)
cor_ts_ogcms_std=np.nanstd(cor_ts_ogcms,axis=0)
cor_ts_ogcms_mean=np.nanmean(cor_ts_ogcms,axis=0)
cor_ts_ogcms_mean=np.ma.masked_invalid(cor_ts_ogcms_mean)
cor_ts_ogcms=np.ma.masked_invalid(cor_ts_ogcms)

#cor_ts_obss=cor_ts_obss[:,:,index_out,:]
#cor_ts_ogcms=cor_ts_ogcms[:,:,index_out,:]
#cor_ts_models=cor_ts_models[:,:,index_out,:]
print(np.shape(cor_ts_obss))
print(np.shape(cor_ts_ogcms))
cor_ts_obss_mean=np.nanmean(cor_ts_obss,axis=0)
cor_ts_ogcms_mean=np.nanmean(cor_ts_ogcms,axis=0)
cor_ts_models_mean=np.nanmean(cor_ts_models,axis=0)

slope_ogcms=np.zeros((nogcm_names,nflags))
slope_models=np.zeros((nmodel_names,nflags))
fig = plt.figure(figsize=(15.0,5.0))
#fig = plt.figure(figsize=(9.0,3.0))
thres=0.0

for iflag in range(0,nflags):
    plt.subplot(1,nflags,iflag+1)
    a=cor_ts_obss[0,iflag,:]

    a=np.ravel(a)

    a_ori=np.copy(a)
    a_ori[np.abs(a)<thres]=np.nan
    print(np.shape(a_ori))
    for i in range(0,nobs):
        b=cor_ts_obss[i,iflag,:]
        b=np.ravel(b)
        b[np.abs(a_ori)<thres]=np.nan
        a=a_ori[~np.isnan(a_ori)&(~np.isnan(b))]
        b=b[(~np.isnan(a_ori))&(~np.isnan(b))]

        z = np.polyfit(a, b, 1)
        print(z)

    print("OGCM")

    for i in range(0,nogcm_names):
        b=cor_ts_ogcms[i,iflag,:]
        b=np.ravel(b)
        b[np.abs(a_ori)<thres]=np.nan
        a=a_ori[~np.isnan(a_ori)&(~np.isnan(b))]
        b=b[(~np.isnan(a_ori))&(~np.isnan(b))]
        #a=np.abs(a);b=np.abs(b)

        z = np.polyfit(a, b, 1)
        print(z[0])
        slope_ogcms[i,iflag]=z[0]
        #     #plt.scatter(a,b,color="red",s=5,alpha=0.02)
    b=cor_ts_ogcms_mean[iflag,:,:]
    b=np.ravel(b)
    a=a_ori[~np.isnan(a_ori)&(~np.isnan(b))]
    b=b[(~np.isnan(a_ori))&(~np.isnan(b))]
    #a=np.abs(a);b=np.abs(b)
    z = np.polyfit(a, b, 1)
    print(z[0])
    plt.scatter(a,b,color="skyblue",marker="x",s=2,alpha=0.2)
    # print(np.mean(slope_ogcms))
    print("CGCM")
    for i in range(0,nmodel_names):
        b=cor_ts_models[i,iflag,:]
        b=np.ravel(b)
        a=a_ori[~np.isnan(a_ori)&(~np.isnan(b))]
        b=b[(~np.isnan(a_ori))&(~np.isnan(b))]
        #a=np.abs(a);b=np.abs(b)

        z = np.polyfit(a, b, 1)
        slope_models[i,iflag]=z[0]
        print(z[0])
        #plt.scatter(a,b,color="lime",marker="x",s=2,alpha=0.02)
        # print(np.mean(slope_models))

    b=cor_ts_models_mean[iflag,:]
    b=np.ravel(b)
    a=a_ori[~np.isnan(a_ori)&(~np.isnan(b))]
    b=b[(~np.isnan(a_ori))&(~np.isnan(b))]
    #a=np.abs(a);b=np.abs(b)
    z = np.polyfit(a, b, 1)
    print(np.shape(b))
    print(z[0])
    plt.scatter(a,b,color="lime",marker="x",s=2,alpha=0.1)


    x=np.linspace(-1,1,11)
    y=np.zeros(11)
    plt.plot(x,y,color="grey",linestyle='dashed')
    x=np.zeros(11)
    y=np.linspace(-1,1,11)
    plt.plot(x,y,color="grey",linestyle='dashed')
    x=np.linspace(-1,1,11)
    y=x
    plt.plot(x,y,color="grey",linewidth=2.0)
    y_model=x*np.mean(slope_models[:,iflag])
    plt.plot(x,y_model,color="green",linewidth=2.0)
    y_ogcm=x*np.mean(slope_ogcms[:,iflag])
    plt.plot(x,y_ogcm,color="blue",linewidth=2.0)
    plt.xlim(-1,1)
    plt.ylim(-1,1)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.xlabel("Cor(SST,SSS) (OBS)",fontsize=15)
    plt.ylabel("Cor(SST,SSS) (Model)",fontsize=15)
    plt.title(captions[iflag],fontsize=18)
    text='%2.2f'%np.mean(slope_ogcms[:,iflag])
    plt.text(-0.9,0.9,"Slope(OFES)="+str(text),fontsize=12,color="blue")
    text='%2.2f'%np.mean(slope_models[:,iflag])
    plt.text(-0.9,0.79,"Slope(CMIP5)="+str(text),fontsize=12,color="green")
plt.tight_layout(rect=[0,0,1,0.96])
if (out_form =="eps"):
    plt.savefig(out_name)
else:
    plt.show()
