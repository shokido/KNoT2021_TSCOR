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

dtrend=True
out_form="X11"
#out_form="eps"

#dtrend=False
lon1=0.0; lon2=359.0; lat1=-60.0; lat2=60.0; area_name="Global"
#lon1=130.0; lon2=240.0; lat1=20.0; lat2=60.0; area_name="KE"
#lon1=250.0; lon2=360.0; lat1=20.0; lat2=60.0; area_name="KE"

dt1=dt.datetime(1800,1,1,0,0,0)
dt2=dt.datetime(2019,12,31,0,0,0)
fflags=["","_move5","_move9","_move13","_move15"]
fflags=["","_move9","_move15"]
fflags=["","_move9"]
fflags=["_move9"]
#fflags=[""]
nfile=len(fflags)
lonname="lon"; latname="lat"
varname_ttend="cor_ttendnhf"
varname_stend="cor_stendfwf"
varname_r1_ttend="cor_nhfnhf"; varname_r2_ttend="cor_ttendttend"
varname_r1_stend="cor_fwffwf"; varname_r2_stend="cor_stendstend"

fname_ttend="../CALC/output_ttendnhf_lag_MOAA_JOFURO3_V1.1_2004_2017_move9.nc"
out_name="figure_3_c_plot_map_cor_sscale_dts_MOAA_JOFURO3_V1.5_move5"+area_name
#fname_ttend="../CALC/output_ttendnhf_lag_RG_JOFURO3_2004_2013.nc"
#out_name="figure_3_c_plot_map_cor_sscale_dts_RG_JOFURO3_move5_"+area_name
nc_in=ncdf.Dataset(fname_ttend,"r")
lags=nc_in.variables["lags"][:]
nc_in.close()
lon,lat,time,cor_tt,miss_tt=select_region_TLL_files([fname_ttend],varname_ttend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)

nlon=len(lon); nlat=len(lat)
cor_ttend_array=np.zeros((nfile+1,nlat,nlon))
cor_stend_array=np.zeros((nfile+1,nlat,nlon))
cor_diff_array=np.zeros((nfile+1,nlat,nlon))

lags_out=0
for i in range(0,nfile):
    fname_ttend="../CALC/output_ttendnhf_lag_MOAA_JOFURO3_V1.1_2004_2017"+fflags[i]+".nc"
    fname_stend="../CALC/output_stendfwf_lag_MOAA_JOFURO3_V1.1_2004_2017"+fflags[i]+".nc"
    #    fname_ttend="../CALC/output_ttendnhf_lag_MOAA_JOFURO3_2004_2013"+fflags[i]+".nc"
    #    fname_stend="../CALC/output_stendfwf_lag_MOAA_JOFURO3_2004_2013"+fflags[i]+".nc"
#    fname_ttend="../CALC/output_ttendnhf_lag_RG_JOFURO3_2004_2013"+fflags[i]+".nc"
#    fname_stend="../CALC/output_stendfwf_lag_RG_JOFURO3_2004_2013"+fflags[i]+".nc"
    cor_ttend,miss_ttend=select_region_TLL_files([fname_ttend],varname_ttend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True)
    cor_stend=select_region_TLL_files([fname_stend],varname_stend,dt1,dt2,lat1,lat2,lon1,lon2)


    index_out=np.where(lags==1)[0][0]
    ntime=12*10
    r1_ttend=select_region_TLL_files([fname_ttend],varname_r1_ttend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
    r2_ttend=select_region_TLL_files([fname_ttend],varname_r2_ttend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
    r1_ttend=r1_ttend[index_out,:,:]
    r2_ttend=r2_ttend[index_out,:,:]
    dof=ntime*(1-r1_ttend*r2_ttend)/(1+r1_ttend*r2_ttend)
    tval=t.ppf(q=0.975,df= dof.astype(int))
    q=(tval**2)/(dof-2)
    r_thres_ttend=np.sqrt(q/(1+q))

    r1_stend=select_region_TLL_files([fname_stend],varname_r1_stend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
    r2_stend=select_region_TLL_files([fname_stend],varname_r2_stend,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=False,return_dims=False)
    r1_stend=r1_stend[index_out,:,:]
    r2_stend=r2_stend[index_out,:,:]
    dof=ntime*(1-r1_stend*r2_stend)/(1+r1_stend*r2_stend)
    tval=t.ppf(q=0.975,df= dof.astype(int))
    q=(tval**2)/(dof-2)
    r_thres_stend=np.sqrt(q/(1+q))
    index_out=np.where(lags==lags_out)[0][0]
    cor_ttend=cor_ttend[index_out,:,:]
    cor_stend=cor_stend[index_out,:,:]

    cor_ttend[np.abs(cor_ttend)<=r_thres_ttend]=np.nan
    cor_stend[np.abs(cor_stend)<=r_thres_stend]=np.nan
    cor_ttend_array[i,:,:]=cor_ttend*(-1.0)
    cor_stend_array[i,:,:]=cor_stend
    cor_diff_array[i,:,:]=-1*cor_ttend-cor_stend

i=nfile
cor_ttend_array[i,:,:]=cor_ttend_array[0,:,:]-cor_ttend_array[i-1,:,:]
cor_stend_array[i,:,:]=cor_stend_array[0,:,:]-cor_stend_array[i-1,:,:]
cor_diff_array[i,:,:]=cor_diff_array[0,:,:]-cor_diff_array[i-1,:,:]
cor_ttend_array=np.ma.masked_invalid(cor_ttend_array)
cor_stend_array=np.ma.masked_invalid(cor_stend_array)
cor_diff_array=np.ma.masked_invalid(cor_diff_array)

wks=Ngl.open_wks(out_form,out_name)
res=Ngl.Resources()
res.cnFillOn=True
res.nglDraw=False
res.nglFrame=False
res.cnLineLabelsOn=False
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
res=setmap(res,lon1,lon2,lat1,lat2)
plot=[]
thicks1=[1.0]*21
thicks1[10]=0.0
dashes1=[0.0]*21
dashes1[0:10]=[1.0]*10
res.cnMonoLineThickness=False
res.cnMonoLineDashPattern=False
res.cnLineThicknesses=thicks1
res.cnLineDashPatterns=dashes1
res.cnLinesOn=False
#res.cnFillMode="RasterFill"
for i in range(0,nfile):
    res.cnFillPalette="MPL_jet"
    res.cnFillPalette="MPL_bwr"
        
#    res.tiMainString="Lags="+str(int(lags[index_out]))
    res=setcint(res,0.0,1.0,0.1)
    res=setcint(res,-1.0,1.0,0.1)
        
    res.tiMainString="TT "+fflags[i]
    plot.append(Ngl.contour_map(wks,cor_ttend_array[i,:,:],res))
    res.tiMainString="ST "+fflags[i]
    plot.append(Ngl.contour_map(wks,cor_stend_array[i,:,:],res))
    res=setcint(res,-1.0,1.0,0.1)
    res=setcint(res,-0.4,0.4,0.04)
    res.cnFillPalette="MPL_bwr"
    res.tiMainString="Diff "+fflags[i]
    plot.append(Ngl.contour_map(wks,cor_diff_array[i,:,:],res))
i=nfile
res.tiMainString="TT "
plot.append(Ngl.contour_map(wks,cor_ttend_array[i,:,:],res))
res.tiMainString="ST "
plot.append(Ngl.contour_map(wks,cor_stend_array[i,:,:],res))
res.tiMainString="Diff "
plot.append(Ngl.contour_map(wks,cor_diff_array[i,:,:],res))

pres=Ngl.Resources()
Ngl.panel(wks,plot[0:nfile*3+3],[nfile+1,3],pres)
Ngl.frame(wks)
Ngl.end()

