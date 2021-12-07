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
out_form1="X11"
out_form2="X11"
#out_form1="eps"
#out_form2="eps"

lon1=0.0; lon2=359.0; lat1=-70.0; lat2=70.0; area_name="Global"
points_lon=[]
points_lat=[]
points_lon.append(150.0); points_lat.append(40.0)
#points_lon.append(144.6); points_lat.append(32.0)
#points_lon.append(300.0); points_lat.append(40.0)
points_lon.append(298.0); points_lat.append(40.0)
#points_lon.append(360.0-145.0); points_lat.append(50.0)
points_lon.append(360.0-145.0); points_lat.append(52.0)

points_lon.append(180.0); points_lat.append(0.0)
#points_lon.append(220.0); points_lat.append(0.0)
points_lon.append(240.0); points_lat.append(-8.0)
points_lon.append(32.0); points_lat.append(-45.0)
points_lon.append(110.0); points_lat.append(-30.0)
points_lon.append(320.0); points_lat.append(-40.0)

dt1=dt.datetime(1800,1,1,0,0,0)
dt2=dt.datetime(2019,12,31,0,0,0)
lonname="lon"; latname="lat"
varname_tt="cor_tt";varname_ts="cor_ts";varname_ss="cor_ss"
data_name="MOAA"
fname_in_1="out_ts_lag_"+data_name+"_2004_2019.nc"
out_name1="figure1_map_cor_ts_"+data_name+"_2004_2019"
out_name2="figure1_box_cor_ts_"+data_name+"_2004_2019"

nc_in=ncdf.Dataset(fname_in_1,"r")
lags=nc_in.variables["lags"][:]
nc_in.close()

lon,lat,time,cor_tt,miss_tt=select_region_TLL_files([fname_in_1],varname_tt,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True,return_dims=True)

cor_ts,miss_ts=select_region_TLL_files([fname_in_1],varname_ts,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True)

cor_ss,miss_ss=select_region_TLL_files([fname_in_1],varname_ss,dt1,dt2,lat1,lat2,lon1,lon2,return_missing=True)
cor_ts=np.ma.masked_invalid(cor_ts)



lags_out=[0]
nlag=len(lags_out)
index_out=np.where(lags==1)[0][0]
r1=cor_tt[index_out,:,:]
r2=cor_ss[index_out,:,:]
ntime=12*(2019-2004+1)
dof=ntime*(1-r1*r2)/(1+r1*r2)
tval=t.ppf(q=0.975,df= dof.astype(int))
q=(tval**2)/(dof-2)
r_thres=np.sqrt(q/(1+q))
cor_ts_t=np.copy(cor_ts)
cor_ts_t[np.abs(cor_ts)<r_thres]=np.nan
cor_ts_t=np.ma.masked_invalid(cor_ts_t)

#ntime*(1-r1*r2)/(1+r1*r2)
nlag=len(lags)
npoint=len(points_lat)
cor_tt_point=np.zeros((npoint,nlag))
cor_ts_point=np.zeros((npoint,nlag))
cor_ss_point=np.zeros((npoint,nlag))
rs_point=np.zeros((npoint,nlag))
for i in range(0,npoint):
    a=np.abs(lon-points_lon[i]); ilon=np.min(np.where(a==np.min(a)))
    a=np.abs(lat-points_lat[i]); ilat=np.min(np.where(a==np.min(a)))
    print(ilon)
    print(ilat)
    cor_tt_point[i,:]=cor_tt[:,ilat,ilon]
    cor_ts_point[i,:]=cor_ts[:,ilat,ilon]
    cor_ss_point[i,:]=cor_ss[:,ilat,ilon]
    rs_point[i,:]=r_thres[ilat,ilon]
    
print(cor_ts_point)
wks1=Ngl.open_wks(out_form1,out_name1)
wks2=Ngl.open_wks(out_form2,out_name2)
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
res=setcint(res,-1.0,1.0,0.1)
#res.cnFillMode="RasterFill"
plot=[]

index_out=np.where(lags==0)[0][0]
res.tiMainString="Lags="+str(int(lags[index_out]))
res.lbOrientation="Vertical"
res.pmLabelBarHeightF=0.6
res.pmLabelBarWidthF=0.06
plot.append(Ngl.contour_map(wks1,cor_ts_t[index_out,:,:],res))


gsres                   = Ngl.Resources()
# Polymarker resources.
gsres.gsMarkerIndex     = 16       # dots
gsres.gsMarkerColor     = "Black"
gsres.gsFillLineThicknessF=3.0
gsres.gsMarkerSizeF     = 0.014    # twice normal size
prim6 = Ngl.add_polymarker(wks1,plot[0],points_lon,points_lat,gsres)

pres=Ngl.Resources()
Ngl.panel(wks1,plot[0:1],[1,1],pres)
#Ngl.frame(wks1)

print("s")
cres=Ngl.Resources()
cres.nglDraw=False
cres.nglFrame=False
cres.xyLineThicknessF=3.0
cres.trYMaxF=1.0
cres.trYMinF=-1.0
cres.tiXAxisString="Lag(Month)"
#cres.tmYMajorGrid= True
#cres.tmXMajorGrid= True
cres.tmYROn=False
cres.tmXTOn=False
cres.tmXBLabelFontHeightF=0.035
cres.tmYLLabelFontHeightF=0.035
cres.tiMainFontHeightF=0.035
lplot1=[]
lplot2=[]
lplot3=[]
lplot4=[]
lplot5=[]
for i in range(0,npoint):
    cres.xyDashPattern=0
    print(i)
    cres.tiMainString=str(points_lon[i])+" "+str(points_lat[i])
    cres.xyLineColor="Black"
    lplot1.append(Ngl.xy(wks2,lags,cor_ts_point[i,:],cres))
#    cres.xyLineColor="Blue"
#    lplot2.append(Ngl.xy(wks2,lags,cor_tt_point[i,:],cres))
#    cres.xyLineColor="Red"
#    lplot3.append(Ngl.xy(wks2,lags,cor_ss_point[i,:],cres))
    cres.xyLineColor="Grey"
    cres.xyDashPattern=2
    lplot4.append(Ngl.xy(wks2,lags,rs_point[i,:],cres))
    lplot5.append(Ngl.xy(wks2,lags,-1*rs_point[i,:],cres))
#    Ngl.overlay(lplot1[i],lplot2[i])
#    Ngl.overlay(lplot1[i],lplot3[i])
    Ngl.overlay(lplot1[i],lplot4[i])
    Ngl.overlay(lplot1[i],lplot5[i])
print("a")
Ngl.panel(wks2,lplot1[0:npoint],[npoint/2,2],pres)
Ngl.frame(wks2)
Ngl.end()

