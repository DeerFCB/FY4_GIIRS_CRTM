# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 16:19:03 2019

@author: Deer
"""

import os
import netCDF4 as nc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import math

'''
读FY4数据
需要：
ES_RealLW
IRLW_Latitude，IRLW_Longitude，
IRLW_SatelliteScanAngle,IRLW_SatelliteZenith,
ES_RealMW
IRMW_Latitude，IRMW_Longitude，
IRMW_SatelliteScanAngle,IRMW_SatelliteZenith,

'''
def read_fy4(fy4_path, cloud_value):
    '''
    input:
        filename(path)
    output:
        LW_data: 'IRLW_Latitude','IRLW_Longitude','IRLW_SatelliteAzimuth','IRLW_SatelliteZenith'
        MW_data: 'IRMW_Latitude','IRMW_Longitude','IRMW_SatelliteAzimuth','IRMW_SatelliteZenith'
        'ES_RealLW',
        'ES_RealMW',
    '''
    #将数据读为字典格式
    with nc.Dataset(fy4_path) as file:
        file.set_auto_mask(False)  # 可选
        fy4_data_ori = {x: file[x][()] for x in file.variables}
    LW_index = {'IRLW_Latitude','IRLW_Longitude','IRLW_SatelliteAzimuth','IRLW_SatelliteZenith'} 
    MW_index = {'IRMW_Latitude','IRMW_Longitude','IRMW_SatelliteAzimuth','IRMW_SatelliteZenith'} 

    #建立子字典
    LW_dic = {key:value for key,value in fy4_data_ori.items()if key in LW_index}     
    MW_dic = {key:value for key,value in fy4_data_ori.items()if key in MW_index}  
    #LW_index = {'LW_Lat','LW_Lon','LW_Azi','LW_Zen'} 
    #MW_index = {'MW_Lat','MW_Lon','MW_Azi','MW_Zen'}
    
    LW_data = pd.DataFrame(LW_dic)
    MW_data = pd.DataFrame(MW_dic)
    #LW_data.info(),将bigending转为float
    LW_data = LW_data*1
    MW_data = MW_data*1
    #计算扫描角
    R = 6378.1
    h = 35790
    LW_data['IRLW_ScanAngle'] = np.arcsin(R/(R+h)*np.sin(LW_data['IRLW_SatelliteZenith']*2*np.pi/360))/np.pi*180
    MW_data['IRMW_ScanAngle'] = np.arcsin(R/(R+h)*np.sin(MW_data['IRMW_SatelliteZenith']*2*np.pi/360))/np.pi*180
    #LW,MW的DN值
    LW_real = np.array(fy4_data_ori['ES_RealLW'])
    MW_real = np.array(fy4_data_ori['ES_RealMW'])
    #可见光通道
    VIS_Lat = np.array(fy4_data_ori['VIS_Latitude'])
    VIS_Lon = np.array(fy4_data_ori['VIS_Longitude'])
    VIS_dn_ori = np.array(fy4_data_ori['ES_ContVIS'])
    Cal_Ta = np.array(fy4_data_ori['ES_CalSTableVIS'])
    VIS_dn = Vis_Cal(Cal_Ta, VIS_dn_ori)
    
    for i in range(0,len(MW_data)):
        #没有经纬度的点设置空值
        if LW_data.iloc[i,0]<-90:
            LW_data.iloc[i,:] =np.nan
            LW_real[:,i] = np.nan
        if MW_data.iloc[i,0]<-90:
            MW_data.iloc[i,:] =np.nan
            MW_real[:,i] = np.nan
        #天顶角大于80度设置空值
        if MW_data.iloc[i,3]>80:
            MW_data.iloc[i,:] =np.nan
            MW_real[:,i] = np.nan
        if MW_data.iloc[i,3]>80:
            MW_data.iloc[i,:] =np.nan
            MW_real[:,i] = np.nan
            
    #删除空值点        
    LW_nan_index = np.array(LW_data.index[LW_data.isnull().T.any()])
    MW_nan_index = np.array(MW_data.index[MW_data.isnull().T.any()])
    LW_data.drop( LW_nan_index , inplace = True)
    MW_data.drop( MW_nan_index , inplace = True)
    LW_real = np.delete(LW_real, LW_nan_index, axis = 1)
    MW_real = np.delete(MW_real, MW_nan_index, axis = 1)
   
    if MW_data['IRMW_Latitude'].empty :
        return LW_data, MW_data,LW_real,MW_real
    #通过可见光判断是否有云
    MW_VIS_cloud, LW_VIS_cloud = check_cloud(VIS_Lon, VIS_Lat, VIS_dn, MW_data, LW_data, Cal_Ta, cloud_value)
    
    #删除有云点
    LW_cloud_index = np.array(LW_VIS_cloud.index[LW_VIS_cloud.isnull().T.any()])
    MW_cloud_index = np.array(MW_VIS_cloud.index[MW_VIS_cloud.isnull().T.any()])
    LW_data.drop( LW_data.iloc[LW_cloud_index].index , inplace = True)
    MW_data.drop( MW_data.iloc[MW_cloud_index].index , inplace = True)

    LW_real = np.delete(LW_real, LW_cloud_index, axis = 1)
    MW_real = np.delete(MW_real, MW_cloud_index, axis = 1)
    
    return LW_data, MW_data,LW_real,MW_real

def Vis_Cal(Cal_Ta, VIS_dn_ori):
    #根据转换表，转换为dn值
    VIS_dn = np.zeros(shape = VIS_dn_ori.shape, dtype = 'float32')
    for i in range(0, VIS_dn_ori.shape[0]):
        for j in range(0, VIS_dn_ori.shape[1]):
            if VIS_dn_ori[i,j]>=83 and VIS_dn_ori[i,j]<=951:
                VIS_dn[i,j] = Cal_Ta[VIS_dn_ori[i,j]]
            else:
                VIS_dn[i,j] = 65535.0 
    VIS_dn[VIS_dn == 65535] = np.nan        
    return VIS_dn

def check_cloud(VIS_Lon, VIS_Lat, VIS_dn, MW_data, LW_data, Cal_Ta, cloud_value):
    #通过可见光判断是否有云
    VIS_points = np.array([VIS_Lon.flatten(), VIS_Lat.flatten()]).T
    MW_VIS = griddata( VIS_points, VIS_dn.flatten(), (MW_data['IRMW_Longitude'],MW_data['IRMW_Latitude']))
    LW_VIS = griddata( VIS_points, VIS_dn.flatten(), (LW_data['IRLW_Longitude'],LW_data['IRLW_Latitude']))
    #MW_VIS_cloud: MW_VIS<=Cal_Ta[cloud_value]晴空，1
    MW_VIS_cloud = pd.DataFrame(np.where(MW_VIS <= Cal_Ta[cloud_value], 1, np.nan))
    LW_VIS_cloud = pd.DataFrame(np.where(LW_VIS <= Cal_Ta[cloud_value], 1, np.nan))
    
    return MW_VIS_cloud, LW_VIS_cloud
'''
读NCEP数据
需要：
?10m风场？
u风场(UGRD_10maboveground)
v风场(VGRD_10maboveground)

地面温度TMP_surface

?用31行平均出30个值？
level_pressure
pressure

不同层TMP(TMP_100mb)
Absorber(H2O)
Absorber(O3MR_100mb)
'''

'''
将当前目录转到wgrib工具目录
os.chdir(r"D:\wxzx\wgrib2")
执行cmd命令，查看具体信息
os.system(r"wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -v")
os.system("wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -d l -text -nh -o test.txt")
取出第n个， 转为.csv
os.system(r"wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -d 2 -header -csv test.csv")
转.nc文件
os.system(r"wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -header -netcdf test.nc")
使用管道
windows
wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 | findstr "TMP" | wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -lola 100:100:0.5 20:100:0.5 out.csv spread  -i
linux
wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 | gred 'TMP'
os.system(r"wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -header -netcdf test.nc")
wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -match "TMP" -match "mb"  -not "above" | wgrib2 D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000 -lola 100:100:0.5 20:100:0.5 out.csv spread  -i

wgrib_path = r"D:\wxzx\wgrib2\wgrib2.exe"
ncep_path = r"D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000"
ncep_nc_path = r"D:\wxzx\testdata\test01.nc"
'''


def grib22nc(wgrib_path, ncep_path, ncep_nc_name):
    #grib2转nc格式
    #速度很慢
    order = wgrib_path + ' ' +  ncep_path + " -header -netcdf " + ncep_nc_name
    os.system(order)
    return

def subgrid(LW_data, MW_data):
    #根据卫星的经纬度切割NCEP数据
    '''
    N_lat = math.ceil(max(LW_data['IRLW_Latitude'].max(),MW_data['IRMW_Latitude'].max()))
    S_lat = math.floor(min(LW_data['IRLW_Latitude'].min(),MW_data['IRMW_Latitude'].min()))
    E_lon = math.ceil(max(LW_data['IRLW_Longitude'].max(),MW_data['IRMW_Longitude'].max()))
    W_lon = math.floor(min(LW_data['IRLW_Longitude'].min(),MW_data['IRMW_Longitude'].min()))
    '''
    N_lat = math.ceil(MW_data['IRMW_Latitude'].max())
    S_lat = math.floor(MW_data['IRMW_Latitude'].min())
    E_lon = math.ceil(MW_data['IRMW_Longitude'].max())
    W_lon = math.floor(MW_data['IRMW_Longitude'].min())
    
    N_lat_ind = 2*(N_lat+90)+1
    S_lat_ind = 2*(S_lat+90)
    E_lon_ind = 2*E_lon+1
    W_lon_ind = 2*W_lon
    
    sub_lola = [W_lon, E_lon, S_lat, N_lat]
    sub_lola_ind = [W_lon_ind, E_lon_ind, S_lat_ind, N_lat_ind]
    return sub_lola, sub_lola_ind

def read_ncep_tmp(ncep_data_ori, sub_lola_ind):
    #读气温
    #筛选key
    tmp_key = []
    ncep_tmp_dic = {}
    ncep_tmp_arr = np.zeros(shape = (31, 361, 720))
    i = 0
    for key in ncep_data_ori.keys():
        if ('TMP' in key and 'mb' in key and 'above' not in key):
            tmp_key.append(key)
            ncep_tmp_dic[key] = ncep_data_ori[key]
            ncep_tmp_arr[i,:,:] =  ncep_data_ori[key]
            i += 1
    #ncep_tmp_dic = {key:value for key,value in ncep_data_ori.items()if key in tmp_key} 
    ncep_tmp_arr_sub = ncep_tmp_arr[:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]]
    return ncep_tmp_dic,ncep_tmp_arr_sub

def read_ncep_rh(ncep_data_ori, sub_lola_ind):
    #读相对湿度
    rh_key = []
    ncep_rh_dic = {}
    ncep_rh_arr = np.zeros(shape = (31, 361, 720))
    i = 0
    for key in ncep_data_ori.keys():
        if ('RH' in key and 'mb' in key and 'above' not in key):
            rh_key.append(key)
            ncep_rh_dic[key] = ncep_data_ori[key]
            ncep_rh_arr[i,:,:] =  ncep_data_ori[key]
            i += 1    
    ncep_rh_arr_sub = ncep_rh_arr[:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]]
    return ncep_rh_dic,ncep_rh_arr_sub

def cal_rh2q(ncep_tmp_arr, ncep_rh_arr, pressure_31, sub_lola_ind):
    '''
    input:
        温度ncep_tmp_arr
        相对湿度ncep_rh_arr
        31层气压pressure_31
    output:
        比湿ncep_q_arr
    note:
        Tetens经验公式计算
        饱和水汽压es
        实际水汽压e
        比湿q
    '''
    p = pressure_31.reshape(31,1,1)*np.ones(shape =(31,361,720))[:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]]
   # es = 6.1078*(np.exp(17.2693882*(ncep_tmp_arr-273.16)/(ncep_tmp_arr-35.86)))
  #  '''
    es = 10**(10.79574*(1-273.16/ncep_tmp_arr) - 
              5.028*np.log10(ncep_tmp_arr/273.16) + 
              1.50475*10**(-4)*(1-10**(-8.2969*(ncep_tmp_arr/273.16-1))) + 
              0.42873*10**(-3)*(10**(4.76955*(1-273.16/ncep_tmp_arr))-1)+0.78614)
  #  '''
    e = es * ncep_rh_arr/100
    ncep_q_arr = 0.622 * e / (p-0.378 * e) * 1000
    return ncep_q_arr

def read_ncep_o3(ncep_data_ori,sub_lola_ind):
    o3_key = []
    ncep_o3_kg_dic = {}
    ncep_o3_kg_arr = np.zeros(shape = (31, 361, 720))
    i = 0
    for key in ncep_data_ori.keys():
        if ("O3MR" in key ):
            o3_key.append(key)
            ncep_o3_kg_dic[key] = ncep_data_ori[key]
            ncep_o3_kg_arr[i,:,:] =  ncep_data_ori[key]
            i += 1
    ncep_o3_kg_arr_sub = ncep_o3_kg_arr[:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]]        
    ncep_o3_ppmv_arr_sub = ncep_o3_kg_arr_sub*(29/48)*(10**(6))
    return ncep_o3_kg_dic,ncep_o3_ppmv_arr_sub  

def read_ncep(ncep_nc_path, pressure_31, sub_lola_ind):
    '''
    input:
        filename(path)
        pressure_31
        sub_lola_ind
    
    '''
    
    
    with nc.Dataset(ncep_nc_path) as file:
        file.set_auto_mask(False)  # 可选
        ncep_data_ori = {x: file[x][()] for x in file.variables}
    #u,v取10m风速
    ncep_u_sub = np.array(ncep_data_ori['UGRD_10maboveground'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_v_sub = np.array(ncep_data_ori['VGRD_10maboveground'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_lt_sub = np.array(ncep_data_ori['TMP_surface'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    
    ncep_tmp_dic,ncep_tmp_arr_sub = read_ncep_tmp(ncep_data_ori, sub_lola_ind)
    
    ncep_rh_dic,ncep_rh_arr_sub = read_ncep_rh(ncep_data_ori, sub_lola_ind)
    ncep_o3_kg_dic,ncep_o3_ppmv_arr_sub = read_ncep_o3(ncep_data_ori,sub_lola_ind)
    
    ncep_q_arr_sub = cal_rh2q(ncep_tmp_arr_sub, ncep_rh_arr_sub, pressure_31, sub_lola_ind)
    ncep_lat_sub = ncep_data_ori["latitude"][sub_lola_ind[2]:sub_lola_ind[3]]
    ncep_lon_sub = ncep_data_ori["longitude"][sub_lola_ind[0]:sub_lola_ind[1]]
    
    return ncep_lt_sub, ncep_u_sub, ncep_v_sub, ncep_tmp_arr_sub, ncep_q_arr_sub, ncep_o3_ppmv_arr_sub, ncep_lat_sub, ncep_lon_sub


def read_ncep_in(ncep_nc_in_path, sub_lola_ind):
    with nc.Dataset(ncep_nc_in_path) as file:
        file.set_auto_mask(False)  # 可选
        ncep_data_ori = {x: file[x][()] for x in file.variables}
    ncep_lt_sub = np.array(ncep_data_ori['ncep_lt'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_u_sub = np.array(ncep_data_ori['ncep_u'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_v_sub = np.array(ncep_data_ori['ncep_v'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_tmp_arr_sub = np.array(ncep_data_ori['ncep_tmp'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_q_arr_sub = np.array(ncep_data_ori['ncep_q'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_o3_ppmv_arr_sub = np.array(ncep_data_ori['ncep_o3_ppmv'][:,sub_lola_ind[2]:sub_lola_ind[3],sub_lola_ind[0]:sub_lola_ind[1]])
    ncep_lat_sub = np.array(ncep_data_ori['lat'][sub_lola_ind[2]:sub_lola_ind[3]])
    ncep_lon_sub = np.array(ncep_data_ori['lon'][sub_lola_ind[0]:sub_lola_ind[1]])
           
    return ncep_lt_sub, ncep_u_sub, ncep_v_sub, ncep_tmp_arr_sub, ncep_q_arr_sub, ncep_o3_ppmv_arr_sub, ncep_lat_sub, ncep_lon_sub
        
        
'''
def main():    
    pressure_31 = np.array([
            1.0000, 2.0000, 3.0000, 5.0000, 7.0000, 10.0000, 20.0000, 30.0000, 
            50.0000, 70.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000,
            350.0000, 400.0000, 450.0000, 500.0000, 550.0000, 600.0000,
            650.0000, 700.0000, 750.0000, 800.0000, 850.0000, 900.0000, 
            925.0000, 950.0000, 975.0000, 1000.0000])
    
    LW_data, MW_data,LW_real,MW_real = read_fy4(fy4_path)
    ncep_lt, ncep_u, ncep_v, ncep_tmp_arr, ncep_q_arr, ncep_o3_ppmv_arr, ncep_lat, ncep_lon = read_ncep(ncep_nc_path, pressure_31)
'''  
    
def flat_lola(ncep_lon_sub, ncep_lat_sub):
    ncep_lat_flat = np.repeat(ncep_lat_sub,len(ncep_lon_sub))
    ncep_lon_flat = np.tile(ncep_lon_sub,len(ncep_lat_sub))
    necp_points = np.array( (ncep_lon_flat, ncep_lat_flat) ).T
    return necp_points


def inter_128_30(necp_points, MW_data, ncep_lt_sub, ncep_u_sub, ncep_v_sub, ncep_tmp_arr_sub, ncep_q_arr_sub, ncep_o3_ppmv_arr_sub):
   
    lt_inter = griddata( necp_points, ncep_lt_sub.flatten(), (MW_data['IRMW_Longitude'],MW_data['IRMW_Latitude']) )
    u_inter = griddata( necp_points, ncep_u_sub.flatten(), (MW_data['IRMW_Longitude'],MW_data['IRMW_Latitude']) )
    v_inter = griddata( necp_points, ncep_v_sub.flatten(), (MW_data['IRMW_Longitude'],MW_data['IRMW_Latitude']) )

    
    k = len(ncep_tmp_arr_sub)
    tmp_inter_31 = np.zeros(shape = (k,len(MW_data)))
    q_inter_31 = np.zeros(shape = (k,len(MW_data)))
    o3_inter_31 = np.zeros(shape = (k,len(MW_data)))
    
    for i in range(0, k):
        tmp_inter_31[i,:] = griddata( necp_points, ncep_tmp_arr_sub[i,:,:].flatten(), (MW_data['IRMW_Longitude'],MW_data['IRMW_Latitude']) )
        q_inter_31[i,:] = griddata( necp_points, ncep_q_arr_sub[i,:,:].flatten(), (MW_data['IRMW_Longitude'],MW_data['IRMW_Latitude']) )
        o3_inter_31[i,:] = griddata( necp_points, ncep_o3_ppmv_arr_sub[i,:,:].flatten(), (MW_data['IRMW_Longitude'],MW_data['IRMW_Latitude']) )
    
    tmp_inter_30 = (tmp_inter_31[0:30,:] + tmp_inter_31[1:,:])/2
    q_inter_30 = (q_inter_31[0:30,:] + q_inter_31[1:,:])/2
    o3_inter_30 = (o3_inter_31[0:30,:] + o3_inter_31[1:,:])/2
    return lt_inter, u_inter, v_inter, tmp_inter_30, q_inter_30, o3_inter_30
        
def write2txt(txt_path, file, MW_data, lt_inter, u_inter, v_inter, pressure_31, pressure_30, tmp_inter_30, q_inter_30, o3_inter_30):
    nt = 30
    temp = 996
    #每个文件写入一个txt
    #with open(txt_path +file[-45:-6]+'.txt', 'w+') as f:
    #写入一个txt
    with open(txt_path +file[-45:-37]+'.txt', 'a+') as f:
        for j in range(0, len(MW_data)):
            f.write("%6d%12.4f%12.4f\n" % (nt, MW_data['IRMW_Latitude'].iloc[j], MW_data['IRMW_Longitude'].iloc[j]))
            #扫描角和高度角
            f.write("%12.4f%12.4f%12.4f%12.4f\n" % (MW_data['IRMW_SatelliteZenith'].iloc[j], MW_data['IRMW_ScanAngle'].iloc[j], temp, temp))
            f.write("%6d%12.4f%12.4f%12.4f%12.4f%12.4f\n" % (j+1, lt_inter[j], temp, u_inter[j], v_inter[j], temp))
            for k in range(0, 31):
                f.write("%12.4f" % (pressure_31[k]))
            f.write("\n")
            for k in range(0, 30):
                f.write("%12.4f" % (pressure_30[k]))
            f.write("\n")
            for k in range(0, 30):
                f.write("%12.4f" % (tmp_inter_30[k, j]))
            f.write("\n")
            for k in range(0, 30):
                f.write("%12.4f" % (q_inter_30[k, j]))
            f.write("\n")
            for k in range(0, 30):
                f.write("%12.4f" % (o3_inter_30[k, j]))
            f.write("\n")
    return 

def write_fy4(txt_path, file, MW_data, MW_real):
    print("write fy4")
    with open(txt_path +file[-45:-37]+'_MW_data.txt', 'a+') as f:
        for j in range(0, len(MW_data)):
        #经度，纬度，高度角，扫描角
            f.write("%12.6f%12.6f%12.6f%12.6f\n" 
                    % (MW_data['IRMW_Latitude'].iloc[j], MW_data['IRMW_Longitude'].iloc[j], 
                       MW_data['IRMW_SatelliteZenith'].iloc[j], MW_data['IRMW_ScanAngle'].iloc[j]))
    with open(txt_path +file[-45:-37]+'_MW_real.txt', 'a+') as f:
        for i in range(0, MW_real.shape[1]):
            for j in range(0, len(MW_real)):
                f.write("%6d%12.6f\n"%(j+1,MW_real[j,i]))
    return

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:                   #判断是否存在文件夹如果不存在则创建为文件夹
        os.makedirs(path)            #makedirs创建文件时如果路径不存在会创建这个路径
        print("---  new folder...  ---")
        print("---  OK  ---") 
    else:
        print("---  There is this folder!  ---")            

    
    
if __name__ == '__main__':
    #文件路径
    #fy4_path = r"D:\wxzx\2\GIIRS\FY4A-_GIIRS-_N_REGX_1047E_L1-_IRD-_MULT_NUL_20190604000000_20190604001044_016KM_019V1.HDF"
    fy4_folder_path = r'D:\wxzx\testdata\test0604'+'\\'
    fy4_paths = []
    for root,dirs,files in os.walk(fy4_folder_path): 
        for dir in dirs: 
            fy4_paths.append(os.path.join(root,dir))
        for file in files: 
            fy4_paths.append(os.path.join(root,file)) 
    
    wgrib_path = r"D:\wxzx\wgrib2\wgrib2.exe"
    ncep_path = r"D:\wxzx\2\gfs.2019060400\gfs.t00z.pgrb2.0p50.f000"
    ncep_nc_path = r'D:\wxzx\testdata\20190604000.nc'
    txt_path = r"D:\wxzx\testdata\test0604txt"+"\\"
    
    #ncep_nc_in_path =  r'D:\wxzx\testdata\test0604_ncep_nc_inter'+'\\'+'060400_f000_or.nc'
    #ncep_nc_in_path =  r'D:\wxzx\testdata\test0604_ncep_nc_inter'+'\\'+'060400_f000_or.nc'
    #气压
    pressure_31 = np.array([
            1.0000, 2.0000, 3.0000, 5.0000, 7.0000, 10.0000, 20.0000, 30.0000, 
            50.0000, 70.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000,
            350.0000, 400.0000, 450.0000, 500.0000, 550.0000, 600.0000,
            650.0000, 700.0000, 750.0000, 800.0000, 850.0000, 900.0000, 
            925.0000, 950.0000, 975.0000, 1000.0000])
    pressure_30 = (pressure_31[0:30]+pressure_31[1:])/2
    #云的阈值
    cloud_value = 110
    
    #NCEP数据格式转换
    #grin转nc    
    flag = 0
    if (flag==1):
        grib22nc(wgrib_path, ncep_path, ncep_nc_path)
    
    for i in range(0,len(fy4_paths)): 
    #for i in range(0,20):    
    
        #读卫星数据
        LW_data, MW_data,LW_real,MW_real = read_fy4(fy4_paths[i], cloud_value)
        if MW_data['IRMW_Latitude'].empty or MW_data['IRMW_Latitude'].isna().any() :
            continue
        
        #根据经纬度切片
        sub_lola, sub_lola_ind = subgrid(LW_data, MW_data)
        print("范围位于：lat:%dN - %dN, lon:%dE - %dE" %(sub_lola[2], sub_lola[3],sub_lola[0],sub_lola[1]))
        
        #读NCEP背景场
        ncep_lt_sub, ncep_u_sub, ncep_v_sub, ncep_tmp_arr_sub, ncep_q_arr_sub, ncep_o3_ppmv_arr_sub, ncep_lat_sub, ncep_lon_sub = read_ncep(ncep_nc_path, pressure_31, sub_lola_ind)
        
        #读插值后的NCEP背景场
        #ncep_lt_sub1, ncep_u_sub1, ncep_v_sub1, ncep_tmp_arr_sub1, ncep_q_arr_sub1, ncep_o3_ppmv_arr_sub1, ncep_lat_sub1, ncep_lon_sub1 = read_ncep_in(ncep_nc_in_path, sub_lola_ind)
        
        #网格插值到观测点
        necp_points = flat_lola(ncep_lon_sub, ncep_lat_sub)
        lt_inter, u_inter, v_inter, tmp_inter_30, q_inter_30, o3_inter_30 = inter_128_30(necp_points, MW_data, ncep_lt_sub, ncep_u_sub, ncep_v_sub, ncep_tmp_arr_sub, ncep_q_arr_sub, ncep_o3_ppmv_arr_sub)
        #写入多个txt文件
        #write2txt(txt_path, files[j], MW_data, lt_inter, u_inter, v_inter, pressure_31, pressure_30, tmp_inter_30, q_inter_30, o3_inter_30)
        #写入一个txt
        write2txt(txt_path, files[i], MW_data, lt_inter, u_inter, v_inter, pressure_31, pressure_30, tmp_inter_30, q_inter_30, o3_inter_30)
        
        write_fy4(txt_path, files[i], MW_data, MW_real)


    