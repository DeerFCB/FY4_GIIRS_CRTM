# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 09:37:15 2019

@author: Deer
"""

import os
import netCDF4 as nc
from netCDF4 import Dataset
import numpy as np
import sys
sys.path.append(r'C:\Users\Deer\Desktop\wxzx')
from C_N_F import*


def grib22nc(wgrib_path, ncep_path, ncep_nc_name):
    #grib2转nc格式
    #速度很慢
    order = wgrib_path + ' %s'%(ncep_path) + " -header -netcdf %s"%(ncep_nc_name)
    os.system(order)
    return

def inter2_2(ncep_1,ncep_2):
    ncep_1_1 = ncep_1+(ncep_2 - ncep_1)/3
    ncep_1_2 = ncep_1+2*(ncep_2 - ncep_1)/3
    return ncep_1_1, ncep_1_2

def write_nc(output_path, ncep_lt, ncep_u, ncep_v, ncep_tmp_arr, ncep_q_arr, ncep_o3_ppmv_arr, ncep_lat, ncep_lon):
    print('start writing:'+output_path)
    da = Dataset(output_path, 'w')
    da.createDimension('lons',720) #创建坐标点
    da.createDimension('lats',361)
    da.createDimension('levs',31)
    da.createDimension('times',1)
    da.createVariable("lon",'f',("lons"))
    da.createVariable("lat",'f',("lats"))
    #da.createVariable("levs",'f',("levs"))
    #da.createVariable("times",'f',("times"))
    da.createVariable("ncep_lt",'f8',("times",'lats','lons'))
    da.createVariable("ncep_u",'f8',("times",'lats','lons'))
    da.createVariable("ncep_v",'f8',("times",'lats','lons'))
    da.createVariable("ncep_tmp",'f8',("levs",'lats','lons'))
    da.createVariable("ncep_q",'f8',("levs",'lats','lons'))
    da.createVariable("ncep_o3_ppmv",'f8',("levs",'lats','lons'))
    da.variables['lat'][:]=ncep_lat  #填充数据
    da.variables['lon'][:]=ncep_lon
    da.variables['ncep_lt'][:]=ncep_lt
    da.variables['ncep_u'][:]=ncep_u
    da.variables['ncep_v'][:]=ncep_v
    da.variables['ncep_tmp'][:]=ncep_tmp_arr
    da.variables['ncep_q'][:]=ncep_q_arr
    da.variables['ncep_o3_ppmv'][:]=ncep_o3_ppmv_arr
    #da.variables['levs'][:]=ncep_lon1
    #da.variables['times'][:]=ncep_lon1
    da.close()
    print('finish')

 

        
if __name__ == '__main__':
    
    txt_path = r"D:\wxzx\testdata\test0604txt"+"\\"

            
    #气压
    pressure_31 = np.array([
            1.0000, 2.0000, 3.0000, 5.0000, 7.0000, 10.0000, 20.0000, 30.0000, 
            50.0000, 70.0000, 100.0000, 150.0000, 200.0000, 250.0000, 300.0000,
            350.0000, 400.0000, 450.0000, 500.0000, 550.0000, 600.0000,
            650.0000, 700.0000, 750.0000, 800.0000, 850.0000, 900.0000, 
            925.0000, 950.0000, 975.0000, 1000.0000])
    pressure_30 = (pressure_31[0:30]+pressure_31[1:])/2
    
    #文件路径
    #fy4_path = r"D:\wxzx\2\GIIRS\FY4A-_GIIRS-_N_REGX_1047E_L1-_IRD-_MULT_NUL_20190604000000_20190604001044_016KM_019V1.HDF"
    date_time ='060400_'
    #wgrib2.exe所在路径
    wgrib_path = r"D:\wxzx\wgrib2\wgrib2.exe"
    #grb2文件路径
    ncep_grb2_path = r'D:\wxzx\testdata\test0604_ncep'+'\\'
    #输出转nc文件路径
    ncep_nc_folder = r"D:\wxzx\testdata\test0604_ncep_nc"+"\\"
    #输出插值后文件路径
    output_folder = r'D:\wxzx\testdata\test0604_ncep_nc_inter'+'\\'

    #根据经纬度切片
    sub_lola = [0,360,0,180]
    sub_lola_ind = [0,720,0,361]
    print("范围位于：lat:%dN - %dN, lon:%dE - %dE" %(sub_lola[2], sub_lola[3],sub_lola[0],sub_lola[1]))
    
    
    ncep_paths = []
    for root,dirs,files in os.walk(ncep_grb2_path): 
        for dir in dirs: 
            ncep_paths.append(os.path.join(root,dir))
        for file in files: 
            ncep_paths.append(os.path.join(root,file)) 
    
    #NCEP数据格式转换
    #grin转nc    
    flag = 1
    if (flag==1):
        ncep_nc_paths = []
        for i in range(0,len(ncep_paths)): 
            print(i)
            ncep_nc_paths.append(ncep_nc_folder + date_time + files[i][20:24] + '.nc')
            grib22nc(wgrib_path, ncep_paths[i], ncep_nc_paths[i])
    
    #读#写插值后的nc
    output_paths = []
    output_path = output_folder+date_time + files[0][20:24] + '_or.nc'
    output_paths.append(output_path)
    
    ncep_lt1, ncep_u1, ncep_v1, ncep_tmp_arr1, ncep_q_arr1, ncep_o3_ppmv_arr1, ncep_lat, ncep_lon \
    = read_ncep(ncep_nc_paths[0], pressure_31, sub_lola_ind)
    
    write_nc(output_paths[0], ncep_lt1, ncep_u1, ncep_v1, ncep_tmp_arr1, ncep_q_arr1, ncep_o3_ppmv_arr1, ncep_lat, ncep_lon)

    ncep_lt2, ncep_u2, ncep_v2, ncep_tmp_arr2, ncep_q_arr2, ncep_o3_ppmv_arr2, ncep_lat, ncep_lon \
    = read_ncep(ncep_nc_paths[1], pressure_31, sub_lola_ind)
               
    for i in range(1,len(ncep_nc_paths)): 
    #for i in range(0,5):        
        #网格插值
        lt_1,lt_2 = inter2_2(ncep_lt1,ncep_lt2)
        u_1, u_2 = inter2_2(ncep_u1,ncep_u2)
        v_1, v_2 = inter2_2(ncep_v1,ncep_v2)
        tmp_arr_1, tmp_arr_2 = inter2_2(ncep_tmp_arr1,ncep_tmp_arr2)
        q_arr_1, q_arr_2 = inter2_2(ncep_q_arr1,ncep_q_arr2)
        o3_ppmv_arr_1, o3_ppmv_arr_2 = inter2_2(ncep_o3_ppmv_arr1,ncep_o3_ppmv_arr2)
        
        output_path1 = output_folder+date_time + 'f%03d' %(3*i-2) + '_in.nc'
        output_path2 = output_folder+date_time + 'f%03d' %(3*i-1) + '_in.nc'
        output_path3 = output_folder+date_time + files[i][20:24] + '_or.nc'
        output_paths.append(output_path1)
        output_paths.append(output_path2)
        output_paths.append(output_path3)

        write_nc(output_path1, lt_1, u_1, v_1, tmp_arr_1, q_arr_1, o3_ppmv_arr_1, ncep_lat, ncep_lon)
        write_nc(output_path2, lt_2, u_2, v_2, tmp_arr_2, q_arr_2, o3_ppmv_arr_2, ncep_lat, ncep_lon)
        write_nc(output_path3, ncep_lt1, ncep_u1, ncep_v1, ncep_tmp_arr1, ncep_q_arr1, ncep_o3_ppmv_arr1, ncep_lat, ncep_lon)
        
        
        if i<(len(ncep_nc_paths)-1):
            #更新NCEP背景场
            ncep_lt1, ncep_u1, ncep_v1, ncep_tmp_arr1, ncep_q_arr1, ncep_o3_ppmv_arr1 \
            = ncep_lt2, ncep_u2, ncep_v2, ncep_tmp_arr2, ncep_q_arr2, ncep_o3_ppmv_arr2
            
            ncep_lt2, ncep_u2, ncep_v2, ncep_tmp_arr2, ncep_q_arr2, ncep_o3_ppmv_arr2, ncep_lat, ncep_lon \
            = read_ncep(ncep_nc_paths[i+1], pressure_31, sub_lola_ind)






