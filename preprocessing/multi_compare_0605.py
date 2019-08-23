# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 14:10:53 2019

@author: Deer
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata

def read_MW_data(MW_data_path):
    MW_data =pd.read_csv(MW_data_path, sep ="\s+",header = None)
    return MW_data

def read_MW_real(MW_real_path,point_num, channel):
    MW_real_or = pd.read_csv(MW_real_path, sep ="\s+",header = None)
    MW_real_or = MW_real_or.values
    MW_real_or = MW_real_or.reshape(point_num,961,2)
    MW_real =np.array( MW_real_or[:,channel,1])
    MW_real = MW_real.reshape(MW_real.shape[0],MW_real.shape[1])
    return MW_real

def read_crtm(crtm_path, point_num, channel):
    #读crtm数据
    crtm_data =pd.read_csv(crtm_path, sep ="\s+",header = None)
    crtm_data = crtm_data.values
    crtm_data = crtm_data.reshape(point_num,961,3)
    crtm_real =np.array( crtm_data[:,channel,1])  
    crtm_real = crtm_real.reshape(crtm_real.shape[0],crtm_real.shape[1])
    return crtm_real

def MW_div_10lat(MW_data, crtm_MW_error):
    MW_data_lat = []
    crtm_MW_error_lat = []
    for i in range(0,9):
        index = (MW_data.iloc[:,0]>=i*10) & (MW_data.iloc[:,0]<(i+1)*10)
        if index.any() == True:
            MW_data_lat.append(MW_data[index])
            crtm_MW_error_lat.append(crtm_MW_error[:,index,:])
    
    return MW_data_lat, crtm_MW_error_lat

def mean_MW_crtm_lat(MW_data_lat, crtm_MW_error_lat):
    point_num = []
    mean_lat = []
    mean_lon = []
    mean_scan = []
    mean_crtm_MW_error = []
    for i in range(0,len(MW_data_lat)):
        point_num.append(MW_data_lat[i].shape[0])
        mean_lat.append(MW_data_lat[i].iloc[:,0].mean())
        mean_lon.append(MW_data_lat[i].iloc[:,1].mean())
        mean_scan.append(MW_data_lat[i].iloc[:,3].mean())
        print('平均纬度：%3fN'%mean_lat[i])
        print('平均经度：%3fE'%mean_lon[i])
        print('平均扫描角：%3f'%mean_scan[i])
        print('点总数：%d\n'%point_num[i])
        
        mean_crtm_MW_error.append(crtm_MW_error_lat[i][2].mean(axis=0))
        
    return mean_lat, mean_lat, mean_lon, mean_scan, mean_crtm_MW_error

#def mean_MW_crtm_channel(MW_data_lat, crtm_MW_error_lat):
    

if __name__ == '__main__':
    date = '0605'
    folder =r'D:\wxzx\testdata'+'\\'+date+'_all' +'\\'

    MW_data_folder = folder+'0605_txt\MW_data_txt' +"\\"
    MW_real_folder = folder+'0605_txt\MW_real_txt'+"\\"
    crtm_folder = folder+ "0605_crtm_txt" + "\\"
    channel_path = r'D:\wxzx\testdata\chan_selected_fy4a_giirs.dat'
    
    MW_data_paths = []
    for root,dirs,files in os.walk(MW_data_folder): 
        for file in files: 
            MW_data_paths.append(os.path.join(root,file)) 
            
    MW_real_paths = []
    for root,dirs,files in os.walk(MW_real_folder): 
        for file in files: 
            MW_real_paths.append(os.path.join(root,file))
            
    crtm_paths = []
    for root,dirs,files in os.walk(crtm_folder): 
        for file in files: 
            crtm_paths.append(os.path.join(root,file))
    
    MW_datas = []
    crtm_MW_errors =[]
    for i in range(0,len(crtm_paths)):
    #for i in range(1,2):
        print(files[i])
        channel = pd.read_csv(channel_path,header = None)
        channel = np.array(channel[channel>=690].dropna()-689, dtype = int)    
        
        MW_data = read_MW_data(MW_data_paths[i])
        MW_real =read_MW_real(MW_real_paths[i], MW_data.shape[0], channel)
        crtm_real = read_crtm(crtm_paths[i], MW_data.shape[0], channel)
        crtm_MW_error = np.array([MW_real, crtm_real, crtm_real - MW_real])
        
        MW_datas.append(MW_data)
        crtm_MW_errors.append(crtm_MW_error)
        
        #MW_data_lat, crtm_MW_error_lat = MW_div_10lat(MW_data, crtm_MW_error)
        #MW_data_lats.append(MW_data_lat)
        #crtm_MW_error.append(crtm_MW_error)
        #mean_lat, mean_lat, mean_lon, mean_scan, mean_crtm_MW_error = \
        #mean_MW_crtm_lat(MW_data_lat, crtm_MW_error_lat)
   
    for j in range(2,247,40):    
        plt.figure(j)
        draw_channel = j
        for i in range(len(MW_datas)-1):
            plt.subplot(3,3,i+1)
            sub_lola = np.array([MW_data.iloc[:,1].min(),MW_data.iloc[:,1].max(),MW_data.iloc[:,0].min(),MW_data.iloc[:,0].max()])
            map1 = Basemap(llcrnrlon=sub_lola[0],llcrnrlat=sub_lola[2],urcrnrlon=sub_lola[1],urcrnrlat=sub_lola[3],resolution='i', \
                     projection='tmerc',lat_0 = (sub_lola[2]+sub_lola[3])/2, lon_0 = (sub_lola[0]+sub_lola[1])/2)
            map1.drawcoastlines()
            map1.scatter(np.array(MW_datas[i].iloc[:,1]), np.array(MW_datas[i].iloc[:,0]), c=crtm_MW_errors[i][2,:,draw_channel-1],s = 15, alpha = .5, latlon=True)
            
            #map1.scatter(np.array(MW_datas[i].iloc[:,1]), np.array(MW_datas[i].iloc[:,0]),s = 15, alpha = .5, latlon=True)
            plt.colorbar()
        plt.title('channel%d'%j)
        plt.savefig(r'D:\wxzx\testdata\0605_all\pic'+'\\'+str(j)+'.png', dpi=300)
        plt.show()
    
    
    '''

    plt.figure(4)
    plt.pcolor(z[2,:,:]) 
    plt.show()
    '''
    
    
    
    