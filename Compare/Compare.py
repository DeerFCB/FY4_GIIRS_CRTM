# -*- coding: utf-8 -*-
"""
Created on Wed Aug 14 16:32:32 2019

@author: Deer
"""


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
    MW_data_path = r"D:\wxzx\testdata\test0604txt\20190604_MW_data.txt"
    MW_real_path = r"D:\wxzx\testdata\test0604txt\20190604_MW_real.txt"
    crtm_path = r"D:\wxzx\testdata\test0604txt\20190604_crtm_giirs.txt"
    channel_path = r'D:\wxzx\testdata\chan_selected_fy4a_giirs.dat'
    
    channel = pd.read_csv(channel_path,header = None)
    channel = np.array(channel[channel>=690].dropna()-689, dtype = int)    
    
    MW_data = read_MW_data(MW_data_path)
    MW_real =read_MW_real(MW_real_path, MW_data.shape[0], channel)
    crtm_real = read_crtm(crtm_path, MW_data.shape[0], channel)
    crtm_MW_error = np.array([MW_real, crtm_real, crtm_real - MW_real])
    
    MW_data_lat, crtm_MW_error_lat = MW_div_10lat(MW_data, crtm_MW_error)
    
    mean_lat, mean_lat, mean_lon, mean_scan, mean_crtm_MW_error = \
    mean_MW_crtm_lat(MW_data_lat, crtm_MW_error_lat)
    
    
    
    
'''
#所有通道的误差
    channel = np.linspace(1,len(channel),len(channel))
    plt.figure(1)
    for i in range(0,len(mean_crtm_MW_error)):
        plt.plot(channel, mean_crtm_MW_error[i]) 
    plt.show()

#一个通道，横坐标为扫描角，按纬度带画线
    draw_channel = 10
    plt.figure(2)
    for i in range(0,len(mean_crtm_MW_error)):
        temp = np.array([crtm_MW_error_lat[i][2,:,draw_channel-1], MW_data_lat[i].iloc[:,3]]).T
        temp = temp[np.lexsort(temp.T)]
        plt.plot(temp[:,1], temp[:,0], label='%dN'%(10*(i+1)))
    plt.legend(loc='upper right')
    plt.show()
    
    k = 1
    for j in range (1,240, 60 ):
        plt.subplot(2,2,k)
        draw_channel = j
        for i in range(0,len(mean_crtm_MW_error)):
            temp = np.array([crtm_MW_error_lat[i][2,:,draw_channel-1], MW_data_lat[i].iloc[:,3]]).T
            temp = temp[np.lexsort(temp.T)]
            plt.scatter(temp[:,1], temp[:,0], s=5, alpha = .5, label='%dN'%(10*(i+1)))
        k+=1
        plt.legend(loc='upper left')
        plt.title('channel%d'%j)
    plt.show()
    
#单一通道的误差空间分布    
 #   draw_channel = 1
    sub_lola = np.array([MW_data.iloc[:,1].min(),MW_data.iloc[:,1].max(),MW_data.iloc[:,0].min(),MW_data.iloc[:,0].max()])
    map1 = Basemap(llcrnrlon=sub_lola[0],llcrnrlat=sub_lola[2],urcrnrlon=sub_lola[1],urcrnrlat=sub_lola[3],resolution='i', \
             projection='tmerc',lat_0 = (sub_lola[2]+sub_lola[3])/2, lon_0 = (sub_lola[0]+sub_lola[1])/2)
    map1.drawcoastlines()
    map1.scatter(np.array(MW_data.iloc[:,1]), np.array(MW_data.iloc[:,0]), cmap = 'gray', c=crtm_MW_error[2,:,draw_channel-1],s = 15, alpha = .5, latlon=True)
    plt.colorbar()
    plt.show()
    
    
    
    k = 1
    for j in range (1,240, 60 ):
        plt.subplot(2,2,k)
        draw_channel = j
        sub_lola = np.array([MW_data.iloc[:,1].min(),MW_data.iloc[:,1].max(),MW_data.iloc[:,0].min(),MW_data.iloc[:,0].max()])
        map1 = Basemap(llcrnrlon=sub_lola[0],llcrnrlat=sub_lola[2],urcrnrlon=sub_lola[1],urcrnrlat=sub_lola[3],resolution='i',
                 projection='tmerc',lat_0 = (sub_lola[2]+sub_lola[3])/2, lon_0 = (sub_lola[0]+sub_lola[1])/2)
        map1.drawcoastlines()
        map1.scatter(np.array(MW_data.iloc[:,1]), np.array(MW_data.iloc[:,0]), c=crtm_MW_error[2,:,draw_channel-1],s = 10, latlon=True)
        plt.title('channel%d'%j)
        k+=1
    plt.show()
    
'''    


    