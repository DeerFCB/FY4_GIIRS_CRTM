# FY4_GIIRS_CRTM
通过CRTM模拟FY4GIIRS的MW波段的radiance， 并于观测数据进行比较  
ncep数据（gfs）的时间分辨率为三小时，FY4两小时完成一个周期的扫描共59景图像  
将ncep数据格式由grb2转为nc，并进行线性插值为1小时，将需要的数据输出为.nc  
读写grb2数据使用grib2.exe进行操作  
通过preprocessing/C_N_F.py中的函数进行预处理和计算，multi_为批处理文件夹下的数据，输出为crtm需要的txt文件  
