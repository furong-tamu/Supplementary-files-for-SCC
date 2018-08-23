% This script is used to create the real and simulated data files in this package.
% These data files have been already included in the subfolder "data" in
% this package. So it is not necessary to run this script to make the data
% files available.

%---------------------------------------------------
% simulated data

coef_linear=0.75; % correlation coefficient between the first two covariates
sim_num=100; % simulation numbers

epsilon_mag=0.1; % standard deviation (sigma) for iid random noise
load('cluster_pre.mat');
[x,y]=SCC_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,0.1);
save('cluster_w.mat','lon','lat','x','y','beta');
[x,y]=SCC_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,0.3);
save('cluster_m.mat','lon','lat','x','y','beta');
[x,y]=SCC_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,1);
save('cluster_s.mat','lon','lat','x','y','beta');

load('smooth_pre.mat');
[x,y]=SCC_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,0.1);
save('smooth_w.mat','lon','lat','x','y','beta');
[x,y]=SCC_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,0.3);
save('smooth_m.mat','lon','lat','x','y','beta');
[x,y]=SCC_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,1);
save('smooth_s.mat','lon','lat','x','y','beta');

load('hybrid_pre.mat');
phi=0.3;% range parameter for the spatially correlated noise
epsilon_mag=sqrt(0.5);% standard deviation for spatially correlated noise
[x,y]=SCC_T_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,0.1,phi);
save('hybrid_w.mat','lon','lat','x','y','beta','phi');
[x,y]=SCC_T_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,0.3,phi);
save('hybrid_m.mat','lon','lat','x','y','beta','phi');
[x,y]=SCC_T_simulation_spatial_data_prescribe(beta,epsilon_mag,sim_num,coef_linear,lon,lat,1,phi);
save('hybrid_s.mat','lon','lat','x','y','beta','phi');

disp(['---------------------------------------------------------------------'])
disp(['As the data are randomly generated, the values of x and y'])
disp(['vary by each running of the script.'])

clear

%---------------------------------------------------
% real data
% Please first download the WOA datasets woa13_decav_t00_04v2.nc and
% woa13_decav_s00_04v2.nc from the following URLs
%https://data.nodc.noaa.gov/thredds/catalog/woa/WOA13/DATAv2/temperature/netcdf/decav/0.25/catalog.html?dataset=woa/WOA13/DATAv2/temperature/netcdf/decav/0.25/woa13_decav_t00_04v2.nc
%https://data.nodc.noaa.gov/thredds/catalog/woa/WOA13/DATAv2/salinity/netcdf/decav/0.25/catalog.html?dataset=woa/WOA13/DATAv2/salinity/netcdf/decav/0.25/woa13_decav_s00_04v2.nc

% The WOA dataset is a three-dimensional (lon, lat, and depth) global 
% dataset containing 1440720102=105,753,600 data points. In this study, 
% we select a subsample of data points (10000 points in total) along the segment 
% 25W, 60S-60N to estimate the T-S relationship based on this dataset. 

temp=ncread('/Users/jingzhao/localdata/BR/data/woa13_decav_t00_04v2.nc','t_an');% read temperatue 
salt=ncread('/Users/jingzhao/localdata/BR/data/woa13_decav_s00_04v2.nc','s_an');% read salinity
lon=ncread('/Users/jingzhao/localdata/BR/data/woa13_decav_t00_04v2.nc','lon');% read longitude
lat=ncread('/Users/jingzhao/localdata/BR/data/woa13_decav_t00_04v2.nc','lat');% read latitude
depth=ncread('/Users/jingzhao/localdata/BR/data/woa13_decav_t00_04v2.nc','depth');% read depth
lat=double(lat);
depth=double(depth);
temp=squeeze(temp(620,:,:)); % take temperature record along 25W, 90S-90N
salt=squeeze(salt(620,:,:)); % take salinity record along 25W, 90S-90N
save('data_Figure5.mat','temp','salt','lat','depth'); %creat data_Figure5.mat

lat=lat(120:601);
temp=temp(120:601,:);% take temperature record along 25W, 60S-60N
salt=salt(120:601,:);% take salinity record along 25W, 60S-60N

latn=lat/(max(lat)-min(lat));%transform into nondimensional coordinate
depthn=depth/max(depth);%transform into nondimensional coordinate
[depthn,latn]=meshgrid(depthn,latn);


index_sea=~isnan(salt);% find the locations in the sea
temp=temp(index_sea);
salt=salt(index_sea);
latn=latn(index_sea);
depthn=depthn(index_sea);

n=length(temp);
p=randperm(n,10000); % randomly select 10000 locations
temp=temp(p);
salt=salt(p);
latn=latn(p);
depthn=depthn(p);
save('TS_data.mat','temp','salt','latn','depthn'); %create TS_data.mat
disp(['---------------------------------------------------------------------'])
disp(['As the 10000 locations are randomly selected, the values of variables'])
disp(['in TS_data.mat varies by each running of the script.'])

clear