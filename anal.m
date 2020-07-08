clear all; clc
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/SAT_data_related_CODE')
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/TCA')
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/mapping_code/')
temp_lat=59.875:-0.25:-59.875;
temp_lon=-179.875:0.25:179.875;
[GLDAS_lon,GLDAS_lat]=meshgrid(temp_lon,temp_lat);
%%
 load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_X_merge_compare_flag.mat')
 load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_c1_merge_compare_flag.mat')
 load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_c2_merge_compare_flag.mat')
% load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_M.mat')
% load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_avg.mat')
%%
M=load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_M.mat');
M_flag=load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_M_flag.mat');
avg=load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_avg.mat');
avg_flag=load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_avg_flag.mat');
%%
M_m=matfile('/project/hydrosense/matlab/mat/AMSR2/Merged/AMSR2_025.mat')

%%
i=17
tt=M.fMSE_M_AM.y;
[a,b]=find(~isnan(tt));
[a(i), b(i)]
[M_flag.fMSE_M_AM_flag.y(a(i),b(i)),M.fMSE_M_AM.y(a(i),b(i))]
%%
nanmean(M.fMSE_M_AM.y(:))
nanmean(M_flag.fMSE_M_AM_flag.y(:))
%%

t=M_flag.fMSE_M_AM_flag.y - fMSE_X_AM_flag.y;

%%
%
t(1)=1;
t(2)=-1;
Statistic_Mapping_NDVI(GLDAS_lat, GLDAS_lon, t, min(t(:)),max(t(:)));
%%
histogram(t(:))