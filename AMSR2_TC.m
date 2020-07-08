clear all; clc
% TC calculations for AMSR2 X,C1,C2, AM/PM with SMAP AM/PM and GLDAS AM/PM
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/SAT_data_related_CODE')
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/TCA')
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/mapping_code/')

ifp_list={'/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/AMSR2/AMSR2_SM_2015_2019_DES.mat',...
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/AMSR2/AMSR2_SM_2015_2019_ASC.mat',...
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/AMSR2/AMSR2_SM_2015_2019_DES_c1.mat',...
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/AMSR2/AMSR2_SM_2015_2019_ASC_c1.mat',...
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/AMSR2/AMSR2_SM_2015_2019_DES_c2.mat',...
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/AMSR2/AMSR2_SM_2015_2019_ASC_c2.mat',...
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/SMAP/SMAP_SM_2015_2019.mat',... % AM
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/SMAP/SMAP_SM_2015_2019_pm.mat',... % PM
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/GLDAS/GLDAS_0_SM_2015_2019.mat',... % DES (1:30 am)
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/GLDAS/GLDAS_12_SM_2015_2019.mat',...; % ASC (13:30 pm)
    '/sfs/qumulo/qproject/hydrosense/matlab/mat/resampled_025/ASCAT/ASCAT_SM_2015_2019.mat'};

AMSR2_X_AM_m=matfile(ifp_list{1});
AMSR2_X_PM_m=matfile(ifp_list{2});

AMSR2_c1_AM_m=matfile(ifp_list{3});
AMSR2_c1_PM_m=matfile(ifp_list{4});

AMSR2_c2_AM_m=matfile(ifp_list{5});
AMSR2_c2_PM_m=matfile(ifp_list{6});

SMAP_AM_m=matfile(ifp_list{7});
SMAP_PM_m=matfile(ifp_list{8});

GLDAS_AM_m=matfile(ifp_list{9});
GLDAS_PM_m=matfile(ifp_list{10});

ASCAT_SM_m=matfile(ifp_list{11});
%%
%cdf
tic
[SNR_X_AM, ~, R_X_AM, fMSE_X_AM, VAR_err_X_AM]=TCbasedNumbers_V3(GLDAS_AM_m.GLDAS_SM, AMSR2_X_AM_m.AMSR2_SM, SMAP_AM_m.SMAP_SM, 2015, 2019,GLDAS_AM_m.GLDAS_SM);
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/X_AM.mat','SNR_X_AM','R_X_AM','fMSE_X_AM','VAR_err_X_AM')
toc

tic
[SNR_X_PM, ~, R_X_PM, fMSE_X_PM, VAR_err_X_PM]=TCbasedNumbers_V3(GLDAS_PM_m.GLDAS_SM, AMSR2_X_PM_m.AMSR2_SM, SMAP_PM_m.SMAP_SM,2015,2019,GLDAS_PM_m.GLDAS_SM);
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/X_PM.mat','SNR_X_PM','R_X_PM','fMSE_X_PM','VAR_err_X_PM')
toc

tic
[SNR_c1_AM, ~, R_c1_AM, fMSE_c1_AM, VAR_err_c1_AM]=TCbasedNumbers_V3(GLDAS_AM_m.GLDAS_SM, AMSR2_c1_AM_m.AMSR2_SM, SMAP_AM_m.SMAP_SM,2015,2019,GLDAS_AM_m.GLDAS_SM);
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/c1_AM.mat','SNR_c1_AM','R_c1_AM','fMSE_c1_AM','VAR_err_c1_AM')
toc

tic
[SNR_c1_PM, ~, R_c1_PM, fMSE_c1_PM, VAR_err_c1_PM]=TCbasedNumbers_V3(GLDAS_PM_m.GLDAS_SM, AMSR2_c1_PM_m.AMSR2_SM, SMAP_PM_m.SMAP_SM,2015,2019,GLDAS_PM_m.GLDAS_SM);
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/c1_PM.mat','SNR_c1_PM','R_c1_PM','fMSE_c1_PM','VAR_err_c1_PM')
toc

tic
[SNR_c2_AM, ~, R_c2_AM, fMSE_c2_AM, VAR_err_c2_AM]=TCbasedNumbers_V3(GLDAS_AM_m.GLDAS_SM, AMSR2_c2_AM_m.AMSR2_SM, SMAP_AM_m.SMAP_SM,2015,2016,GLDAS_AM_m.GLDAS_SM);
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/c2_AM.mat','SNR_c2_AM','R_c2_AM','fMSE_c2_AM','VAR_err_c2_AM')
toc

tic
[SNR_c2_PM, ~, R_c2_PM, fMSE_c2_PM, VAR_err_c2_PM]=TCbasedNumbers_V3(GLDAS_PM_m.GLDAS_SM, AMSR2_c2_PM_m.AMSR2_SM, SMAP_PM_m.SMAP_SM,2015,2019,GLDAS_PM_m.GLDAS_SM);
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/c2_PM.mat','SNR_c2_PM','R_c2_PM','fMSE_c2_PM','VAR_err_c2_PM')
toc
%%
  temp_lat=59.875:-0.25:-59.875;
  temp_lon=-179.875:0.25:179.875;
  [GLDAS_lon,GLDAS_lat]=meshgrid(temp_lon,temp_lat);
  t=fMSE_c2_PM.z;
  Statistic_Mapping_NDVI(GLDAS_lat, GLDAS_lon, t, min(t(:)),max(t(:)));