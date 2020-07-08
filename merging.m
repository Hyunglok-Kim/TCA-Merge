clear all; clc
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/SAT_data_related_CODE')
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/TCA')
addpath('/sfs/qumulo/qproject/hydrosense/matlab/libs/mapping_code/')

load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/X_AM.mat')
load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/c1_AM.mat')
load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/c2_AM.mat')
% x=X, y=c1, z=c2
% Ylimaz et al.
% var_x=VAR_err_X_AM.y;
% var_y=VAR_err_c1_AM.y;
% var_z=VAR_err_c2_AM.y;

% Gruber et al
var_x=SNR_X_AM.y;
var_y=SNR_c1_AM.y;
var_z=SNR_c2_AM.y;
%
var_check=sum(~isnan(cat(3, var_x, var_y, var_z)),3);
var_xy=var_x.*var_y;
var_xz=var_x.*var_z;
var_yz=var_y.*var_z;

wx=var_yz./(var_xy+var_xz+var_yz);
wy=var_xz./(var_xy+var_xz+var_yz);
wz=var_xy./(var_xy+var_xz+var_yz);
var_m=wx.^2.*var_x+wy.^2.*var_y+wz.^2.*var_z;

% valid_var=find(var_check==2);
% for i=1:numel(valid_var)
%
%     t_x=var_x(valid_var(i));
%     t_y=var_y(valid_var(i));
%     t_z=var_z(valid_var(i));
%
%     t_sum=nansum([t_x, t_y, t_z]);
%
%     t_w_x=~isnan(t_x)*nansum([t_y,t_z])/t_sum;
%     t_w_y=~isnan(t_y)*nansum([t_x,t_z])/t_sum;
%     t_w_z=~isnan(t_z)*nansum([t_x,t_y])/t_sum;
%
%     wx(valid_var(i))=t_w_x;
%     wy(valid_var(i))=t_w_y;
%     wz(valid_var(i))=t_w_z;
%
%     if t_w_x==0
%        t_var_m=(t_w_y^2)*var_y(valid_var(i))+((1-t_w_y)^2)*var_z(valid_var(i));
%     elseif t_w_y==0
%        t_var_m=(t_w_x^2)*var_x(valid_var(i))+((1-t_w_x)^2)*var_z(valid_var(i));
%     elseif t_w_z==0
%        t_var_m=(t_w_x^2)*var_x(valid_var(i))+((1-t_w_x)^2)*var_y(valid_var(i));
%     end
%     %[var_x(valid_var(i)),var_y(valid_var(i)),var_z(valid_var(i)),t_var_m]
%     var_m(valid_var(i))=t_var_m;
% end
% %wx(var_check<2)=nan;wy(var_check<2)=nan;wz(var_check<2)=nan;
%
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
AMSR2_c1_AM_m=matfile(ifp_list{3});
AMSR2_c2_AM_m=matfile(ifp_list{5});
SMAP_AM_m=matfile(ifp_list{7});
GLDAS_AM_m=matfile(ifp_list{9});

X_SM=AMSR2_X_AM_m.AMSR2_SM;%(:,:,100:300);
c1_SM=AMSR2_c1_AM_m.AMSR2_SM;%(:,:,100:300);
c2_SM=AMSR2_c2_AM_m.AMSR2_SM;%(:,:,100:300);
GLDAS_SM=GLDAS_AM_m.GLDAS_SM;%(:,:,100:300);
% 
% disp('cdf 1...')
% [i1,i2]=find(var_check>=1);
% for i=1:numel(i1)
%     %plot(squeeze(GLDAS_SM(i1(i), i2(i),:)))
%     %hold on
%     %plot(squeeze(X_SM(i1(i), i2(i),:)))
%     %hold on
%     %plot(CDF_match([squeeze(GLDAS_SM(i1(i), i2(i),:)), squeeze(X_SM(i1(i), i2(i),:))]),'k')
%     
%     X_SM(i1(i), i2(i),:)=CDF_match([squeeze(GLDAS_SM(i1(i), i2(i),:)), squeeze(X_SM(i1(i), i2(i),:))]);
%     c1_SM(i1(i), i2(i),:)=CDF_match([squeeze(GLDAS_SM(i1(i), i2(i),:)), squeeze(c1_SM(i1(i), i2(i),:))]);
%     c2_SM(i1(i), i2(i),:)=CDF_match([squeeze(GLDAS_SM(i1(i), i2(i),:)), squeeze(c2_SM(i1(i), i2(i),:))]);
%     
% end
% 
% disp('cdf 2...')
% [i1,i2]=find(var_check<1);
% for i=1:numel(i1)
%     X_SM(i1(i), i2(i),:)=nan;
%     c1_SM(i1(i), i2(i),:)=nan;
%     c2_SM(i1(i), i2(i),:)=nan;
% end

%save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_CDF.mat','X_SM','c1_SM','c2_SM','-v7.3')
%load('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_CDF.mat')
%%
% disp('combining data...')
% SM_m=nan(size(X_SM));
% FLAG=SM_m;
% SM_avg=SM_m;
% % when three w factors are available
% disp('three w case')
% [i1,i2]=find(var_check==3);
% for i=1:numel(i1)    
%     t_SM_m=[];t_FLAG=[]; t_SM_avg=[];
%     for k=1:size(X_SM,3)
%         
%         t_x=X_SM(i1(i), i2(i),k);t_y=c1_SM(i1(i), i2(i),k);t_z=c2_SM(i1(i), i2(i),k);
%         t_xyz=[t_x, t_y, t_z];
%         tt_xyz=~isnan(t_xyz);
%         t_sm_avg=nanmean(t_xyz);
%         
%         switch sum(tt_xyz)
%             case 0
%                 t_sm_m=nan;
%                 t_flag=12;
%             case 1
%                 t_sm_m=t_sm_avg;
%                 tt_flag=[9 10 11];
%                 t_flag=tt_flag(tt_xyz);
%             case 2
%                 t_var_x=var_x(i1(i), i2(i));
%                 t_var_y=var_y(i1(i), i2(i));
%                 t_var_z=var_z(i1(i), i2(i));
%                 t_var_xyz=[t_var_x, t_var_y, t_var_z];
%                 t_var_xyz=t_var_xyz.*tt_xyz;
%                 
%                 t_sum=nansum(t_var_xyz);
%                 t_w_x=tt_xyz(1)*nansum([t_var_xyz(2),t_var_xyz(3)])/t_sum;
%                 t_w_y=tt_xyz(2)*nansum([t_var_xyz(1),t_var_xyz(3)])/t_sum;
%                 t_w_z=tt_xyz(3)*nansum([t_var_xyz(1),t_var_xyz(2)])/t_sum;
%                 t_w_xyz=[t_w_x, t_w_y, t_w_z];
%                 
%                 t_var_m=(t_w_xyz.^2)*t_var_xyz';
%                 t_sm_m=nansum(t_w_xyz.*t_xyz);
%                 tt_flag=[5 4 3];
%                 t_flag=tt_flag(~tt_xyz);
%             case 3
%                 t_var_x=var_x(i1(i), i2(i));
%                 t_var_y=var_y(i1(i), i2(i));
%                 t_var_z=var_z(i1(i), i2(i));
%                 t_var_xyz=[t_var_x, t_var_y, t_var_z];
%                 t_w_xyz=[wx(i1(i), i2(i)), wy(i1(i), i2(i)), wz(i1(i), i2(i))];
%                 t_var_m=(t_w_xyz.^2)*t_var_xyz';
%                 t_sm_m=nansum(t_w_xyz.*t_xyz);
%                 t_flag=1;
%         end
%         t_SM_m(k)=t_sm_m;
%         t_FLAG(k)=t_flag;
%         t_SM_avg(k)=t_sm_avg;
%     end
%     SM_m(i1(i), i2(i),:)=t_SM_m;
%     FLAG(i1(i), i2(i),:)=t_FLAG;
%     SM_avg(i1(i), i2(i),:)=t_SM_avg;
%     
%     % check flags
%      % [squeeze(X_SM(i1(i), i2(i),:)),squeeze(c1_SM(i1(i), i2(i),:)),squeeze(c2_SM(i1(i), i2(i),:)),t_SM_m',t_FLAG']
% %
% %      plot(t_SM_m,'o-r')
% %      hold on
% %      plot(t_SM_avg,'g')
% %      hold on
% %      plot(squeeze(GLDAS_SM(i1(i), i2(i),:)),'b')
% end
% 
% %%
% % when two w factors are available
% disp('two w case')
% [i1,i2]=find(var_check==2);
% for i=1:numel(i1)
%     t_SM_m=[];t_FLAG=[]; t_SM_avg=[];
%     for k=1:size(X_SM,3)
%         
%         t_x=X_SM(i1(i), i2(i),k);t_y=c1_SM(i1(i), i2(i),k);t_z=c2_SM(i1(i), i2(i),k);
%         t_xyz=[t_x, t_y, t_z];
%         tt_xyz=~isnan(t_xyz);
%         t_sm_avg=nanmean(t_xyz);
%         
%         switch sum(tt_xyz)
%             case 0
%                 t_sm_m=nan;
%                 t_flag=12;
%             case 1
%                 t_sm_m=t_sm_avg;
%                 tt_flag=[9 10 11];
%                 t_flag=tt_flag(tt_xyz);
%             case 2
%                 t_var_x=var_x(i1(i), i2(i));
%                 t_var_y=var_y(i1(i), i2(i));
%                 t_var_z=var_z(i1(i), i2(i));
%                 t_var_xyz=[t_var_x, t_var_y, t_var_z];
%                 t_var_xyz=t_var_xyz.*tt_xyz;
%                 
%                 t_sum=nansum(t_var_xyz);
%                 t_w_x=tt_xyz(1)*nansum([t_var_xyz(2),t_var_xyz(3)])/t_sum;
%                 t_w_y=tt_xyz(2)*nansum([t_var_xyz(1),t_var_xyz(3)])/t_sum;
%                 t_w_z=tt_xyz(3)*nansum([t_var_xyz(1),t_var_xyz(2)])/t_sum;
%                 t_w_xyz=[t_w_x, t_w_y, t_w_z];
%                 
%                 t_var_m=(t_w_xyz.^2)*t_var_xyz';
%                 t_sm_m=nansum(t_w_xyz.*t_xyz);
%                 tt_flag=[5 4 3];
%                 t_flag=tt_flag(~tt_xyz);
%             case 3
%                 t_sm_m=t_sm_avg;
%                 t_flag=2;
%         end
%         t_SM_m(k)=t_sm_m;
%         t_FLAG(k)=t_flag;
%         t_SM_avg(k)=t_sm_avg;
%     end
%     SM_m(i1(i), i2(i),:)=t_SM_m;
%     FLAG(i1(i), i2(i),:)=t_FLAG;
%     SM_avg(i1(i), i2(i),:)=t_SM_avg;
%     % check flags
%     % [squeeze(X_SM(i1(i), i2(i),:)),squeeze(c1_SM(i1(i), i2(i),:)),squeeze(c2_SM(i1(i), i2(i),:)),t_SM_m',t_FLAG']
%     % plot(t_SM_m,'o-r')
%     % hold on
%     % plot(t_SM_avg,'g')
%     % hold on
%     % plot(squeeze(GLDAS_SM(i1(i), i2(i),:)),'b')
% end
% %% when just one w factor is available
% disp('one w case')
% [i1,i2]=find(var_check==1);
% for i=1:1%numel(i1)
%     [i1(i), i2(i)]
%     t_SM_m=[];t_FLAG=[]; t_SM_avg=[];
%     for k=1:size(X_SM,3)
%         
%         t_x=X_SM(i1(i), i2(i),k);t_y=c1_SM(i1(i), i2(i),k);t_z=c2_SM(i1(i), i2(i),k);
%         t_xyz=[t_x, t_y, t_z];
%         tt_xyz=~isnan(t_xyz);
%         t_sm_avg=nanmean(t_xyz);
%         
%         switch sum(tt_xyz)
%             case 0
%                 t_sm_m=nan;
%                 t_flag=12;
%             case 1
%                 t_sm_m=t_sm_avg;
%                 tt_flag=[9 10 11];
%                 t_flag=tt_flag(tt_xyz);
%             case 2
%                 t_sm_m=t_sm_avg;
%                 tt_flag=[8 7 6];
%                 t_flag=tt_flag(~tt_xyz);
%             case 3
%                 t_sm_m=t_sm_avg;
%                 t_flag=2;
%         end
%         t_SM_m(k)=t_sm_m;
%         t_FLAG(k)=t_flag;
%         t_SM_avg(k)=t_sm_avg;
%     end
%     SM_m(i1(i), i2(i),:)=t_SM_m;
%     FLAG(i1(i), i2(i),:)=t_FLAG;
%     SM_avg(i1(i), i2(i),:)=t_SM_avg;
%     % check flags
%     % t_check= [squeeze(X_SM(i1(i), i2(i),:)),squeeze(c1_SM(i1(i), i2(i),:)),squeeze(c2_SM(i1(i), i2(i),:)),t_SM_m',t_FLAG']
%     % [var_x(i1(i), i2(i)), var_y(i1(i), i2(i)), var_z(i1(i), i2(i))]
%     % [i1(i), i2(i)]    
%     %
%     % plot(t_SM_m,'o-r')
%     % hold on
%     % plot(t_SM_avg,'g')
%     % hold on
%     % plot(squeeze(GLDAS_SM(i1(i), i2(i),:)),'b')
% end
% %clearvars c1_SM c2_SM X_SM
% 
% %disp('saving merged data')
% %save('/project/hydrosense/matlab/mat/AMSR2/Merged/AMSR2_025.mat','SM_m','SM_avg','FLAG','-v7.3')

%% TC for both M and avg data
SMAP_SM=SMAP_AM_m.SMAP_SM;
load('/project/hydrosense/matlab/mat/AMSR2/Merged/AMSR2_025.mat')
clearvars FLAG
check_SM=SM_m + SM_avg + X_SM + c1_SM + c2_SM;
nan_val=isnan(check_SM);
clearvars check_SM
SM_m(nan_val)=nan;
SM_avg(nan_val)=nan;
X_SM(nan_val)=nan;
c1_SM(nan_val)=nan;
c2_SM(nan_val)=nan;

disp('TC calculation X')
tic
[SNR_X_AM, ~, R_X_AM, fMSE_X_AM, VAR_err_X_AM]=TCbasedNumbers_V3(GLDAS_SM, X_SM, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving X TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_X_merge_compare.mat','SNR_X_AM','R_X_AM','fMSE_X_AM','VAR_err_X_AM')
toc

disp('TC calculation c1')
tic
[SNR_c1_AM, ~, R_c1_AM, fMSE_c1_AM, VAR_err_c1_AM]=TCbasedNumbers_V3(GLDAS_SM, c1_SM, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving c1 TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_c1_merge_compare.mat','SNR_c1_AM','R_c1_AM','fMSE_c1_AM','VAR_err_c1_AM')
toc

disp('TC calculation c2')
tic
[SNR_c2_AM, ~, R_c2_AM, fMSE_c2_AM, VAR_err_c2_AM]=TCbasedNumbers_V3(GLDAS_SM, c2_SM, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving c2 TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_c2_merge_compare.mat','SNR_c2_AM','R_c2_AM','fMSE_c2_AM','VAR_err_c2_AM')
toc

disp('TC calculation M')
tic
[SNR_M_AM, ~, R_M_AM, fMSE_M_AM, VAR_err_M_AM]=TCbasedNumbers_V3(GLDAS_SM, SM_m, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving merged TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_M.mat','SNR_M_AM','R_M_AM','fMSE_M_AM','VAR_err_M_AM')
toc

disp('TC calculation avg')
tic
[SNR_avg_AM, ~, R_avg_AM, fMSE_avg_AM, VAR_err_avg_AM]=TCbasedNumbers_V3(GLDAS_SM, SM_avg, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving avg TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_avg.mat','SNR_avg_AM','R_avg_AM','fMSE_avg_AM','VAR_err_avg_AM')
toc

%% TC for both M and avg data (flag applied)
load('/project/hydrosense/matlab/mat/AMSR2/Merged/AMSR2_025.mat')

disp(nanmean(SM_m(:)))

SM_m(FLAG==2 | FLAG==6 | FLAG==7 | FLAG==8)=nan;
check_SM=SM_m + SM_avg + X_SM + c1_SM + c2_SM;
nan_val=isnan(check_SM);
clearvars check_SM
SM_m(nan_val)=nan;
SM_avg(nan_val)=nan;
X_SM(nan_val)=nan;
c1_SM(nan_val)=nan;
c2_SM(nan_val)=nan;

disp(nanmean(SM_m(:)))

disp('TC calculation X')
tic
[SNR_X_AM_flag, ~, R_X_AM_flag, fMSE_X_AM_flag, VAR_err_X_AM_flag]=TCbasedNumbers_V3(GLDAS_SM, X_SM, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving X TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_X_merge_compare_flag.mat','SNR_X_AM_flag','R_X_AM_flag','fMSE_X_AM_flag','VAR_err_X_AM_flag')
toc

disp('TC calculation c1')
tic
[SNR_c1_AM_flag, ~, R_c1_AM_flag, fMSE_c1_AM_flag, VAR_err_c1_AM_flag]=TCbasedNumbers_V3(GLDAS_SM, c1_SM, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving c1 TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_c1_merge_compare_flag.mat','SNR_c1_AM_flag','R_c1_AM_flag','fMSE_c1_AM_flag','VAR_err_c1_AM_flag')
toc

disp('TC calculation c2')
tic
[SNR_c2_AM_flag, ~, R_c2_AM_flag, fMSE_c2_AM_flag, VAR_err_c2_AM_flag]=TCbasedNumbers_V3(GLDAS_SM, c2_SM, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving c2 TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_c2_merge_compare_flag.mat','SNR_c2_AM_flag','R_c2_AM_flag','fMSE_c2_AM_flag','VAR_err_c2_AM_flag')
toc

disp('TC calculation M')
tic
[SNR_M_AM_flag, ~, R_M_AM_flag, fMSE_M_AM_flag, VAR_err_M_AM_flag]=TCbasedNumbers_V3(GLDAS_SM, SM_m, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving merged TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_M_flag.mat','SNR_M_AM_flag','R_M_AM_flag','fMSE_M_AM_flag','VAR_err_M_AM_flag')
toc

disp('TC calculation avg')
tic
[SNR_avg_AM_flag, ~, R_avg_AM_flag, fMSE_avg_AM_flag, VAR_err_avg_AM_flag]=TCbasedNumbers_V3(GLDAS_SM, SM_avg, SMAP_SM, 2015, 2019, GLDAS_SM);
disp('saving avg TC ressults')
save('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_avg_flag.mat','SNR_avg_AM_flag','R_avg_AM_flag','fMSE_avg_AM_flag','VAR_err_avg_AM_flag')
toc
%% 3/81
% i1=3; i2=81;
% 
% %tt_sm=squeeze(SM_m(i1,i2,:));
% 
% %t_mat=matfile('/project/hydrosense/matlab/mat/AMSR2/Merged/AMSR2_025.mat');
% t_sm=squeeze(t_mat.SM_m(i1,i2,100:300));
% t_flag=squeeze(t_mat.FLAG(i1,i2,100:300));
% 
% %sm_check=[tt_sm, t_sm];
% 
% 
% SM_cdf_m=matfile('/project/hydrosense/matlab/mat/TCresults_025_CDF_AMSR2/AMSR2_CDF.mat');
% 
% t_X_SM=squeeze(SM_cdf_m.X_SM(i1,i2,100:300));
% t_c1_SM=squeeze(SM_cdf_m.c1_SM(i1,i2,100:300));
% t_c2_SM=squeeze(SM_cdf_m.c2_SM(i1,i2,100:300));
% 
% t_check=[t_X_SM, t_c1_SM, t_c2_SM, t_sm, t_flag]
% [var_x(i1, i2), var_y(i1, i2), var_z(i1, i2)]
%%

% % %%
% % temp_lat=59.875:-0.25:-59.875;
% % temp_lon=-179.875:0.25:179.875;
% % [GLDAS_lon,GLDAS_lat]=meshgrid(temp_lon,temp_lat);
% % t=nanmean(SM_avg,3) - nanmean(SM_m,3);
% % Statistic_Mapping_NDVI(GLDAS_lat, GLDAS_lon, t, min(t(:)),max(t(:)));
%%
system('ls')

