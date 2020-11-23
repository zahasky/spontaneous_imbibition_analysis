% streamtube_rad_change_comp
% Christopher Zahasky
% 2/17/2017 updated 5/27/2020
% This code is for analyzing the December PET imbibition experiment
clear all
% close all
set(0,'DefaultAxesFontSize',16, 'defaultlinelinewidth', 1.1,...
    'DefaultAxesTitleFontWeight', 'normal')

% adjust paths depending on computer
current_folder = pwd;
str_index = strfind(pwd, '\Dropbox');

% Path to colorbrewer
addpath([current_folder(1:str_index),'Dropbox\Matlab\high_res_images'])
flux_cc = flipud(cbrewer('div', 'RdYlBu', 100 , 'linear'));
scc = cbrewer('seq', 'Purples', 7 , 'linear');
% permeability
gray_map = [1 1 1; (cbrewer('seq', 'Greys', (100) , 'linear'))];

% Load PET data
load('SI_concat_PET_4D_22x22')
load('BSS_c1_2ml_2_3mm_vox')

% Add path to perm maps and load them
load('2D_perm_data_input.mat')


% Crop matrix
PET_matrix = SI_concat_PET_4D(2:end-1, 2:end-1, 1:45,:);
% saturated
PET_matrix_sat = PET_4D_coarse(2:end-1, 2:end-1, :,:);
pet_size = size(PET_matrix_sat);

% times between which change in radioactivity will be measured
t1 = 10;
t2 = 19;
t3 = 36;
t4 = 74;

vox_size = [0.2329 0.2329 0.2388];


%% Calculate moments of imbibition
[M0, Xc, Sx]= streamtube_moment_calc_function(PET_matrix, vox_size(3));

% Calculate zero moment difference in each streamtube
DM1 = (M0(:,:,t2)-M0(:,:,t1))./M0(:,:,t1).*100;
DM2 = (M0(:,:,t3)-M0(:,:,t1))./M0(:,:,t1).*100;
DM3 = (M0(:,:,t4)-M0(:,:,t1))./M0(:,:,t1).*100;
% DM(DM==0)=nan;

% plot measured change
figure('position', [214   313   787   565])
% figure('position', [203         558        1715         420])
s1 = subplot(3,2,1);
h1 = imagesc(DM1);
title(['PV 0.16 to 0.25'])
axis equal
axis tight
axis off
shading flat
% colorbar('fontsize', 15)
colormap(gca, flux_cc)
set(h1,'alphadata',~isnan(DM1))
caxis([-50 150])

s1 = subplot(3,2,3);
h1 = imagesc(DM2);
title(['PV 0.16 to 0.40'])
axis equal
axis tight
axis off
shading flat
y1 = colorbar;
colormap(gca, flux_cc)
set(h1,'alphadata',~isnan(DM1))
caxis([-50 150])
ylabel(y1, 'Change in radioactivity [%]');

s1 = subplot(3,2,5);
h1 = imagesc(DM3);
title(['PV 0.16 to 0.57'])
axis equal
axis tight
axis off
shading flat
% colorbar('fontsize', 16)
colormap(gca, flux_cc)
set(h1,'alphadata',~isnan(DM1))
caxis([-50 150])


% plot flux-calculate change
hg = subplot(3,2,2)
mean_perm = nanmean(Kmd_mat(2:end-1, 2:end-1,:),3);
% standard error
std_err = nanstd(Kmd_mat(2:end-1, 2:end-1,:),0, 3)./sqrt(5);
h2 = imagesc(mean_perm);
title(['Permeability [mD]'])
axis equal
axis tight
axis off
shading flat
% colorbar('fontsize', 15);
caxis([10  32])

set(h2,'alphadata',~isnan(DM1))
colormap(gca, gray_map)
hold on
% draw box around area of interest
plot([10.5 10.5], [0.5 20.5], 'r', 'linewidth', 1) 
plot([9.5 9.5], [0.5 20.5], 'r', 'linewidth', 1) 
plot([10.5 9.5], [0.5 0.5], 'r', 'linewidth', 1) 
plot([10.5 9.5], [20.5 20.5], 'r', 'linewidth', 1) 

pos1 = get(hg, 'Position')
new_pos1 = pos1 +[-0.04 -0.09 0 0.07]
set(hg, 'Position',new_pos1 )

%% Crossplot
h1 = subplot(3,2,[4,6])
% scatter(mean_perm(:), DM1(:), std_err(:).*80, 'MarkerEdgeColor', scc(3,:))
hold on
% scatter(mean_perm(:), DM2(:), std_err(:).*80, 'MarkerEdgeColor', scc(5,:))
% scatter(mean_perm(:), DM3(:), std_err(:).*80, 'MarkerEdgeColor', scc(7,:))
scatter(mean_perm(:), DM3(:), std_err(:).*80, 'k')
% % Linear regression
% xdata = mean_perm(:);
% ydata = DM(:);
% % remove potential nans
% xdata(isnan(ydata))=[];
% ydata(isnan(ydata))=[];
% 
% X = [ones(length(xdata),1) xdata];
% % b(2) is the slope and b(1) is the intercept
% b = X\ydata;
% % sythenthic vector of x-values for plotting regression results
% xreg = linspace(min(xdata), max(xdata), 100);
% Plot when regression intercept is not equal to o
% plot(xreg, xreg.*b(2)+b(1), ':r', 'linewidth', 2)
% Plot when regression intercept is equal to 0
%     plot(xreg, xreg.*b(1), ':r', 'linewidth', 2)

% fitted linear regression model calculation
% yCalc1 = X*b;

% nonlinear regression
% p = polyfit(xdata,ydata,2);
% yreg = polyval(p, xreg);
% plot(xreg, yreg, ':r', 'linewidth', 2)
% yCalc1 = polyval(p, xdata);
% 
% % r^2 coefficient of determination (see wiki page for details)
% % fit_data = 
% residuals = (ydata - yCalc1);
% % The sum of squares of residuals, also called the residual sum of squares
% ssres = sum(residuals.^2);
% % The total sum of squares (proportional to the variance of the data):
% sstot = sum((ydata- mean(ydata)).^2);
% % coefficient of determination
% R2 = 1- (ssres/sstot)


axis([15 32 -50 155])
box on
title('Cross-plot')
xlabel('Permeability [mD]')
ylabel('Change in radioactivity [%]')
set(gca, 'Color', 'none');
set(gca,'linewidth',1.1)

% figure
% hold on
% sumM0 = squeeze(nansum(nansum(nansum(PET_matrix_sat))));
% plot([1:pet_size(4)], sumM0./max(sumM0), '*b')
% % plot([0, 300], [mean_peak_activity, mean_peak_activity], 'r')
% % Finally plot indication of times between when rad changes are calculated
% % plot([T(t1), T(t1)], [0 1], '--k', 'linewidth', 1)
% % plot([T(t2), T(t2)], [0 1], '--k', 'linewidth', 1)
% xlabel('Time (minutes)', 'fontsize', 14)
% ylabel('Total activity in core (-)', 'fontsize', 14)
% box on