% streamtube_rad_change_comp
% Christopher Zahasky
% 2/17/2017 updated 5/27/2020
% This code is for analyzing the December PET imbibition experiment
clear all
close all
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

% Set voxel size input
vox_size = [0.2329 0.2329 0.2388];

%% Calculate moments of imbibition
[M0, Xc, Sx]= streamtube_moment_calc_function(PET_matrix, vox_size(3));

% Calculate zero moment difference in each streamtube
DM1 = (M0(:,:,t2)-M0(:,:,t1))./M0(:,:,t1).*100;
DM2 = (M0(:,:,t3)-M0(:,:,t1))./M0(:,:,t1).*100;
DM3 = (M0(:,:,t4)-M0(:,:,t1))./M0(:,:,t1).*100;

% plot measured change
figure('position', [214   313   787   565])
subplot(3,2,1);
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

subplot(3,2,3);
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

subplot(3,2,5);
h1 = imagesc(DM3);
title(['PV 0.16 to 0.57'])
axis equal
axis tight
axis off
shading flat
colormap(gca, flux_cc)
set(h1,'alphadata',~isnan(DM1))
caxis([-50 150])

%% Along-axis permeability
subplot(3,2,2)
mean_perm = nanmean(Kmd_mat(2:end-1, 2:end-1,:),3);
% standard error
std_err = nanstd(Kmd_mat(2:end-1, 2:end-1,:),0, 3)./sqrt(5);
h2 = imagesc(mean_perm);
title(['Permeability [mD]'])
axis equal
axis tight
axis off
shading flat
caxis([10  32])

set(h2,'alphadata',~isnan(DM1))
colormap(gca, gray_map)
hold on
% draw box around area of interest
plot([10.5 10.5], [0.5 20.5], 'r', 'linewidth', 1) 
plot([9.5 9.5], [0.5 20.5], 'r', 'linewidth', 1) 
plot([10.5 9.5], [0.5 0.5], 'r', 'linewidth', 1) 
plot([10.5 9.5], [20.5 20.5], 'r', 'linewidth', 1) 

% pos1 = get(hg, 'Position')
% new_pos1 = pos1 +[-0.04 -0.09 0 0.07]
% set(hg, 'Position',new_pos1 )

%% Crossplot
subplot(3,2,[4,6])
hold on
scatter(mean_perm(:), DM3(:), std_err(:).*80, 'k')
% Format plot
axis([15 32 -50 155])
box on
title('Cross-plot')
xlabel('Permeability [mD]')
ylabel('Change in radioactivity [%]')
set(gca, 'Color', 'none');
set(gca,'linewidth',1.1)
