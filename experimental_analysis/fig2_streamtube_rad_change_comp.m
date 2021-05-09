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
addpath([current_folder(1:str_index),'Dropbox\Codes\high_res_images'])
flux_cc = flipud(cbrewer('div', 'RdYlBu', 100 , 'linear'));
scc = cbrewer('seq', 'Reds', 7 , 'linear');
% permeability
gray_map = [1 1 1; (cbrewer('seq', 'Greys', (100) , 'linear'))];

% Load PET data
load('SI_concat_PET_4D_22x22')
% Load PET single phase data
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
figure('position', [429   248   900   643])
subplot(3,2,2);
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

subplot(3,2,4);
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
ylabel(y1, 'Change in radioactivity [%]', 'fontsize', 16);

subplot(3,2,6);
h1 = imagesc(DM3);
title(['PV 0.16 to 0.57'])
axis equal
axis tight
axis off
shading flat
colormap(gca, flux_cc)
set(h1,'alphadata',~isnan(DM1))
caxis([-50 150])

subplot(3,2,[1,3, 5])
load('imbibition_rate_data')

qw_v = qw_si_dry.*(pi*2.54^2);

loglog(t_si_dry(1:end-5), qw_v, 'o', 'color', [0.4 0.4 0.4], 'MarkerSize', 5)
hold on

q01 = 0.56.*tt.^(-0.5);
% q01 = 0.607.*tt.^(-0.5);
plot(tt, q01, 'k')
load('tough2_imbibe_model_data')
plot(imbibe_model_data(:,1), imbibe_model_data(:,2), 'k--')

% mean(imbibe_model_data(:,2)./imbibe_model_data(:,1).^(-0.5))

ct_frames = [10,20,34,54];
for i = 1:length(ct_frames)
plot(t_si_dry(ct_frames(i)), qw_v(ct_frames(i)), 'o', 'MarkerEdgeColor',[0.4 0.4 0.4], 'MarkerFaceColor', scc(i+2,:) , 'MarkerSize', 8)
end

% label plots
xlabel('Time [min]')
ylabel('Volumetric rate [mL/min]')
xticks([1 10 100])
yticks([0.01 0.025 0.05 0.1 0.25 0.5 1])
% save('imbibition_rate_data', 'tt', 'q01', 'qw_si_dry', 't_si_dry')
title('Imbibition Rate')

%% Along-axis permeability
figure('position', [876   111   560   713])
subplot(2,1,1)
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

yp = colorbar;


%% Crossplot

subplot(2,1,2)
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


