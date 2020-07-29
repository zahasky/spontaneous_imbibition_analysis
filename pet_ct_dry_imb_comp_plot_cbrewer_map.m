% pet_ct_dry_imb_comp_plot_cbrewer_map.m
% Christopher Zahasky
% 5/27/2020
clear all
close all
set(0,'DefaultAxesFontSize',14, 'defaultlinelinewidth', 0.5,...
    'DefaultAxesTitleFontWeight', 'normal')
% adjust paths depending on computer
current_folder = pwd;
str_index = strfind(pwd, '\Dropbox');

% Path to colorbrewer maps and PATCH_3Darray
addpath([current_folder(1:str_index),'Dropbox\Matlab\high_res_images'])

% plot PET data first
load('SI_concat_PET_4D_22x22')
comp_PET_4D_coarse = SI_concat_PET_4D(:,:,1:41,:)./0.32;
clim = [0 1];
% PET voxel size
vox_size = [0.2329 0.2329 0.2388];
% Define pet grid
s = size(comp_PET_4D_coarse);
gridX = ([1:s(1)].*vox_size(1) - vox_size(1)/2);
gridY = ([1:s(2)].*vox_size(2) - vox_size(2)/2);
gridZ = ([1:s(3)].*vox_size(3) - vox_size(3)/2);

% Setup custom colormaps
% Set create custom PET colormap
ind_trans = 62;
reds = cbrewer('div', 'RdGy', 120 , 'linear');
greys = ones(1,3).*reds(ind_trans,:);
white_jet = flipud([reds(1:ind_trans, :); greys]);

% Create custom saturation colormap
blues = cbrewer('seq', 'Blues', 70 , 'linear');
grey_map = flipud(cbrewer('seq', 'Greys', 30 , 'linear'));
sat_blue = [greys(2:end, :); blues(2:end,:)];

% Create figure
shp = figure;
set(shp, 'Position', [750    64   804   930])
% Plot PET data first
n=2;
subplot(4,2,n)
title('PET')

for i=[10,19,36,74]
    subplot(4,2,n)
    slice_plane = squeeze(comp_PET_4D_coarse(:,:,:,i));

    % Crop to half-core (NANs are ignored in PATCH_3Darray plot function)
    slice_plane(1:end,11:end,:) = nan;
    slice_plane = flip(slice_plane);
    slice_plane = flip(slice_plane,2);
    slice_plane = permute(slice_plane,[3 2 1]);
    
    PATCH_3Darray(slice_plane, gridZ, gridY, gridX, white_jet, clim, 'col')
    
    % format plot
    axis([0 max(gridZ) 0 max(gridY) 0 max(gridX)])
    grid on
    axis equal
    axis tight
    xticks([0:2:10])
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
    view(-15,32)
    set(gca,'color','none')
   
    n=n+2;
end
subplot(4,2,7)
xlabel('Distance from inlet [cm]')
subplot(4,2,3)
y1 = colorbar('Position',...
    [0.453865336658354 0.144086021505376 0.0187032418952618 0.772043010752688]);
ylabel(y1, 'Water Saturation', 'fontsize', 16);


% plot CT data 
% Load primary imbibition
load SI_CT_coarse_sat.mat
SAT = SAT(2:end-1, 2:end-1,1:end-2,:);
s= size(SAT);
% Define CT grid
gridX = ([1:s(1)].*2.5 - 2.5/2)./10;
gridY = ([1:s(2)].*2.5 - 2.5/2)./10;
gridZ = ([1:s(3)].*2.5 - 2.5/2)./10;

satlim = [0.0 0.75];
n = 1;
subplot(4,2,n)
title('X-ray CT')
for i=[10,20,34,54]
   
    mean_sat = nanmean(nanmean(nanmean(SAT(:,:,:,i))));
    
    subplot(4,2,n)
    slice_plane = squeeze(SAT(:,:,:,i));
    slice_plane(1:end,10:end,:) = nan;
    slice_plane = flip(slice_plane);
    slice_plane = flip(slice_plane,2);
    slice_plane = permute(slice_plane,[3 2 1]);
    
    PATCH_3Darray(slice_plane, gridZ, gridY, gridX, sat_blue, satlim, 'col')
    %     title(['Imbibition time ', num2str(round(tt(nn))), ' minutes'], 'FontWeight', 'Normal')
    % title('Water saturation')
    axis([0 max(gridZ) 0 max(gridY) 0 max(gridX)])
    grid on
    axis equal
    axis tight
    xticks([0:2:10])
    set(gca,'YTickLabel',[]);
    set(gca,'ZTickLabel',[]);
    view(-15,32)
    set(gca,'color','none')
    
    subplot(4,2,n)
    zlabel(sprintf('PV: %.2f', mean_sat))

    n=n+2;
end

subplot(4,2,8)
xlabel('Distance from inlet [cm]')
y2 = colorbar('Position',...
    [0.897720686970357 0.144086021505376 0.0187381658974735 0.774193548387097]);
ylabel(y2, 'Normalized radioactivity concentration', 'fontsize', 16);
colormap(gca, white_jet)
caxis([0 1])
