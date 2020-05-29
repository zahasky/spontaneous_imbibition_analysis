% particle_results_comparison.m
% Christopher Zahasky
% 11/6/2019

clear all
close all

% set a few new plot defaults
set(0,'DefaultAxesFontSize',15, 'defaultlinelinewidth', 2,...
    'DefaultAxesTitleFontWeight', 'normal')

% colorbrewer setup for better colors
% addpath to folder containing colorbrewer scripts and data
addpath('C:\Users\czahas\Dropbox\Matlab\high_res_images')
% define colormap, see demo colormap in colorbrewer folder for examples
% Example of different shades of blue light to dark. To change colormap
% change 'Blues' to other colorbrewer map. Depending on map choosen you may
% need to change 'seq' to 'div' or 'qual'. See colorbrewer documentation
% for details
% bcc = cbrewer('seq', 'RdPu', 9 , 'linear');
% pcc = cbrewer('qual', 'Set1', 9 , 'linear');
% rcc = cbrewer('seq', 'YlOrRd', 9 , 'linear');
gcc = cbrewer('seq', 'YlOrRd', 9 , 'linear');
% If you instead want dark to light use 'flipud' to flip color matrix
% bluecc = cbrewer('seq', 'Blues', 4 , 'linear');

%% Load perm data
load('ss_stream_perm_velocity_field_5psi_dp', 'perm_profile_md')
load('100_particles_SS_stream_perm_start_w_diff_2e10_fixed_bc')

% Visual check that velocity vectors appear correct
figure('position', [133 558  1680 420])
subplot(1,2,1)
gridX = [1:Grid.nx].*Grid.dx - (Grid.dx/2);
gridY = [1:Grid.ny].*Grid.dy - (Grid.dy/2);
imagesc([0 10], gridY.*100, perm_profile_md)
colormap(gray)
xlabel('Distance from inlet [cm]')
axis equal
axis tight
axis([0 10 0 max(Grid.ye).*100])
set(gca, 'Ydir', 'reverse')

hold on
plot(P.xsave(:,1).*100, P.ysave(:,1).*100, '.', 'color', gcc(8,:), 'markersize', 12)
% box on, set(gca,'linewidth',1.5)


for i = 1:P.total_particles
    % Find index of last saved particle
    ind = find(P.xsave(i, :)>0, 1, 'last');
    
    plot(P.xsave(i, 1:ind).*100, P.ysave(i, 1:ind).*100, '-', 'color', gcc(5,:), 'linewidth', 1)
    plot(P.xsave(i,ind).*100, P.ysave(i,ind).*100, '.', 'color', gcc(2,:), 'markersize', 12)
    
end

% Now plot imbibition data
load('100_imbibe_particles_frame_76_start_rand_dist_w_diff1e12_t150x')
subplot(1,2,2)
imagesc([0 10], gridY.*100, perm_profile_md)
colormap(gray)
xlabel('Distance from inlet [cm]')
axis equal
axis tight
axis([0 10 0 max(Grid.ye).*100])
set(gca, 'Ydir', 'reverse')

hold on
plot(P.xsave(:,1).*100, P.ysave(:,1).*100, '.', 'color', gcc(8,:), 'markersize', 12)
% box on, set(gca,'linewidth',1.5)


for i = 1:P.total_particles
    % Find index of last saved particle
    ind = find(P.xsave(i, :)>0, 1, 'last');
    
    plot(P.xsave(i, 1:ind).*100, P.ysave(i, 1:ind).*100, '-', 'color', gcc(5,:), 'linewidth', 1)
    plot(P.xsave(i,ind).*100, P.ysave(i,ind).*100, '.', 'color', gcc(2,:), 'markersize', 12)
    
end
            
h = colorbar('Position',...
    [0.916765872829609 0.202380952380952 0.0112103176465811 0.628142142857145]);
ylabel(h, 'Permeability [mD]');
