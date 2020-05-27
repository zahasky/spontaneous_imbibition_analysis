% imibibe_flux_solution.m
% Christopher Zahasky
% 2/17/2017
% This code is for analyzing the December PET imbibition experiment
clear all
% close all
set(0,'DefaultAxesFontSize',15, 'defaultlinelinewidth', 1,...
    'DefaultAxesTitleFontWeight', 'normal')

% adjust paths depending on computer
current_folder = pwd;
str_index = strfind(pwd, '\Dropbox');

% Path to colorbrewer
addpath([current_folder(1:str_index),'Dropbox\Matlab\high_res_images'])
flux_cc = flipud(cbrewer('div', 'RdYlBu', 100 , 'linear'));
% perm_cc = flipud(cbrewer('div', 'RdYlBu', 100 , 'linear'));


addpath([current_folder(1:str_index),'Dropbox\Matlab\moment_analysis'])
addpath('C:\Users\zahas\Dropbox\Research\Experiment stuff\Data\Stanford_data\BSS_c1\PET_SI_compiled_data')
% cd([current_folder(1:str_index),'Dropbox\Research\Experiment stuff\Data\BSS_c1\PET_SI_compiled_data'])
load('SI_concat_PET_4D_22x22')

% Add path to perm maps and load them
% addpath([current_folder(1:str_index),'Dropbox\Research\Experiment stuff\Data\BSS_c1\june_17_pet\6_12_single_phase'])
load('2D_perm_data_input.mat')

pet_size = size(SI_concat_PET_4D);
% PET_matrix = SI_concat_PET_4D;
% Crop matrix
PET_matrix = SI_concat_PET_4D(2:end-1, 2:end-1, 1:45,:);
% PET_matrix = flip(PET_matrix,3);

sumM0 = squeeze(nansum(nansum(nansum(PET_matrix))));
% Normalize total activity
PET_matrix = PET_matrix./max(sumM0);

% times between which change in radioactivity will be measured
t1 = 10;
t2 = 74;
% t1 = 70;
% t2 = 118;

vox_size = [0.2329 0.2329 0.2388];


%% change in activity between defined timesteps

% crop edges
% PET_matrix = PET_matrix(2:end-1, 2:end-1, 1:45,:);

% sumM0 = squeeze(nansum(nansum(nansum(PET_matrix))));
% mean_peak_activity = mean(sumM0(25:30));

figure(3)
hold on
% norm_act = max(total_activity1);
plot(T, sumM0./max(sumM0), '*b')
% plot([0, 300], [mean_peak_activity, mean_peak_activity], 'r')
xlabel('Time (minutes)', 'fontsize', 14)
ylabel('Total activity in core (-)', 'fontsize', 14)
box on

%% scale activity levels so that total activity is constant during entire
% imbibition
% for i=30:97
%     % calculate scale factor
%     sf = mean_peak_activity./squeeze(nansum(nansum(nansum(PET_matrix(:,:,:,i)))));
%     PET_matrix(:,:,:,i) = PET_matrix(:,:,:,i).*sf;
%     
% end
% PET_matrix(:,:,:,i+1:end) = PET_matrix(:,:,:,i+1:end).*sf;

% Scanner-derived concentration
% Cm = mean_peak_activity/vol_injected;

[M0, Xc, Sx]= streamtube_moment_calc_function(PET_matrix, vox_size(3));

% Determine system dimenstions
number_cells = length(find(M0(:,:,1)>0));
PET_dim = size(PET_matrix);

% sumM03 = squeeze(nansum(nansum(nansum(PET_matrix))));
% % Plot new scaled activity sums
% plot(T, sumM03, '*g')

% Finally plot indication of times between when rad changes are calculated
plot([T(t1), T(t1)], [0 1], '--k', 'linewidth', 1)
plot([T(t2), T(t2)], [0 1], '--k', 'linewidth', 1)



%% Analyze how rad in each streamtube changes with time

dm_dt = zeros(number_cells, (t2-t1));
n=1;
for i=t1:1:t2
    Diff_m = (M0(:,:,i)-M0(:,:,i-1))';
    diff_mv = Diff_m(:);
    diff_mv(diff_mv==0)=[];
    if n >1
        dm_dt(:,n) = dm_dt(:,n-1)+diff_mv;
    else
        dm_dt(:,n) = diff_mv;
    end
    n=n+1;
end

figure
set(gca, 'ColorOrder', jet(233), 'NextPlot', 'replacechildren');
plot(T(t1:t2), dm_dt./max(max(dm_dt)))
xlabel('Time (minutes)', 'fontsize', 14)
ylabel('Total change in activity (-)', 'fontsize', 14)
box on

%% Invert data to solve for flux (this is done with least norm solution
% for underdetermined matrices

% Preallocate A matrix
A_columns = ((PET_dim(1)+1)* PET_dim(2)) + (PET_dim(1)* (PET_dim(2) +1));
A = zeros(number_cells, A_columns);
% Keep track of location of Qx (left right fluxes). May not be neccessary
Qx_ind = zeros(PET_dim(1), PET_dim(2));
Qy_ind = zeros(PET_dim(1), PET_dim(2));
% Grid index locations
G = zeros(PET_dim(1), PET_dim(2));

DM = (M0(:,:,t2)-M0(:,:,t1))./(T(t2)-T(t1));
DM(DM==0)=nan;
% convert activity change to volume change
% DM = DM./Cm;
% normalize DM
% DM_norm = DM - nanmean(nanmean(DM));
% DM = DM_norm;
n=1;
nq = 0;
q = 1;
nend = 1;
for i = 1:PET_dim(1)
    q_in_row = length(find(~isnan(DM(i,:)))) - 1;
    nstart = nend;
    for j = 1:PET_dim(2)
        if ~isnan(DM(i,j))
            if j+1 <= PET_dim(2) && ~isnan(DM(i,j+1))
                Qx_ind(i,j) = q;
                A(n,q) = -1;
                q = q+1;
            end
            if j-1 > 0 && Qx_ind(i,j-1)~=0
                A(n, Qx_ind(i,j-1)) = 1;
            end
            % assign grid cell numbers
            G(i,j) = n;
            n=n+1;
        end
    end
    
    n=nstart;
    % assign y interfaces
    for j = 1:PET_dim(2)
        
        if ~isnan(DM(i,j)) && i+1 <= PET_dim(1) && ~isnan(DM(i+1,j))
            Qy_ind(i,j) = q;
            A(n,q) = -1;
            q = q+1;
        end
        
        if ~isnan(DM(i,j)) && i-1 >0 && ~isnan(DM(i-1,j))
            A(n, Qy_ind(i-1,j)) = 1;
        end
        
        if ~isnan(DM(i,j))
            n = n+1;
        end
       
    end
    nend = n;
end

% crop A
A = A(1:n-1, 1:q-1);

% figure
% imagesc(A)
% axis equal
% axis tight
% 
% figure
% imagesc(G)
% axis equal
% axis tight

% find b vector
DM_trans = DM';
b = DM_trans(:);
b(isnan(b))=[];

% Different options for solving for x, the matlab function built in seems
% to work well. The BACKSLASH SHOULD NOT BE USED TO CALCULATE THE SOLUTION
% x = pinv(A)*b;
x = lsqminnorm(A,b);
% x = A'*inv(A*A')*b;
% x = A'/(A*A')*b;

%% Verify the solution by comparing the calculated change in rad to the 
% measured change in rad
dm = A*x; 

DM_calc = nan(PET_dim(1), PET_dim(2));
xp = nan(PET_dim(1), PET_dim(2));
xm = nan(PET_dim(1), PET_dim(2));
yp = nan(PET_dim(1), PET_dim(2));
ym = nan(PET_dim(1), PET_dim(2));
for i = 1:PET_dim(1)
    for j = 1:PET_dim(2)
        if Qx_ind(i,j) ~=0
            xp(i,j) = x(Qx_ind(i,j));
        end
        
        if G(i,j) ~=0
            DM_calc(i,j) = dm(G(i,j));
        end
%         if Qy_ind(i,j) ~= 0
%             ym(i,j) = x(Qy_ind(i,j));
%         end
        
%         if j>1 && Qx_ind(i,j-1) ~=0
%             xm(i,j) = x(Qx_ind(i,j-1));
%         end
%         
        if i>1 && Qy_ind(i-1,j) ~= 0
            yp(i,j) = x(Qy_ind(i-1,j));
        end
    end
end

% plot measured change
figure('position', [618 558  1195 420])
s1 = subplot(1,2,1);
h1 = imagesc(DM);
title(['Change in normalized radioactivity'])
axis equal
axis tight
axis off
shading flat
colorbar
colormap(gca, flux_cc)
set(h1,'alphadata',~isnan(DM))
clim_m = caxis;



% plot flux-calculate change
s2 = subplot(1,2,2);
% h2 = imagesc(DM_calc);

mean_perm = nanmean(Kmd_mat,3);
h2 = imagesc(mean_perm(2:end-1, 2:end-1));
title(['Calculated radiotracer flux'])
axis equal
axis tight
axis off
shading flat
y1 = colorbar;
caxis([16.1484   30.5185])

set(h2,'alphadata',~isnan(DM))
colormap(gca, gray)
ylabel(y1, 'Permeability [mD]');
hold on
% draw box around area of interest
plot([10.5 10.5], [0.5 20.5], 'r', 'linewidth', 1) 
plot([9.5 9.5], [0.5 20.5], 'r', 'linewidth', 1) 
plot([10.5 9.5], [0.5 0.5], 'r', 'linewidth', 1) 
plot([10.5 9.5], [20.5 20.5], 'r', 'linewidth', 1) 

quiver(xp,yp, 'k', 'color', 'y')



% figure
% subplot(1,2,1)
% h1 = imagesc(xp);
% title(['Positive x flux'])
% axis equal
% axis tight
% axis off
% shading flat
% colorbar
% set(h1,'alphadata',~isnan(DM))
% 
% subplot(1,2,2)
% h1 = imagesc(yp);
% title(['Positive y flux'])
% axis equal
% axis tight
% axis off
% shading flat
% colorbar
% set(h1,'alphadata',~isnan(DM))
% hold on
% quiver(xp,yp, 'k')



% for i=30
%    figure(1)
%     center_slice = squeeze(squeeze(comp_PET_4D_coarse(:,9,:,i)));
%     center_slice(center_slice==0)=nan;
%    h = imagesc(center_slice);
%     axis equal
%     axis tight
%     axis off
%     shading flat
%     title(['Center slice activity location'])
%     colorbar
%     set(h,'alphadata',~isnan(center_slice))
%     caxis([0 0.3])
%     drawnow
% %     pause(0.1)
% %     subplot(2,1,2)
% figure(2);
%     percent_diff = (M0(:,:,i)-M0(:,:,19))./M0(:,:,19).*100;
%     change_mass = M0(:,:,i)-M0(:,:,i-2);
%     change_mass(change_mass==0)=nan;
%     h1 = imagesc(change_mass);
%     title(['Change in radioactivity, frame: ', num2str(i)])
%     axis equal
%     axis tight
%     axis off
%     shading flat
% %     title(['Percent change in radioactivity, frame: ', num2str(i)])
%     colorbar
%     set(h1,'alphadata',~isnan(change_mass))
% %     caxis([-40 30])
% %     drawnow
% %     pause(0.8)
% end
%
% figure
% hist(change_mass(:), 20)