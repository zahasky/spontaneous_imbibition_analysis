% MRST_imbibe_model
clear all
close all

% formating standards
set(0,'DefaultAxesFontSize',14, 'defaultlinelinewidth', 2,...
    'DefaultAxesTitleFontWeight', 'normal')

% NOTE that this script requires both colorbrewer: https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
% and the Matlab reservoir simulation toolbox: https://www.sintef.no/projectweb/mrst/download/
% Define colormaps for plotting
color_inc = 1000;
% pressure
pubar = cbrewer('div', 'Spectral', 1000 , 'PCHIP');
cbar = [1 1 1; pubar(500:1000,:)];
% effective permeability
gray_map = [1 1 1; (cbrewer('seq', 'Greys', (color_inc-1) , 'PCHIP'))];
% velocity
vbar = [1 1 1; flipud(cbrewer('div', 'RdYlGn', (color_inc-1) , 'PCHIP'))];
% capillary pressure colorbar
bbar = jet(5);
    
%% Add paths to functions and data
addpath('C:\Users\zahas\Dropbox\Matlab\high_res_images')
% add path to MRST codes
addpath('C:\Users\zahas\Dropbox\Matlab\MRST_2020a\mrst-2020a\modules\incomp\fluid\incompressible')
addpath('C:\Users\zahas\Dropbox\Matlab\MRST_2020a\mrst-2020a\modules\incomp')
addpath('C:\Users\zahas\Dropbox\Matlab\MRST_2020a\mrst-2020a\modules\incomp\fluid\utils')
% Path for grid utilities
addpath('C:\Users\zahas\Dropbox\Matlab\MRST_2020a\mrst-2020a\modules\co2lab\utils\grid')
% Path to startup
% addpath('C:\Users\zahas\Dropbox\Matlab\MRST_2020a\mrst-2020a')
%  run startup

% CT saturation data
load('ct_scan_imbibe_times')
% Load model permeability data
load('2D_perm_data_input')

%% Experiment input
inlet_pressure = 0; % psi

% Number of grid cells along core axis
Grid.nx = 500;

% porosity
Core.phi = 0.2;
% relative permeability parameters
Core.swir = 0;
Core.krwmax = 1;
Core.lam_w = 7.6;

% Set up intrinsic permeability map
Core.k = repmat(nanmean(Km2_mat(2:end-1,11,:),3), [1, Grid.nx]);
[Grid.ny, Grid.nx] = size(Core.k);

%% Set up model geometry
area = pi*2.5^2;
equiv_h = area/(0.2329*Grid.ny);

[nx, ny, nz] = deal(Grid.ny,  Grid.nx, 1);
[Dx, Dy, Dz] = deal(0.2329/100*Grid.ny, 0.1, equiv_h/100);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);
G.size = [Dx, Dy, Dz];

%% Define fluid stuff
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
gravity off

%% Set up imbibition loop
% Load CT wet data for calculating saturation
load full_sat_sa_profile.mat
Saturation.CT_wr = squeeze(squeeze(SAVG));
% Load dry air CT data
load dry_sa_profile.mat
Saturation.CT_ar = squeeze(squeeze(SAVG));

nn = 1;
for scan_n = [46 76 118 178]
    % Set scan number
    Saturation.scan_num = scan_n;
    % Calculate absolute scan number (CT scanner counts every 3 and the
    % starting scan for the experiment was scan 19)
    abs_scan_num = length(19:3:scan_n);
    % Extract scan time in seconds
    Saturation.scan_time = t_si_dry(abs_scan_num)*60;
    % Set polyfit order for fitting water saturation profile
    polyfit_order = 5;
    
    figure(1) % plot for 1D saturation profile
    hold on
    % Call saturation calculation function and plot profile
    [Saturation, p, G] = imb_sat_profile_MRST_function(Saturation, G);
    drawnow
    
    % Calculate the product of porsity times saturation
    sat_porosity = repmat(Saturation.Sw, [Grid.ny, 1]).* Core.phi;
    
    % Build matrix from saturation profile
    Sat_matrix = repmat(Saturation.Sw, [Grid.ny, 1]);
    % Define relative permeability based on saturation map
    Sw_eff = (Sat_matrix - Core.swir)/(1-Core.swir);
    kw = Core.krwmax.*(Sw_eff).^Core.lam_w;
    % Effective permeability map
    k_effective = Core.k .*kw;
    
    %% Define rock model
    % For this problem we use a rock of constant compressibility.  The
    % pore-volume is therefore a simple analytic function of fluid pressure.
    rock = makeRock(G, k_effective(:), sat_porosity(:));
    % Visual check
    % figure
    % subplot(1,2,1)
    % plotCellData(G, rock.perm/(milli*darcy)), colorbar, axis equal, axis tight
    % title('Effective permeability [mD]')
    %
    % subplot(1,2,2)
    % plotCellData(G, rock.poro), colorbar, axis equal, axis tight
    % title('Saturated porosity')
    
    
    %% Imbibition Conditions
    % apply imbibition conditions, this is calculated from function
    % identified for this sample in Zahasky and Benson (2019)
    % https://onlinelibrary.wiley.com/doi/abs/10.1029/2019GL084532
    ds_dlamda = -2.2677E3.*Saturation.Sw.^5 + 2.6112E3.*Saturation.Sw.^4 + -0.2813E3.*Saturation.Sw.^3 ...
        -0.5815E3.*Saturation.Sw.^2 + 0.1810E3.*Saturation.Sw;
    % fluid imbibing only where there the max saturation isn't reached
    ds_dlamda(Saturation.Sw == Saturation.max_sat) = 0;
    % scale ds_dlamba back to dSw_dt
    ds_dt = ds_dlamda.*G.d2 ./4.*Saturation.scan_time^(-1.25);
    % convert ds_dt to m^3/sec
    imbibe_along_core = ds_dt .*((Dy/ny)*(Dx/nx)*(Dz/nz)*Core.phi);
    % total flow rate at this time
    inlet_imbibition_rate_ml_min = sum(imbibe_along_core)*100^3*60*Grid.ny;
    
    %% Boundary Conditions
    % Inlet
    bc  = pside([], G, 'YMin', inlet_pressure.*psia());
    
    % Find bottom face index
    facesbottom = boundaryFaceIndices(G, 'Bottom', [], [], []);
    % Set imbibition rate as sink in out of plane dimension to replicate
    % flow rate conditions
    for row_from_inlet = 1:G.cartDims(2)
        if imbibe_along_core(row_from_inlet)>0
            spec_faces = facesbottom(1+(row_from_inlet-1)*G.cartDims(1): ...
                (row_from_inlet)*G.cartDims(1));
            bc = addBC(bc, spec_faces, 'flux', -imbibe_along_core(row_from_inlet));
        end
    end
    
    %% Assemble and solve the linear system
    % To solve the flow problem, we use the standard two-point
    % flux-approximation method (TPFA), which for a Cartesian grid is the same
    % as a classical seven-point finite-difference scheme for Poisson's
    % equation. This is done in two steps: first we compute the
    % transmissibilities and then we assemble and solve the corresponding
    % discrete system.
    T   = computeTrans(G, rock);
    sol = incompTPFA(initResSol(G, 0.0, 1), G, T, fluid, 'bc', bc); 
    
    %% VELOCITY
    % Convert from flux velocity to advection velocity
    % Direction parallel to imbibition
    % Extract data from MRST data structures
    facesf = boundaryFaceIndices(G, 'Front', [], [], []);
    facesb = boundaryFaceIndices(G, 'Back', [], [], []);
    U_flux_faces = [facesb(1): facesf(end)];
    % calculate area of face
    uface_area = (Dy/ny)*(Dz/nz);
    % Save and plot in particle tracking code format
    VEL.U = reshape(sol.flux(U_flux_faces)./uface_area, [G.cartDims(1), G.cartDims(2)+1]);
    sat_porosity(sat_porosity < 0.00001) = nan;
    VEL.U = VEL.U(:,1:end)./[sat_porosity, zeros(G.cartDims(1), 1)];
    
    % Direction perpendicular to imbibition
    % Extract data from MRST data structures
    facesr = boundaryFaceIndices(G, 'Right', [], [], []);
    V_flux_faces = [1:facesr(end)];
    vface_area = (Dx/nx)*(Dz/nz);
    VEL.V = reshape(sol.flux(V_flux_faces)./vface_area, [G.cartDims(1)+1, G.cartDims(2)]);
    VEL.V = VEL.V(1:end,:)./[zeros(1,G.cartDims(2)); sat_porosity]; % meters/sec
    % Reenforce no flow at top and bottom boundaries 
    VEL.V(1,:) = 0;
    VEL.V(end,:) = 0;
    
    %% PLOT RESULTS
    % cell centered grid in y direction
    G.d1 = linspace(G.size(1) /G.cartDims(1)/2, G.size(1), G.cartDims(1))*100;
    % Define coarsened meshgrid for plotting velocity vectors
    xinc = 40;
    % Define meshgrid
    [X, Y] = meshgrid(G.d2(1:xinc:end-xinc), G.d1);
    
    % Calculate y velocity range for colormaps
    velocity_y = (VEL.V(1:end-1, :)+VEL.V(2:end, :))/2;
    velocity_yc = velocity_y(:, 1:xinc:end-xinc);
    max_vy = max(abs(velocity_yc(:)));
    vy_range = linspace(-max_vy, max_vy, color_inc);
    
    % Calculate x velocity range for colormaps
    velocity_x = (VEL.U(:, 1:end-1)+VEL.U(:, 2:end))/2;
    velocity_x(isinf(velocity_x)) = 0;
    velocity_xc = velocity_x(:, 1:xinc:end-xinc);
    
    f2 = figure(2);
    % Effective permeability
    subplot(4,2,2*(nn-1)+1)
    imagesc(G.d2, G.d1, k_effective/(milli*darcy));
    hold on
    axis equal
    axis tight
    colormap(gca, gray_map)
    caxis([0 0.7])
    axis_limits = axis;
    if nn == 1
        title('Advective velocity')
    end
    
    % Extract size of coarsend velocity field for looping
    [di, dj] = size(velocity_yc);
    % Loop to plot quiver vectors so that they can by colored by orthogonal
    % velocity compoenent
    for i = 1:di
        for j = 1:dj
            if isnan(velocity_yc(i,j)) ~= 1
                color_index = find(vy_range > velocity_yc(i,j), 1, 'first');
                quiver(X(i,j), Y(i,j), velocity_xc(i,j)./500, velocity_yc(i,j), 3.6E6, 'color', vbar(color_index,:), 'linewidth', 1.2, 'MaxHeadSize', 1);
            end
        end
    end
    
    % Plot corresponding 2D pressure field
    subplot(4,2,2*(nn-1)+2)
    pressure_psi = reshape(convertTo(sol(end).pressure, psia()), [G.cartDims(1), G.cartDims(2)]);
    imagesc(G.d2, G.d1, pressure_psi.*6.89)
    axis equal
    axis tight
    if nn == 1
        title('Fluid pressure')
    end
    caxis([-300, 0])
    set(gca,'ytick',[])
    colormap(gca, cbar)
    
    % Plot capillary pressure as a function of saturation
    figure(4)
    hold on
    plot(Sat_matrix(:), -pressure_psi(:), 'o', 'color', bbar(nn,:), 'linewidth', 0.5)
    
    % Update counter
    nn = nn + 1;
end
% Format capillary pressure plot
nn = nn +1;
ylim([0 50])
xlim([0 1])
xlabel('Water saturation')
ylabel('Capillary Pressure [psi]')
title('Capillary pressure calculated from model')
box on


% Lots of plot formating
figure(2)
set(f2, 'position', [734    63   964   916])
subplot(4,2,7)
xlabel('Distance from inlet [cm]')
y1 = colorbar('Position',...
    [0.0993472796460112 0.141921397379913 0.0180425011762578 0.779666150161994]);
ylabel(y1, 'Effective permeability [mD]', 'fontsize', 16);

subplot(4,2,1)
h3 = text('FontSize',14,'Rotation',90,...
    'String','Orthogonal advective velocity [m/s]',...
    'Position',[12.6183272949099 16.5016505744996 0]);

ax2 = axes('Position',[0.13 0.11 0.292284644194757 0.815]);
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax2,vbar)
colorbar(ax2, 'Position',...
    [0.475343320848939 0.14 0.019038701622971 0.778888888888889])
caxis([-max_vy, max_vy])

subplot(4,2,8)
y2 = colorbar('Position',...
    [0.91010151522202 0.144086021505376 0.0187381658974735 0.774193548387097]);
ylabel(y2, '[kPa]', 'fontsize', 16);
xlabel('Distance from inlet [cm]')

