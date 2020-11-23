% MRST_SS_model
% Christopher Zahasky
% 11/23/2020
clear all
close all

% formating standards
set(0,'DefaultAxesFontSize',14, 'defaultlinelinewidth', 2,...
    'DefaultAxesTitleFontWeight', 'normal')

% NOTE that this script requires both colorbrewer: https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
% and the Matlab reservoir simulation toolbox: https://www.sintef.no/projectweb/mrst/download/

% add path to MRST codes
addpath('C:\Users\zahas\Dropbox\Matlab\MRST_2020a\mrst-2020a')
% % run startup
% startup

% Load model permeability data
load('2D_perm_data_input')
Kmap_streamtube = repmat(nanmean(Km2_mat(2:end-1,11,:),3), [1, 42]);
ksize = size(Kmap_streamtube);

% porosity 
porosity = 0.20;

%% Experiment input
back_pressure = 7.2569; % psi
injection_rate = 2; % mL/min

%% Set up model geometry
area = pi*2.5^2;
equiv_h = area/(0.2329*ksize(1));

[nx, ny, nz] = deal( ksize(1),  ksize(2), 1);
[Dx, Dy, Dz] = deal(0.2329/100*ksize(1), 0.2329/100*ksize(2), equiv_h/100);
G = cartGrid([nx, ny, nz], [Dx, Dy, Dz]);
G = computeGeometry(G);
G.size = [Dx, Dy, Dz];

%% Define rock model
% For this problem we use a rock of constant compressibility.  The
% pore-volume is therefore a simple analytic function of fluid pressure.
rock = makeRock(G, Kmap_streamtube(:), porosity);

%% Define fluid stuff
fluid = initSingleFluid('mu', 1*centi*poise, 'rho', 1014*kilogram/meter^3);
gravity off

%% Boundary Conditions
% Outlet
bc  = pside([], G, 'YMax', back_pressure.*psia());
% Inlet
bc  = pside(bc, G, 'YMin', back_pressure.*psia().*3);

%% Assemble and solve the linear system
% To solve the flow problem, we use the standard two-point
% flux-approximation method (TPFA), which for a Cartesian grid is the same
% as a classical seven-point finite-difference scheme for Poisson's
% equation. This is done in two steps: first we compute the
% transmissibilities and then we assemble and solve the corresponding
% discrete system.
T   = computeTrans(G, rock);
sol = incompTPFA(initResSol(G, 0.0, 1), G, T, fluid, 'bc', bc);

%% Extract advection velocity
% extract face locations for faces parallel to axis of core
facesf = boundaryFaceIndices(G, 'Front', [], [], []);
facesb = boundaryFaceIndices(G, 'Back', [], [], []);
U_flux_faces = [facesb(1): facesf(end)];
% calculate area of face
uface_area = (Dy/ny)*(Dz/nz);
% Save and plot in particle tracking code format
VEL.U = reshape(sol.flux(U_flux_faces)./uface_area, [G.cartDims(1), G.cartDims(2)+1]);
% Convert from flux velocity to advection velocity
VEL.U = VEL.U./porosity; % meters/sec

% extract face locations for faces perpendicular to axis of core
facesr = boundaryFaceIndices(G, 'Right', [], [], []);
V_flux_faces = [1:facesr(end)];
vface_area = (Dx/nx)*(Dz/nz);
VEL.V = reshape(sol.flux(V_flux_faces)./vface_area, [G.cartDims(1)+1, G.cartDims(2)]);
% Convert from flux velocity to advection velocity
VEL.V = VEL.V./porosity; % meters/sec
% Re-enforce no flux boundary boundaries at core walls
VEL.V(1,:) = 0;
VEL.V(end,:) = 0;

%% PLOT RESULTS
color_inc = 1000;
% Set colormaps
% pressure
pubar = cbrewer('div', 'Spectral', 1000 , 'PCHIP');
cbar = [pubar(500:1000,:)];
% permeability
gray_map = [1 1 1; (cbrewer('seq', 'Greys', (color_inc-1) , 'linear'))];
% velocity
vbar = [1 1 1; flipud(cbrewer('div', 'RdYlGn', (color_inc-1) , 'PCHIP'))];

% cell centered grid in y direction
G.d1 = linspace(G.size(1) /G.cartDims(1)/2, G.size(1), G.cartDims(1))*100;
% define cell centered grid in dimension 2
G.d2 = linspace(G.size(2) /G.cartDims(2)/2, G.size(2), G.cartDims(2))*100;
% Define coarsened meshgrid for plotting velocity vectors
xinc = 4;
[X, Y] = meshgrid(G.d2(1:xinc:end-xinc), G.d1);

% Calculate y velocity range for colormaps
velocity_y = (VEL.V(1:end-1, :)+VEL.V(2:end, :))/2;
velocity_yc = velocity_y(:, 1:xinc:end-xinc);
max_vy = max(abs(velocity_yc(:)));

% Calculate x velocity range for colormaps
velocity_x = (VEL.U(:, 1:end-1)+VEL.U(:, 2:end))/2;
velocity_x(isinf(velocity_x)) = 0;
velocity_xc = velocity_x(:, 1:xinc:end-xinc);
max_vx = max(abs(velocity_xc(:)));
vy_range = linspace(-max_vx./10000, max_vx./10000, color_inc);

figure('position', [34, 558, 1367, 420])
% Effective permeability
subplot(1,2,1)
imagesc(G.d2, G.d1, Kmap_streamtube/(milli*darcy));
hold on
axis equal
axis tight
colormap(gca, gray_map)
caxis([10 32])
axis_limits = axis;
title('Advective velocity')
y1 = colorbar;
ylabel(y1, 'Permeability [mD]', 'fontsize', 16);
xlabel('Distance from inlet [cm]')

% Extract size of coarsend velocity field for looping
[di, dj] = size(velocity_yc);
% Plot velocity field, the loop enables the color of every vector to be set
% by the orthogonal component of the velocity
for i = 1:di
    for j = 1:dj
        if isnan(velocity_yc(i,j)) ~= 1
            color_index = find(vy_range > velocity_yc(i,j), 1, 'first');
            if isempty(color_index)
               quiver(X(i,j), Y(i,j), velocity_xc(i,j)./500, velocity_yc(i,j), 2E6, 'color', vbar(end,:), 'linewidth', 1.2, 'MaxHeadSize', 1); 
            else
                quiver(X(i,j), Y(i,j), velocity_xc(i,j)./500, velocity_yc(i,j), 2E6, 'color', vbar(color_index,:), 'linewidth', 1.2, 'MaxHeadSize', 1);
            end
        end
    end
end

% Pressure
subplot(1,2,2)
% Plot 2D pressure field
pressure_psi = reshape(convertTo(sol(end).pressure, psia()), [G.cartDims(1), G.cartDims(2)]);
imagesc(G.d2, G.d1, pressure_psi.*6.89)
axis equal
axis tight
title('Fluid pressure')
colormap(gca, cbar)
y2 = colorbar;
ylabel(y2, '[kPa]', 'fontsize', 16);
xlabel('Distance from inlet [cm]')

