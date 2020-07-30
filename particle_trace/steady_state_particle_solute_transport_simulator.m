% steady_state_particle_solute_transport_simulator
% Christopher Zahasky
% 1/9/2019

close all
clear all
set(0,'DefaultAxesFontSize',14, 'defaultlinelinewidth', 1,...
    'DefaultAxesTitleFontWeight', 'normal')

% This is the Matlab code simulating particle transport in 2D given a
% velocity field. This version does not account for dispersion or diffusion

% adjust paths depending on computer
current_folder = pwd;
str_index = strfind(pwd, '\Dropbox');

% Path to colorbrewer
addpath([current_folder(1:str_index),'Dropbox\Matlab\high_res_images'])
% Path to functions
addpath([current_folder,'\functions'])

% load perm map
load('streamtube_perm_field')

%% Load velocity map input data
load('ss_10kpa_stream_perm_velocity_field')
% Divide to account for porosity
VEL.U = VEL.U./0.2;
VEL.V = VEL.V./0.2;
% Define grid
Grid.xe = [0: Grid.nx-3].*Grid.dx;
Grid.ye = [0: Grid.ny]'.*Grid.dy;
VEL.U = VEL.U(:,1:end-1);
%% Convergence and timestepping input
% start time of simulation (seconds)
T.start_time = 0;

% Particle simulation will terminate if particle reaches boundary, the
% end_time is reached, or the maximum number of timesteps is exceeded
% maximum time
T.end_time = 50000;
% timestep limit
T.max_timesteps = 1000;
% save times (seconds)
T.save_times = [100:100:T.end_time];
tic
% model tolerance
model_tol = eps*100;
% if particle is within boundary tolerance of the model edge then terminate
boundary_tolerance = min([Grid.dx, Grid.dy])/10;

%% Generate random particle
% Number of particles to be generated
P.total_particles = 100;
% particle distribution options
% if set equal to '1' then the particles will be evenly distributed
% in region defined by P.xreg and P.yreg
% if set to '2' then particles will be evenly distributed in
% y-direction and as defined by the pdf in the x-direction
P.particle_dist = 1;

% x region (min max coordinates)
P.xreg = [0 0.001];
% y region (min max coordinates)
P.yreg = [0 Grid.dy*Grid.ny];
% Call particle distribution function
P = particle_distribution_function(P, Grid, 0);

% preallocate saved particle location
P.xsave = zeros(P.total_particles, length(T.save_times)+1);
P.ysave = zeros(P.total_particles, length(T.save_times)+1);
% save initial location
P.xsave(:,1) = P.x;
P.ysave(:,1) = P.y;

% Visual check of initial particle distribution
figure
hold on
gridX = [1:Grid.nx].*Grid.dx - (Grid.dx/2);
gridY = [1:Grid.ny].*Grid.dy - (Grid.dy/2);
imagesc(gridX, gridY, perm_profile_md)
colormap(gray)
hold on
xlabel('Distance from inlet [m]')
axis equal
axis tight
plot(P.x, P.y, 'or')
% Option to plot grid
% for gx = 1:Grid.nx+1
%     plot([Grid.xe(gx) Grid.xe(gx)], [0 max(Grid.ye)], '-k', 'linewidth', 1)
% end
% for gy = 1:Grid.ny+1
%     plot([0 max(Grid.xe)], [Grid.ye(gy) Grid.ye(gy)], '-k', 'linewidth', 1)
% % end
axis equal
axis([0 0.1 0 max(Grid.ye)])
box on
set(gca, 'Ydir', 'reverse')
drawnow

%% Calculate velocity gradient constants (A terms in MODPATH guide)
% right face minus left face
VEL.Ax = (VEL.U(:, 2:end) - VEL.U(:, 1:end-1))./Grid.dx;
% bottom face minus top face (DOWN IS POSITIVE)
VEL.Ay = (VEL.V(2:end, :) - VEL.V(1:end-1, :))./Grid.dy;

%% Calculate particle velocity
for i = 1:P.total_particles
    % reset time tracking stuff
    T.timestep_counter = 0;
    save_time_index = 1;
    T.particle_run_time = 0;
    % reset while loop switch
    particle_active = 0;
    
    while particle_active == 0
        % Find grid cell where particle is located
        [P, VEL, particle_active] = locate_particle_function ...
            (P, VEL, Grid, i, particle_active, model_tol);
        
        % If particle is in last grid cell then terminate that particle
        if P.gj(i) > Grid.nx
            continue
        end
            
        if particle_active == 1
            % end current loop
            continue
        end
        
        % Calculate interpolated particle velocity and how long it will
        % take for particle to reach grid cell edge
        [P, VEL, dte] = time_to_edge_function(P, VEL, Grid, i);
        
        if isreal(dte) == 0
            wtf=0
        end
        
        % track particle time
        [T, dte, save_time_index, particle_active] = update_time_function...
            (T, dte, save_time_index, particle_active);

        % calculate new location
        [P, VEL] = new_particle_location_function(P, VEL, i, dte);
        
        % transverse diffusion component
        random_vec = normrnd(0,1);
        diff_y = random_vec*sqrt(1e-10*dte);
        % Displace particle by small diffusive component
        pye_temp = P.ye(i) + diff_y;
        % Check if new particle location is outside of model
        % if particle is above model then don't add diffusion and
        % dispersion terms
        if pye_temp > 0 && pye_temp < Grid.ye(end)
            P.ye(i) = pye_temp;
        else % if outside model then reflect motion
            P.ye(i) = P.ye(i) - diff_y;
        end
        
        % longitudinal diffusion component
        random_vec = normrnd(0,1);
        diff_x = random_vec*sqrt(1e-10*dte);
        % Displace particle by small diffusive component
        pxe_temp = P.xe(i) + diff_x;
        % Check if new particle location is outside of model
        % if particle is past outlet then set to outlet grid location
        if pxe_temp > 0 && pxe_temp < Grid.xe(end)
            P.xe(i) = pxe_temp;
        elseif pxe_temp > Grid.xe(end)
            P.xe(i) = Grid.xe(end);
        end
        
        % update location
        P.x(i) = P.xe(i);
        P.y(i) = P.ye(i);
        
        % Check timestep limits
        if T.timestep_counter > T.max_timesteps
            particle_active = 1;
            T.save_data = 1;
            i
            warning('Particle trajectory incomplete, reached max timesteps')
        end
        
        % Save some data
        if T.save_data == 1
            plot([P.xsave(i, save_time_index-1), P.x(i)], ...
                [P.ysave(i, save_time_index-1),P.y(i)], '-b')
            % save initial location
            P.xsave(i, save_time_index) = P.x(i);
            P.ysave(i, save_time_index) = P.y(i);
        end
        
        T.timestep_counter = T.timestep_counter +1;
        
        % Test if particle is near a boundary
%         ybound = ismembertol(P.y(i), [Grid.ye(1), Grid.ye(end)], boundary_tolerance);
        xbound = ismembertol(P.x(i), [Grid.xe(1), Grid.xe(end)], boundary_tolerance);
        % if particle is at the edge of the model then terminate now
        if xbound
            particle_active = 1;
        end
    end
end

% Plot final location
plot(P.x, P.y, 'og', 'Markersize', 8)
axis([0 max(Grid.xe) 0 max(Grid.ye)])

% Save data
% savefilename = [num2str(P.total_particles),'_',...
%     'particles_SS_stream_perm_1e10_advect'];
% 
% save(savefilename, 'P', 'T', 'Grid')

% Display run time length
toc