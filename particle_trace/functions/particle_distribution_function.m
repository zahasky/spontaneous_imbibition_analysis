function [P] = particle_distribution_function(P, Grid, C)
% Christopher Zahasky
% Stanford University
% 1/9/2019
% Calculate particle distrubution

% particle distribution options
% if set equal to '1' then the particles will be evenly distributed
% in region defined by P.xreg and P.yreg
% if set to '2' then particles will be evenly distributed in
% y-direction and as defined by the pdf in the x-direction
% In this case the C input is a saturation profile
% if set to '3' then particles will be distributed based on a 2D
% concentration map defined by the C input.
% In this case the C input is a matrix of x,y locations and concentrations

if P.particle_dist == 1
    % particle x coordinates
    P.x = rand(P.total_particles, 1).*diff(P.xreg)+ P.xreg(1);
    % particle y coordinates
    P.y = rand(P.total_particles, 1).*diff(P.yreg)+ P.yreg(1);
    
elseif P.particle_dist == 2
    
    area_under_curve = sum(C(1:Grid.nx).*Grid.dx);
    particle_pdf = zeros(1, Grid.nx);
    for i = 1:Grid.nx
        particle_pdf(i) = sum(C(1:i).*Grid.dx)/area_under_curve;
    end
    figure
    hold on
    plot(particle_pdf, Grid.xe(1:end-1)+Grid.dx/2)
    p = polyfit(particle_pdf,Grid.xe(1:end-1)+Grid.dx/2, 6);
    % generate random sampling of pdf
    rand_points = rand(P.total_particles, 1);
    % particle x coordinates
    P.x = polyval(p,rand_points);
    P.x(P.x<0) = Grid.dx/2;
    plot(rand_points, P.x, 'ok')
    
    % particle y coordinates
    P.y = rand(P.total_particles, 1).*diff(P.yreg)+ P.yreg(1);
    
elseif P.particle_dist == 3
    
    
else
    error('P.particle_dist must be defined as equal to 1 (square pulse) or 2 (custom pdf)')
end

% preallocate coordinate matrix
P.gi = ones(P.total_particles, 1);
P.gj = ones(P.total_particles, 1);
% preallocate next particle location
P.xe = zeros(P.total_particles, 1);
P.ye = zeros(P.total_particles, 1);

