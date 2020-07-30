function [P, VEL, particle_active] = locate_particle_function ...
    (P, VEL, Grid, i, particle_active, model_tol)
% Christopher Zahasky
% Stanford University
% 1/9/2019
% Calculate particle distrubution
% Find grid cell where particle is located, if on a line then choose
% cell that particle is moving into, if on a corner then choose the
% cell that the particle is moving into.
% Test if particle is near a boundary
[ygrid_face, grid_row] = ismembertol(P.y(i), Grid.ye, model_tol);
[xgrid_face, grid_col] = ismembertol(P.x(i), Grid.xe, model_tol);

% set particle location at exactly the boundary if found to be close,
% except if it is near the boundary
if ygrid_face && grid_row ~=1 && grid_row ~= Grid.ny+1 
    P.y(i) = Grid.ye(grid_row);
end

if xgrid_face
    P.x(i) = Grid.xe(grid_col);
end

% if particle is on a corner
if ygrid_face && xgrid_face
    % if particle is at the edge of the model then terminate now
    if grid_row == Grid.ny+1 || grid_col == Grid.nx+1
        particle_active = 1;
        % end current loop
        return
    end
    
    VEL.px(i) = VEL.U(grid_row, grid_col);
    VEL.py(i) = VEL.V(grid_row, grid_col);
    
    
    if VEL.px(i)>0 && VEL.py(i)>0
        % Particle in correct location
        P.gi(i) = grid_row;
        % Particle in correct location
        P.gj(i) = grid_col;
        
    elseif VEL.px(i)<0 && VEL.py(i)>0
        % Particle in correct location
        P.gi(i) = grid_row;
        % Particle moves left
        P.gj(i) = grid_col - 1;
        
    elseif VEL.px(i)>0 && VEL.py(i)<0
        % Particle moves up a row
        P.gi(i) = grid_row - 1;
        % Particle in correct location
        P.gj(i) = grid_col;
        
    elseif VEL.px(i)<0 && VEL.py(i)<0
        % Particle moves up a row
        P.gi(i) = grid_row - 1;
        % Particle moves left
        P.gj(i) = grid_col - 1;
    end
    % if particle lies on y-face of grid cell
elseif ygrid_face
    %         grid_row = find(P.y(i)==Grid.ye);
    
    % if particle is at the edge of the model then terminate now
%     if grid_row == Grid.ny+1
%         particle_active = 1;
%         % end current loop
%         return
%     end
    
    % Find column
    P.gj(i) = find(P.x(i) > Grid.xe, 1, 'last');
    VEL.py(i) = VEL.V(grid_row, P.gj(i));
    
    % Find row
    if VEL.py(i)>0
        % Particle in correct location
        P.gi(i) = grid_row;
        
    elseif VEL.py(i)<0
        % Particle moves up a row
        P.gi(i) = grid_row - 1;
    end
    
    % if particle lies on x-face of grid cell
elseif xgrid_face
    % if particle is at the edge of the model then terminate now
    if grid_col == Grid.nx+1
        particle_active = 1;
        % end current loop
        return
    end
    
    % Find row
    P.gi(i) = find(P.y(i) > Grid.ye, 1, 'last');
    
    VEL.px(i) = VEL.U(P.gi(i), grid_col);
    
    % Find column
    if VEL.px(i)>0
        % Particle in correct location
        P.gj(i) = grid_col;
        
    elseif VEL.px(i)<0
        % Particle moves left
        P.gj(i) = grid_col - 1;
    end
    
    % else particle lies in a grid cell
else
    % Find row
    P.gi(i) = find(P.y(i) > Grid.ye, 1, 'last');
    % Find column
    P.gj(i) = find(P.x(i) > Grid.xe, 1, 'last');
end

