function [P, VEL, dte] = time_to_edge_function(P, VEL, Grid, i)
% Christopher Zahasky
% Stanford University
% 1/9/2019
% Calculate particle velocity and how much time it will take to reach edge
% of current grid cell

% Calculate interpolated particle velocity
VEL.px(i) = VEL.Ax(P.gi(i), P.gj(i)) *(P.x(i) - Grid.xe(P.gj(i))) + ...
    VEL.U(P.gi(i), P.gj(i));

VEL.py(i) = VEL.Ay(P.gi(i), P.gj(i)) *(P.y(i) - Grid.ye(P.gi(i))) + ...
    VEL.V(P.gi(i), P.gj(i));

% calculate time to reach grid cell face (depends on whether velocity
% is positive or negative)
if VEL.px(i)>0
    % dx = 1/Ax ln(Vx2/Vxp)
    dtx = (1/VEL.Ax(P.gi(i), P.gj(i)))* ...
        log(VEL.U(P.gi(i), P.gj(i)+1) / VEL.px(i));
    % If there is no change in gradient then just linearly interpolate
    if isnan(dtx) && VEL.Ax(P.gi(i), P.gj(i))== 0
        dtx = (Grid.xe(P.gj(i)+1) -P.x(i)) /VEL.px(i);
    end
    
%     if isreal(dtx) == 0
%         dtx = (Grid.xe(P.gj(i)+1) -P.x(i)) /VEL.px(i);
%     end
        
else
    % dx = 1/Ax ln(Vx1/Vxp)
    dtx = (1/VEL.Ax(P.gi(i), P.gj(i)))*...
        log(VEL.U(P.gi(i), P.gj(i)) / VEL.px(i));
    % If there is no change in gradient then just linearly interpolate
    if isnan(dtx) && VEL.Ax(P.gi(i), P.gj(i))== 0
        dtx = (Grid.xe(P.gj(i)) - P.x(i)) /VEL.px(i);
    end
%     if isreal(dtx) == 0
%         dtx = (Grid.xe(P.gj(i)) - P.x(i)) /VEL.px(i);
%     end
end

if VEL.py(i)>0
    % dty = 1/Ay ln(Vy2/Vyp)
    dty = (1/VEL.Ay(P.gi(i), P.gj(i))) *...
        log(VEL.V(P.gi(i)+1, P.gj(i)) / VEL.py(i));
    if isnan(dty) && VEL.Ay(P.gi(i), P.gj(i))== 0 || isinf(dty)
        %%% This is used for SS
        dty = (Grid.ye(P.gi(i)+1) - P.y(i)) /VEL.py(i);
        %%% This is used for Imbibition
%         dty = (Grid.ye(P.gi(i)) - P.y(i)) /VEL.py(i);
    end
else
    % dty = 1/Ay ln(Vy1/Vyp)
    dty = (1/VEL.Ay(P.gi(i), P.gj(i))) *...
        log(VEL.V(P.gi(i), P.gj(i)) / VEL.py(i));
    if isnan(dty) && VEL.Ay(P.gi(i), P.gj(i))== 0 || isinf(dty)
        dty = (Grid.ye(P.gi(i)) - P.y(i)) /VEL.py(i);
    end
end

% need to build a check and timestep cut if dty is imaginary, for now
% take the norm
dty = abs(dty);
dtx = abs(dtx);
% assign exit time to be the shorter of the two times
if dtx > dty
    % particle will enter cell i+-1, j
    dte = dty;
elseif dtx < dty
    % particle will enter cell i, j+-1
    dte = dtx;
else
    % particle will enter cell i+-1, j+-1
    dte = dtx;
end