function [P, VEL] = new_particle_location_function(P, VEL, i, dte)
% Christopher Zahasky
% Stanford University
% 1/9/2019
% Calculate new particle location

% calculate new location
if VEL.Ax(P.gi(i), P.gj(i))== 0
    % If no gradient use linear calculation
    P.xe(i) = P.x(i) + dte.*VEL.px(i);
else
    % Use Goode formulation
    P.xe(i) = P.x(i)+ (VEL.px(i)/VEL.Ax(P.gi(i), P.gj(i)))*...
        (exp(VEL.Ax(P.gi(i), P.gj(i))* dte) - 1);
end

if VEL.Ay(P.gi(i), P.gj(i))== 0
    % If no gradient use linear calculation
    P.ye(i) = P.y(i) + dte.*VEL.py(i);
else
    % Use Goode formulation
    P.ye(i) = P.y(i)+ (VEL.py(i)/VEL.Ay(P.gi(i), P.gj(i)))* ...
        (exp(VEL.Ay(P.gi(i), P.gj(i))* dte) - 1);
end