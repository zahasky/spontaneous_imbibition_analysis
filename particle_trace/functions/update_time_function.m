function [T, dte, save_time_index, particle_active] = update_time_function...
    (T, dte, save_time_index, particle_active)
% Christopher Zahasky
% Stanford University
% 1/9/2019
% Check if timestep needs to be cut to save data or end simulation

initial_save_time = save_time_index;

next_time = T.particle_run_time + dte;
% if save time occures before end of current timestep then cut current
% timestep
if save_time_index <= length(T.save_times) && ...
        next_time > T.save_times(save_time_index)
    % cut current timestep
    dte = T.save_times(save_time_index) - T.particle_run_time;
    % update timestep save index
    save_time_index = save_time_index + 1;
    T.save_data = 1;
else
    T.save_data = 0;
end

% if end time occures before end of current timestep then cut current
% timestep
if next_time > T.end_time
    % cut current timestep
    dte = T.end_time - T.particle_run_time;
%     if initial_save_time == initial_save_time
        % update timestep save index
        save_time_index = save_time_index + 1;
%     end
    % set to save data
    T.save_data = 1;
    % terminate particle sim after velocity is calculated
    particle_active = 1;
end

T.particle_run_time = T.particle_run_time +dte;