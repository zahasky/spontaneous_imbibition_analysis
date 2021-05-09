% streamtube_moment_calc_function
function [M0, Xc, Sx]= streamtube_moment_calc_function(PET_4D_cc, ...
    vox_length, varargin)
% Christopher Zahasky
% 10/5/2016
% This script is used to calculate the advection dispersion pulse or
% continous solution given the follow core data

PET_dim = size(PET_4D_cc);
% preallocate zero moment matrix
M0 = zeros(PET_dim(1), PET_dim(2), PET_dim(4));
% preallocate first moment matrix
Xc = zeros(PET_dim(1), PET_dim(2), PET_dim(4));
% preallocate second moment matrix
Sx = zeros(PET_dim(1), PET_dim(2), PET_dim(4));

x = [0.5:PET_dim(3)]'.*vox_length;
for t_frame=1:PET_dim(4)
    for i=1:PET_dim(1)
        for j=1:PET_dim(2)
            % isolate streamtube
            streamtube_n = squeeze(squeeze(PET_4D_cc(i,j,:,t_frame)));
            % check to make sure tracer is in tube
            if sum(streamtube_n) > 0
                % calculate zero moment of streamtube
                m0 = trapz(x, streamtube_n); % 'x' is the voxel center locations, 'streamtube_n' is the concentration profile
%                 m0 = sum(streamtube_n);
                M0(i,j,t_frame) = m0;
                
                % calculate first moment of streamtube
                m1 = trapz(x, streamtube_n.*x);  % 'x' is the voxel center locations, 'streamtube_n' is the concentration profile
                % calculate center of mass of streamtube
                xc = m1/m0;
                Xc(i,j,t_frame) = xc;
                
                % calculate second moment of streamtube
                m2 = trapz(x, streamtube_n.*x.^2);
                % calculate spread of tracer
                sxx = m2/m0 - xc^2;
                Sx(i,j,t_frame) = sxx;
            end
        end
    end
end

% if addition variables are listed this triggers the plot of percent change
% in streamtube moments
extra_var = nargin-2;
if extra_var > 0 && length(varargin{1})>1
    % plot of 3D center of mass with time
    figure
    slice = (M0(:,:,varargin{1}(end))-M0(:,:,varargin{1}(1)))./M0(:,:,varargin{1}(end)).*100;
    h3 = imagesc(slice);
    set(h3,'alphadata',~isnan(slice))
    title(['Percent change in zero moment'])
    axis equal
    axis tight
    axis off
    % caxis([-10 10])
    colorbar
end

