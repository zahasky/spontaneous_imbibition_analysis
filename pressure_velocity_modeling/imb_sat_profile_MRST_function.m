function [Saturation, p, G]= imb_sat_profile_MRST_function...
    (Saturation, G)
% imb_sat_profile_MRST_function
% Christopher Zahasky
% 12/3/2018 updated October 2020 to work with MRST

% Set polyfit order for fitting water saturation profile
polyfit_order = 5;

vox_length_z = 0.125; % [cm]
length_along_core_cm = [1:length(Saturation.CT_wr)].*vox_length_z - (vox_length_z/2);

% Assigne file name and load CT scan 
filename=dir([num2str(Saturation.scan_num),'_*']);
load(filename.name)
CT_exp = squeeze(squeeze(SAVG));

% Calculate saturation
slice_sat = ((CT_exp -Saturation.CT_ar)./(Saturation.CT_wr - Saturation.CT_ar));
% Plot slice average saturation
plot(length_along_core_cm, slice_sat, 'color', 'r')

% Calculate maximum saturation
Saturation.max_sat = max(slice_sat)*0.99;

% find where front crosses max_sat threshold, this will be starting point
% of fitting
back_slice_sat = fliplr(slice_sat');
backward_ind = find(back_slice_sat > Saturation.max_sat, 1, 'first');
start_ind = length(slice_sat)-backward_ind+1;
start_x = length_along_core_cm(start_ind); % end of flat saturation
% find where front goes to zero
end_ind = find(slice_sat < 0.03, 1, 'first')-1;

% fit saturation profile
x = length_along_core_cm(start_ind:end_ind);
y = slice_sat(start_ind:end_ind)';
plot(x(:), y(:), '.k')
p = polyfit(x,y,polyfit_order);

% define cell centered grid in dimension 2
dx = G.size(2) /G.cartDims(2); %[m | ft]
G.d2 = linspace(dx/2, G.size(2), G.cartDims(2))*100;
Sw = polyval(p,G.d2);

% Set upper and lower bounds of analytical profile
upper_ind = Sw > Saturation.max_sat;
if isempty(upper_ind)
    upper_ind = Sw > max(Sw);
end
Sw(upper_ind) = Saturation.max_sat;
Sw(G.d2<start_x) = Saturation.max_sat;
% Set everything past the front equal to zero
ind = find(Sw<0);
Sw(ind:end) = eps*10;
% plot results
plot(G.d2, Sw, 'g')

xlabel('Distance from inlet [cm]')
ylabel('Water saturation')
title(['Slice-average water saturation profiles'])
box on

Saturation.Sw = Sw;
