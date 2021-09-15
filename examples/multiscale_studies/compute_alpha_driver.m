clc
clear all

% pressure_angles = pi/2 - logspace(-5, pi, 10)/2;
% 
% for i = 1:length(pressure_angles)
%    [alpha(i), alpha_saffman(i), infiltration_angle(i)] = compute_alpha(pressure_angles(i));
% end

pressure_angle = [pi/2-1e-2, pi/2-1e-4];
%interface_offset = [-4+1e-4, -3+1e-4, -2+1e-4, -1+1e-4, 1e-4, 1e-3, 1e-2, 1e-1, 1+1e-4, 2+1e-4];
interface_offset = 1e-4;
volume_fraction = 0.4;
theta = pi/6;

alpha = zeros(size(pressure_angle));
alpha_saffman = zeros(size(alpha));
alpha_bar = zeros(size(alpha));
alpha_saffman_bar = zeros(size(alpha));
infiltration_angle = zeros(size(alpha));
interface = zeros(size(alpha));
U = zeros(50,length(alpha));
K = zeros(2,2,length(alpha));

for i = 1:length(alpha)
   [alpha(i), alpha_saffman(i), alpha_bar(i), alpha_saffman_bar(i), U(:,i),...
       K(:,:,i), interface(i), infiltration_angle(i)] = compute_alpha(pressure_angle(i),...
            interface_offset, volume_fraction, theta);
end