function [rrr, dotrrr] = cart2efix(rr, dotrr, t)
% Description: transform position and velocity in a spacefixed system into
% position and velocity in an Earth-fixed system
%
% Inputs: position and velocity in an inertial spacefixed system and time
%
% Outputs: the position and velocity in an Earth-fixed system


%% Parameters: 

% Earth rotation period [rad/s]
omega_E = 2*pi/86164;

% Sideral angle
sideral = (3 + 29/60) * 15;

nsteps = length(t);

rrr = zeros(3, nsteps);
dotrrr = zeros(3, nsteps);

%% Time loop:

for j=1:nsteps

    theta_0 = omega_E .* t(j) + deg2rad(sideral);

    R3_theta_0 = [cos(theta_0) sin(theta_0) 0;
        -sin(theta_0) cos(theta_0) 0; 
        0 0 1];

    rrr(:, j) = R3_theta_0 * rr(:, j);
    dotrrr(:, j) = R3_theta_0 * dotrr(:, j);
end

end 



