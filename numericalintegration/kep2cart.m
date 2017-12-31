function [rr, dotrr] = kep2cart(a, e, i, raan, omega, t_0,  t)
% Description: this function transforms the keplerian elements to position
% and velocity in an inertial (space-fixed) system
%
% Inputs: semi mayor axis, eccentricity, inclination, right ascension of
% the asecendent node, argument of the perigee, perigee passing time and
% time
%
% Outputs: the position and velocity in and inertial space-fixed system
%

raan = deg2rad(raan);
omega = deg2rad(omega);
i = deg2rad(i);

[r, v, ~, ~] = kep2orb(a, e, t_0, t);

nsteps = length(t);

rr = zeros(3, nsteps);
dotrr = zeros(3, nsteps);

rr(1, :) = r(:) .* cos(v(:));
rr(2, :) = r(:) .* sin(v(:));

u = 3.986004418e+14;  % Geocentric gravitational constant: u = GM
k = sqrt(u ./ (a .* (1 - e^2)));

dotrr(1, :) = -k .*  sin(v(:));
dotrr(2, :) = k .* (e + cos(v(:)));
 
%% Rotations

% First rotation of RAAN in z axis
R3_raan = [cos(raan) sin(raan) 0;
    -sin(raan) cos(raan) 0;
    0 0 1];

% Second rotation of the inclination in x axis
R1_i = [1 0 0
    0 cos(i) sin(i)
    0 -sin(i) cos(i)];

% Third rotation of the argument of the perigee in z axis
R3_omega = [cos(omega) sin(omega) 0; 
    -sin(omega) cos(omega) 0;
    0 0 1];

% ACHTUNG: R(-omega) = R(omega).T

% Position rotation
rr = R3_raan' * R1_i' * R3_omega' * rr;

% Velocity rotation
dotrr = R3_raan' * R1_i' * R3_omega' * dotrr;

end 

