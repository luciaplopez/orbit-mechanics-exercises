function [rrrr, azim, elev] = efix2topo(rrr, t)
% Description: transforms position and velocity into a topocentric system
% centered at the Wettzell
%
% Inputs: position and velocity in an Earth-fixed system
%
% Outputs: position and velocity 
%

%% Parameters

% Postion vector in an Earth-fixed system given by:
r_w = [4075.53022; 931.78130; 4801.61819]*1e3;  % [m]

%% Translated vector
r_tran = rrr - r_w;

% Topocentric vector

% Latitude and longitude of Wettzell
lon_w = atan2(r_w(2), r_w(1));
lat_w = atan2(r_w(3), sqrt(r_w(1).^2 + r_w(2).^ 2));

% Right handed to left handed system converter matrix
Q = [-1 0 0; 
    0 1 0; 
    0 0 1];

R2ang = pi/2 - lat_w;
R3ang = lon_w;

R2_lat = [cos(R2ang) 0 -sin(R2ang);
          0 1 0;
          sin(R2ang) 0 cos(R2ang)];
            
R3_lon = [cos(R3ang) sin(R3ang) 0;
          -sin(R3ang) cos(R3ang) 0;
          0 0 1];

rrrr = Q * R2_lat * R3_lon * r_tran;


%% Azimuth and elevaton calculation

nsteps = length(t);

azim = zeros(1, nsteps);
elev = zeros(1, nsteps);


azim = atan2(rrrr(2, :), rrrr(1, :));
elev = atan2(rrrr(3, :), ...
        sqrt(rrrr(1, :).^2 + rrrr(2, :).^ 2));

end




