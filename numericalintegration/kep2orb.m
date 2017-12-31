function [r, v, M, E] = kep2orb(a, e, t_0, t)
% Description: compute polar coordinates r (radius) and v(true anomaly)
% based on input orbital elements
%
% Inputs: semi mayor axis, eccentricity, perigee passign time and time
%
% Outputs: radius, true anomally, mean motion and eccentric anomaly
%

%% Parameters

% Geocentric gravitational constant: u = GM
u = 3.986004418e+14;  

% Mean motion
n = sqrt(u/(a^3));  % Kepler's 3rd law: u = (n^2*a^3) 

% Tolerance for eccentric anomaly difference: diff = (E_1 - E_0)
tol = 1e-6;

nsteps = length(t);

% Allocate arrays for radius and true anomally.
% One for each time step
r = zeros(nsteps, 1);
v = zeros(nsteps, 1);
M = zeros(nsteps, 1);
E = zeros(nsteps, 1);


%% Time iteration

% For each time step compute the position of the satellite
for i=1:nsteps
    
    time = t(i);
    
    if (time - t_0) < 0
        % Compute true anomaly for current time [rad]
         M_1 = n * (time - t_0) + 2*pi;
    else
        M_1 = n * (time - t_0);
    end
    
    M_1 = mod(M_1, 2*pi);
        
    M(i, 1) = M_1;
    
    % Initialize eccentric anomaly [rad]
    E_0 = M_1;
    
    % Initialize difference [rad]
    diff = 1e10;
    
    % Newton iterarion
    while diff >= tol
        % Eccentric anomaly [rad]
        del_E = ((M_1 + e * sin(E_0)) - E_0) / (1 - (e * cos(E_0))); 
        E_1 = E_0 + del_E;
        diff = E_1 - E_0;
        E_0 = E_1;
    end
    
    E(i, 1) = E_1;
    
    % True anomaly [rad]
    v(i, 1) = 2 * atan2(sqrt(1 + e) * sin(E_1/2), sqrt(1 - e) * cos(E_1/2));
    
    % Radius of satellite's orbit [m]
    r(i, 1) = a * (1 - e*cos(E_1));  
end

% Outputs r and v have being assigned during time fordwarding for loop

end 



