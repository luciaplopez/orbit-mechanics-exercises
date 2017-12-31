%% Main

clear
clf
clc

%% Parameters

a = 7192 * 1e3; % (m) Semi mayor axis 
e = 0.004;  % (-) Eccentricity
i = 98.3;  % Inclination [deg]
raan = 257.7;  % RAAN [deg]
omega = 144.2;  % Argument of the perigee [deg]
t_0 = 0;  % Perigee passing time

omega_E = 2*pi/86164;  % Earth rotation period [rad/s]
period = 48*60*60;  % Orbit period [s]
u = 3.986004418e+14;  % Geocentric gravitational constant: u = GM
R = 6.371e+6; % Earth radius [m]
J_2 = 0.00108263;

%% Task 1:

dt = 60 * 1;
t = 0:dt:period;

% Load kep2orb
[r, v, ~, ~] = kep2orb(a, e, t_0, t);

% From polar to cartesian coordinates
x = r .* cos(v);
y = r .* sin(v);

% Load kep2cart
[rr, dotrr] = kep2cart(a, e, i, raan, omega, t_0,  t);

% Load cart2efix
[rrr, dotrrr] = cart2efix(rr, dotrr, t);


% Plot 
plot3(rrr(1), rrr(2), rrr(3))
hold on
grid on

legend('Sentinel-3')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Trayectory of the Sentinel-3 in Earth-Fixed System')


Earth_coast(3)

%% Task 2: write program yprime

% functions ode23, ode45 and ode113
