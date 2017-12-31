%% Main

clear
clf
clc

addpath('../keplerianorbits');

%% Parameters

a = 7192 * 1e3; % (m) Semi mayor axis 
e = 0.004;  % (-) Eccentricity
i = 98.3;  % Inclination [deg]
raan = 257.7;  % RAAN [deg]
omega = 144.2;  % Argument of the perigee [deg]
t_0 = 0;  % Perigee passing time

omega_E = 2*pi/86164;  % Earth rotation period [rad/s]
u = 3.986004418e+14;  % Geocentric gravitational constant: u = GM
period = 3*2*pi*sqrt(a^3/u);  % Orbit period [s]
R = 6.371e+6; % Earth radius [m]
J_2 = 0.00108263;

% Compute initial conditions (t=0)

t = 0;

% Load kep2cart
[rr, dotrr] = kep2cart(a, e, i, raan, omega, t_0,  t);

%% Task 2: Write program yprime

% Time computation interval
times = [0 period];

% Initial conditions
y0 = [rr(:,1); dotrr(:,1)];

% ODE integration options
options = odeset('InitialStep',5,'MaxStep',5);

% ODE integration
[t, y] = ode23(@(t,y) yprime(t,y,u), times, y0, options);

% Analitical solution
[rr, dotrr] = kep2cart(a, e, i, raan, omega, t_0, t);
[rrr, dotrrr] = cart2efix(rr, dotrr, t);

% Plot 
figure(1)
plot3(rrr(1, :), rrr(2, :), rrr(3, :))
hold on
grid on

legend('Sentinel-3')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Trayectory of the Sentinel-3 in Earth-Fixed System')

figure(2)
hold on
plot(t, y(:, 1) - rr(1, :)')
plot(t, y(:, 2) - rr(2, :)')
plot(t, y(:, 3) - rr(3, :)')

legend('')
xlabel('Time sequence (s)')
ylabel('Difference (m)')
title('Position with ODE23')