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
%J_2 = 0;
J_2 = 0.00108263;
ae = 6378 * 1e3;

% Compute initial conditions (t=0)
t = 0;
[rr, dotrr] = kep2cart(a, e, i, raan, omega, t_0,  t);

%% Task 2: Undisturbed Keppler problem

% Time computation interval
times = [0 period];

% Initial conditions
y0 = [rr(:,1); dotrr(:,1)];

% ODE integration options
options = odeset('InitialStep',5,'MaxStep',5);

% ODE integration
[t, y] = ode45(@(t,y) yprime_u(t,y), times, y0,options);

% Analitical solution (Task 1)
[rr, dotrr] = kep2cart(a, e, i, raan, omega, t_0, t);
[rrr, dotrrr] = cart2efix(rr, dotrr, t);

% Plot the analytical solution
figure(1)
plot3(rr(1, :), rr(2, :), rr(3, :),'-', 'LineWidth', 1.2)
hold on
grid on

legend('Sentinel-3')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Trayectory of the Sentinel-3 in Space-Fixed System')

Earth_coast(3)
 
%% Task 3: Compare the analitycal and numerical

% Plot analytical vs. numerical undisturbed ODE45
figure(2)

subplot(2,1,1) 
hold on 
grid on
plot(t, (y(:, 1) - rr(1, :)') ./ 1000)
plot(t, (y(:, 2) - rr(2, :)') ./ 1000)
plot(t, (y(:, 3) - rr(3, :)') ./ 1000)

legend('x','y','z')
xlabel('Time sequence (s)')
ylabel('Difference (km)')
title('Position with ODE45')

subplot(2,1,2) 
hold on
grid on
plot(t, (y(:, 4) - dotrr(1, :)') ./ 1000)
plot(t, (y(:, 5) - dotrr(2, :)') ./ 1000)
plot(t, (y(:, 6) - dotrr(3, :)') ./ 1000)

legend('vx','vy','vz')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Velocity with ODE45')


%% Task 4: Disturbed forces

% Disturbed case
options = odeset('InitialStep',5,'MaxStep',5);
[tt, yy] = ode23(@(tt,yy) yprime_d(tt,yy), times, y0, options);

% Plot disturbed vs undisturbed using ODE23
figure(4)
subplot(2,1,1) 
hold on
grid on
plot(tt, (yy(:,1) - y(:,1)) ./ 1000);
plot(tt, (yy(:,2) - y(:,2)) ./ 1000);
plot(tt, (yy(:,3) - y(:,3)) ./ 1000);

legend('x','y','z')
xlabel('Time sequence (s)')
ylabel('Difference (km)')
title('Position with ODE23')

subplot(2,1,2) 
hold on
grid on
plot(tt, (yy(:,4) - y(:,4)) ./ 1000);
plot(tt, (yy(:,5) - y(:,5)) ./ 1000);
plot(tt, (yy(:,6) - y(:,6)) ./ 1000);

legend('vx','vy','vz')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Velocity with ODE23')


[tt, yy] = ode45(@(tt,yy) yprime_d(tt,yy), times, y0, options);

% Plot disturbed vs undisturbed using ODE45
figure(5)
subplot(2,1,1) 
hold on
grid on
plot(tt, (yy(:,1) - y(:,1)) ./ 1000);
plot(tt, (yy(:,2) - y(:,2)) ./ 1000);
plot(tt, (yy(:,3) - y(:,3)) ./ 1000);

legend('x','y','z')
xlabel('Time sequence (s)')
ylabel('Difference (km)')
title('Position with ODE45')

subplot(2,1,2) 
hold on
grid on
plot(tt, (yy(:,4) - y(:,4)) ./ 1000);
plot(tt, (yy(:,5) - y(:,5)) ./ 1000);
plot(tt, (yy(:,6) - y(:,6)) ./ 1000);

legend('vx','vy','vz')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Velocity with ODE45')


[tt, yy] = ode23(@(tt,yy) yprime_d(tt,yy), times, y0, options);

% Compute the base unit vectors of the RSW system based on the 
% unperturbed orbit, i.e. starting from position and velocity vectors:
%  - The base unit vectors are different for each of the positions on the
%    orbit.
%  - The R (r_) vector is tangent to the radius vector.
%  - The S (s_) vector is the normal to the R on the orbiting plane
%    and pointing in the direction of the flight.
%  - The W (w_) vector is perpendicular to the other two axes and
%    perpendicular to the orbital plane. Parallel to the orbital 
%    angular momentum vector.
%  - We compute first R for each position on the orbit.
%  - Then we compute W using the vectorial product of two different
%    position vectors of the satellite.

% Undisturbed position and velocity
up = y(:, 1:3);
uv = y(:, 4:6);

% Compute the norm of each position vector
r_norm = sqrt(sum(up.^2,2));
% Normalize each postion vector in order to get an unitary vector
% pointing parallel to the radius vector
r_ = up ./ r_norm;
r_ = r_';


[n, m] = size(r_);

% To compute W we use the postion vector and the velocity vector,
% since both lay on the orbit plane
pv = up(1, :)';
vv = uv(1, :)';

% Compute vectorial product of the two vectors to obtain an orthogonal
% vector w.r.t. the orbit plane and make it unitary.
% W is the same for all the positions of the trajectory
ov = cross(pv, vv);
ov_norm = norm(ov, 2);
w_ = ov ./ ov_norm;

% Repeat w_ to have a matrix with the dimensions of r_
w_ = repmat(w_, [1 m]);

% Compute S as the vectorial product of W and R.
% S is different for all the positions of the trajectory
% (as R was).
% Compute the cross product of the column j of each matrix
s_ = cross(w_, r_);

% Decompose the difference vector between perturbed and unperturbed 
% orbit into this base.
delta_r = yy(:,1:3) - y(:,1:3);
delta_r = delta_r';

delta_v = yy(:,4:6) - y(:,4:6);
delta_v = delta_v';

% Project the difference vector onto the components of the RSW.
% dot treats the columns of A and B as vectors and calculates 
% the dot product of corresponding columns
delta_r_RSW = [dot(r_, delta_r); dot(s_, delta_r); dot(w_, delta_r)];

delta_v_RSW = [dot(r_, delta_v); dot(s_, delta_v); dot(w_, delta_v)];

figure(5)

subplot(1,2,1) 
hold on
grid on
plot(tt, delta_r_RSW(1, :) ./ 1000);
plot(tt, delta_r_RSW(2, :) ./ 1000);
plot(tt, delta_r_RSW(3, :) ./ 1000);
legend('R','S','W')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Position disturbed case 0DE45, RWS system')

subplot(1,2,2) 
hold on
grid on
plot(tt, delta_v_RSW(1, :) ./ 1000);
plot(tt, delta_v_RSW(2, :) ./ 1000);
plot(tt, delta_v_RSW(3, :) ./ 1000);
legend('R','S','W')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Velocity disturbed case ODE45, RWS system')



%% Task 4

[t,ye] = eulerode('yprime_d', t, y0);

[t,yk] = rungekuta('yprime_d', t, y0);

% Plot
figure(6)
subplot(2,1,1) 
hold on 
grid on
plot(t, ye(:,1)-rr(1, :)');
plot(t, ye(:,2)-rr(2, :)');
plot(t, ye(:,3)-rr(3, :)');
legend('x','y','z')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Position in Euler integrator')

subplot(2,1,2) 
hold on 
grid on
plot(t, ye(:,4)-dotrr(1, :)');
plot(t, ye(:,5)-dotrr(2, :)');
plot(t, ye(:,6)-dotrr(3, :)');
legend('vx','vy','vz')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Velocity in Euler integrator')



figure(7)
subplot(2,1,1) 
hold on 
grid on
plot(t, yk(:, 1)-rr(1, :)');
plot(t, yk(:, 2)-rr(2, :)');
plot(t, yk(:, 3)-rr(3, :)');
legend('x','y','z')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Position in Runge-Kutta integrator')


subplot(2,1,2) 
hold on 
grid on
plot(t, yk(:,4)-dotrr(1, :)');
plot(t, yk(:,5)-dotrr(2, :)');
plot(t, yk(:,6)-dotrr(3, :)');
legend('vx','vy','vz')
xlabel('Time sequence (s)')
ylabel('Difference (km/s)')
title('Velocity in Runge-Kutta integrator')
