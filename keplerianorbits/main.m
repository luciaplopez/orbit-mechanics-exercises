%% Main

clear
clf
clc

a = xlsread('data.xlsx','A1:A5');  % Semi mayor axis
e = xlsread('data.xlsx','B1:B5');  % Eccentricity
i = xlsread('data.xlsx','C1:C5');  % Inclination
raan = xlsread('data.xlsx','D1:D5');  % RAAN
omega = xlsread('data.xlsx','E1:E5');  % Argument of the perigee
t_0 = xlsread('data.xlsx','F1:F5');  % Perigee passing time

omega_E = 2*pi/86164;  % Earth rotation period [rad/s]
period = 24*60*60;  % Orbit period [s]
u = 3.986004418e+14;  % Geocentric gravitational constant: u = GM
R = 6.371e+6; % Earth radius [m]

%% Task 1:
dt = 60 * 1;
t = 0:dt:period;

[r1, v1, M1, E1] = kep2orb(a(1), e(1), t_0(1), t);
[r2, v2, M2, E2] = kep2orb(a(2), e(2), t_0(2), t);
[r3, v3, M3, E3] = kep2orb(a(3), e(3), t_0(3), t);
[r4, v4, M4, E4] = kep2orb(a(4), e(4), t_0(4), t);
[r5, v5, M5, E5] = kep2orb(a(5), e(5), t_0(5), t);

% Plot for the 5 satellites in the orbital plane for one orbital
% revolution


figure(1)
hold on 
grid on 

x1 = r1 .* cos(v1);
y1 = r1 .* sin(v1);

x2 = r2 .* cos(v2);
y2 = r2 .* sin(v2);

x3 = r3 .* cos(v3);
y3 = r3 .* sin(v3);

x4 = r4 .* cos(v4);
y4 = r4 .* sin(v4);

x5 = r5 .* cos(v5);
y5 = r5 .* sin(v5);

plot(x1, y1)
plot(x2, y2)
plot(x3, y3)
plot(x4, y4)
plot(x5, y5)

legend('GOCE', 'GPS', 'MOLNIYA', 'GEO', 'MICHIBIKI')
xlabel('x(m)')
ylabel('y(m)')
title('Orbit of the satellites in a 2D plane')

% Plot M, E, v, v-M for the GPS satellite
figure(2)
hold on
grid on

t_plot = t;
t_plot = t_plot ./ 3600;

plot(t_plot, rad2deg(M2))
plot(t_plot, rad2deg(E2))
plot(t_plot, rad2deg(v2))
plot(t_plot, rad2deg(v2-M2))

legend('M','E','v','v-M')
xlabel('time(s)')
ylabel('angle(degre)')
title('M, E, v, v-M of the GPS satellite')

% Plot M, E, v, v-M for the MOLNIYA satellite

figure(3)
hold on
grid on

t_plot = t;
t_plot = t_plot ./ 3600;

plot(t_plot, rad2deg(M3))
plot(t_plot, rad2deg(E3))
plot(t_plot, rad2deg(v3))
plot(t_plot, rad2deg(v3-M3))

legend('M','E','v','v-M')
xlabel('time(s)')
ylabel('angle(degre)')
title('M, E, v, v-M of the MOLNIYA satellite')

%% Task 2:
dt = 60 * 1;
t = 0:dt:period;
nsteps = length(t);

rr = zeros(3, nsteps, 5);
dotrr = zeros(3, nsteps, 5);



for j=1:5
    [rr(:, :, j), dotrr(:, :, j)] = kep2cart(a(j), e(j), i(j), raan(j), omega(j), t_0(j),  t);
end


figure(4)


for j=1:5
    plot3(rr(1, :, j), rr(2, :, j), rr(3, :, j))
    hold on
end

legend('GOCE', 'GPS', 'MOLNIYA', 'GEO', 'MICHIBIKI')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Orbits of 5 satellites in 3D plane')
grid on


Earth_coast(3)


figure(5)


for j=1:5
    plot(rr(1, :, j), rr(2, :, j))
    hold on
end

legend('GOCE', 'GPS', 'MOLNIYA', 'GEO', 'MICHIBIKI')
xlabel('x(m)')
ylabel('y(m)')
title('Orbits of 5 satellites in XY plane')
grid on



figure(6)

for j=1:5
    plot(rr(1, :, j), rr(3, :, j))
    hold on
end

legend('GOCE', 'GPS', 'MOLNIYA', 'GEO', 'MICHIBIKI')
xlabel('x(m)')
ylabel('z(m)')
title('Orbits of 5 satellites in XZ plane')
grid on


figure(7)

for j=1:5
    plot(rr(2, :, j), rr(3, :, j))
    hold on
end

legend('GOCE', 'GPS', 'MOLNIYA', 'GEO', 'MICHIBIKI')
xlabel('y(m)')
ylabel('z(m)')
title('Orbits of 5 satellites in YZ plane')
grid on


figure(8)

for j=1:5
    velocity = sqrt(dotrr(1, :, j).^2 + dotrr(2, :, j).^2 + dotrr(3, :, j).^2);
    plot(t, velocity)
    hold on
end

legend('GOCE', 'GPS', 'MOLNIYA', 'GEO', 'MICHIBIKI')
xlabel('time(s)')
ylabel('v(m/s)')
title('Time serires magnitude of velocity')
grid on

%% Task 3:
dt = 60 * 1;
t = 0:dt:period;
nsteps = length(t);

rrr = zeros(3, nsteps, 5);
dotrrr = zeros(3, nsteps, 5);


for j=1:5
    [rrr(:, :, j), dotrrr(:, :, j)] = cart2efix(rr(:, :, j), dotrr(:, :, j), t);
end

figure(9)

for j=[1 2 3 5]
    plot3(rrr(1, :, j), rrr(2, :, j), rrr(3, :, j))
    hold on
end

scatter3(rrr(1, 1, 4), rrr(2, 1, 4), rrr(3, 1, 4),'filled')
legend('GOCE', 'GPS', 'MOLNIYA', 'MICHIBIKI','GEO')
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
title('Trayectory of 5 satellites in Earth-Fixed System')
grid on

Earth_coast(3)

latitude = zeros(nsteps, 5);
longitude = zeros(nsteps, 5);

for j=1:5
    longitude(:, j) = atan2(rrr(2, :, j), rrr(1, :, j));
    latitude(:, j) = atan2(rrr(3, :, j), ...
        sqrt(rrr(1, :, j).^2 + rrr(2, :, j).^ 2));
end



figure(10)

for j=1:5
    scatter(rad2deg(longitude(:, j)), rad2deg(latitude(:, j)))
    hold on
end

legend('GOCE', 'GPS', 'MOLNIYA','GEO','MICHIBIKI')
xlabel('longitude')
ylabel('latitude')
title('Satellite ground-tracks on the Earth surface')
grid on

Earth_coast(2)

%% Task 4:
dt = 60 * 1;
t = 0:dt:period;
nsteps = length(t);

rrrr = zeros(3, nsteps, 5);
azim = zeros(nsteps, 5);
elev = zeros(nsteps, 5);

figure(11)

for j=1:5
    [rrrr(:, :, j), azim(:, j), elev(:, j)] = efix2topo(rrr(:, :, j), t);
    skyplot(rad2deg(azim(:, j)), rad2deg(elev(:, j)), '+');
    hold on
end

title('Orbits in the topocentric system')
legend('GOCE', 'GPS', 'MOLNIYA','GEO','MICHIBIKI')
visibility = zeros(nsteps, 5);

for j=1:5
    visibility(:, j) = rad2deg(elev(:, j)) > 0.0;
end



figure(12)


subplot(5,1,1) 
bar(t./3600, visibility(:, 1))
xticks(0:1:24)
yticks(0:1)
title('GOCE')

subplot(5,1,2) 
bar(t./3600, visibility(:, 2))
xticks(0:1:24)
yticks(0:1)
title('GPS')

subplot(5,1,3) 
bar(t./3600, visibility(:, 3))
xticks(0:1:24)
yticks(0:1)
title('MOLNIYA')

subplot(5,1,4) 
bar(t./3600, visibility(:, 4))
xticks(0:1:24)
yticks(0:1)
title('GEO')

subplot(5,1,5) 
bar(t./3600, visibility(:, 5))
xticks(0:1:24)
yticks(0:1)
title('MICHIBIKI')











