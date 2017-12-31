function [] = Earth_coast(dim)

% Draw the coast line of the earth in 2D or 3D
% dim = 2 for 2D
% dim = 3 for 3D


RE = 6371*1000;
lat = load('coast_lat.dat');
long = load('coast_long.dat');

if dim == 2
    plot(long,lat)
    xlim([-180 180])
    ylim([-90 90])
    set(gca,'Xtick',-180:30:180)
    set(gca,'Ytick',-90:30:90)
    
elseif dim == 3
    lat = lat*pi/180;
    long = long*pi/180;
      
    r = zeros(length(long(:)),3);
    for k=1:length(long(:))
        r(k,:)=[cos(lat(k))*cos(long(k)) cos(lat(k))*sin(long(k)) sin(lat(k))*(1-0.5*1/298.257223563)]*RE; 
    end
    
    [x,y,z] = sphere;
    mesh( RE*x, RE*y, RE*z, 'EdgeColor', 'None' );
    hold on;
    plot3(r(:,1),r(:,2),r(:,3));
    axis equal;
    hold on;
    %plot3([0:100000:50e6],zeros(1,501),zeros(1,501))

else
    error('Input not vaild!')
end