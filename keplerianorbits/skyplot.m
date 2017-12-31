function hsky = skyplot(azim,elev,line_style)

% SKYPLOT Polar coordinate plot using azimuth and elevation data.
%	SKYPLOT(AZIM,ELEV) makes a polar plot of the azimuth
%       AZIM [deg] versus the elevation ELEV [deg].
%       Negative elevation is allowed.
%       Azimuth is counted clock-wise from the North.
%	SKYPLOT(AZIM,ELEV,S) uses the linestyle, specified
%       in the string S (Default: '*').
%	See PLOT for a description of legal linestyles.
%
%	See also POLAR

% This function is modelled after polar.m. 
% Nico Sneeuw                           Munich                       05/01/96

% Checks and stuff
if isstr(azim) | isstr(elev)
	error('Input arguments must be numeric.');
end
if any(size(azim) ~= size(elev))
	error('AZIM and ELEV must be the same size.');
end
if any(elev > 90  | elev < -90)
   error('ELEV within [-90;90] [deg].')
end
if nargin < 2
   error('Requires 2 or 3 input arguments.')
elseif nargin == 2 
   line_style = '*';
end

% get hold state
cax = newplot;
next = lower(get(cax,'NextPlot'));
hold_state = ishold;

% get x-axis text color so grid is in same color
tc = get(cax,'xcolor');

% only do grids if hold is off
if ~hold_state

% make a radial grid
   hold on;
% check radial limits and ticks
   zenmax = max(90-elev(:));
   zenmax = 15*ceil(zenmax/15);
   elmax  = 90;
% define a circle
   az     = 0:pi/50:2*pi;
   xunit  = sin(az);
   yunit  = cos(az);

% make solid circles each 30 deg, and annotate.
% The horizon (elev=0) is made thicker.
% Inbetween, and below horizon only dotted lines.
   for i=[30 60]
      plot(xunit*i,yunit*i,'-','color',tc,'linewidth',1);
   end
   i=90; plot(xunit*i,yunit*i,'-','color',tc,'linewidth',2);
   for i=[15:30:75 105:15:zenmax]
      plot(xunit*i,yunit*i,':','color',tc,'linewidth',1);
   end
   for i=30:30:zenmax
      text(0,i,[' ' num2str(90-i)],'verticalalignment','bottom');
   end

% plot spokes
   az = (1:6)*2*pi/12;                     % define half circle
   caz = cos(az); 
   saz = sin(az);
   ca = [-caz; caz];
   sa = [-saz; saz];
   plot(elmax*ca,elmax*sa,'-','color',tc,'linewidth',1);
   if zenmax > elmax
      plot(zenmax*ca,zenmax*sa,':','color',tc,'linewidth',1);
   end

% annotate spokes in degree
   rt = 1.1*elmax;
   for i = 1:length(az)
      loc1 = int2str(i*30);
      if i == length(az)
         loc2 = int2str(0);
      else
         loc2 = int2str(180+i*30);
      end
      text( rt*saz(i), rt*caz(i),loc1,'horizontalalignment','center');
      text(-rt*saz(i),-rt*caz(i),loc2,'horizontalalignment','center');
   end

% brush up axes
   view(0,90);
   axis(max(zenmax,elmax)*[-1 1 -1.1 1.1]);
   set(cax,'position',[.05 .05 .9 .9])
end


% transform data to Cartesian coordinates.
yy = (90-elev).*cos(azim/180*pi);
xx = (90-elev).*sin(azim/180*pi);

% plot data on top of grid 
q = plot(xx,yy,line_style);
if nargout > 0, hsky = q;end

if ~hold_state
   axis('equal')
   axis('off')
end

% reset hold state
if ~hold_state, set(cax,'NextPlot',next); end
