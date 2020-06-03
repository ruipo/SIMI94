
%% Element Locations
% Cross Array Only
N = 32;

z = zeros(1,N);

x = [15 25.981 30 25.981 15 0 -15 -25.981 -30 -25.981 -15 0 30 60 100.582 51.962 103.923 207.846 -36 -60 -108.533 -51.962 -103.923 -197.696 151.901 76.845 -144.685 -110.618 -76.125 153.432 346.314 -151.987];
y = [25.981 15 0 -15 -25.981 -30 -25.981 -15 0 15 25.981 30 51.962 103.923 179.6 -30 -60 -120 -62.354 -103.923 -188.186 30 60 117.580 -87.7 -121.796 -5.315 148.248 -131.853 80.531 42.595 87.75];

p = [x; y; z];

%% Xarray without bad elements
N = 32;

z = zeros(1,N);

x = [15 25.981 30 25.981 15 0 -15 -25.981 -30 -25.981 -15 0 30 60 100.582 51.962 103.923 NaN -36 NaN -108.533 -51.962 -103.923 -197.696 151.901 NaN -144.685 NaN NaN 153.432 346.314 -151.987];
y = [25.981 15 0 -15 -25.981 -30 -25.981 -15 0 15 25.981 30 51.962 103.923 179.6 -30 -60 NaN -62.354 NaN -188.186 30 60 117.580 -87.7 NaN -5.315 NaN NaN 80.531 42.595 87.75];

p = [x; y; z];

%% Plot array Geometry
plot3(x/100,y/100,z','o','MarkerFaceColor', 'b');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('SIMI94 Array Geometry');
set(gca,'Fontsize',30);
grid on

%% Plot Beampattern

phi_vec = 0:pi/100:2*pi;
theta_vec = -pi/2:pi/100:pi/2;

phis = deg2rad(60);
thetas = deg2rad(-90);
weight = 'simi_xarray_hanning';

c = 1430;
f = 50;
lambda = c/f;
wls = 0;

[B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls,f);

[theta,phi] = meshgrid(theta_vec,phi_vec);

B_norm = abs(B)./max(max(abs(B)));
[X,Y,Z]=sph2cart(phi+pi/2,theta+pi/2,B_norm);

figure
h = surf(X,Y,Z,20*log10(B_norm));
set(h,'edgecolor','none');
caxis([-20 0]);
colorbar
colormap jet
xlabel('X')
ylabel('Y')
xlim([-1 1])
ylim([-1 1])
title('Cross Array Beampattern; f = 50 Hz; \phi = 150 deg.');
view(2)
set(gca,'Fontsize',30);