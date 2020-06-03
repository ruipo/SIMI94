
%% Element Locations
% 7m spacing in vert array
N = 32;
z2 = zeros(N,1);
d = 7;
for n = 0:N-1
    z2(n+1) = -(n-(N-1)/2)*d; 
end

z1 = z2(1)*ones(length(z2),1);

x1 = [15 25.981 30 25.981 15 0 -15 -25.981 -30 -25.981 -15 0 30 60 100.582 51.962 103.923 207.846 -36 -60 -108.533 -51.962 -103.923 -197.696 151.901 76.845 -144.685 -110.618 -76.125 153.432 346.314 -151.987];
y1 = [25.981 15 0 -15 -25.981 -30 -25.981 -15 0 15 25.981 30 51.962 103.923 179.6 -30 -60 -120 -62.354 -103.923 -188.186 30 60 117.580 -87.7 -121.196 -5.315 148.248 -131.853 80.531 42.595 87.75];
x2 = zeros(1,length(z2));
y2 = zeros(1,length(z2));

x = [x1 x2];
y = [y1 y2];
z = [z1;z2];

p = [x; y; z'];

%% Plot array Geometry
plot3(x,y,z','o','MarkerFaceColor', 'b');
xlabel('x (m)');
ylabel('y (m)');
zlabel('z (m)');
title('SIMI94 Array Geometry');
set(gca,'Fontsize',30);
grid on

%% Plot Beampattern

phi_vec = 0:pi/100:2*pi;
theta_vec = 0:pi/100:pi;

phis = deg2rad(0);
thetas = deg2rad(90);
weight = 'uniform';

c = 1435;
f = 90;
lambda = c/f;
wls = 0;

[B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls);

[theta,phi] = meshgrid(theta_vec,phi_vec);

B_norm = abs(B)./max(max(abs(B)));
[X,Y,Z]=sph2cart(phi+pi/2,theta+pi/2,B_norm);

figure
surf(X,Y,Z,20*log10(B_norm));

%% Plot Vert Array Beampattern vs frequency
clear beampattern_mat
set(0,'DefaultFigureVisible','off')
N = 32;
z = zeros(N,1);
d = 7;
for n = 0:N-1
    z(n+1) = -(n-(N-1)/2)*d; 
end

p = [zeros(1,N) ; zeros(1,N) ; z'];
phi_vec = 0;
theta_vec = 0:pi/500:pi;
phis = 0;
thetas = pi/2;
weight = 'hanning';
wls = 'none';

f = 80;
c = 1430;

for i = 1:length(f) 
    disp(i)
    lambda = c/f(i);
    [B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls,f(i));
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    beampattern_mat(i) = getframe(gcf);
end

v_beampattern = VideoWriter('VertArray_Beampattern.avi','Uncompressed AVI');
v_beampattern.FrameRate = 3;
open(v_beampattern)
writeVideo(v_beampattern,beampattern_mat)
close(v_beampattern)

%% Plot Beampattern width
N = 32;
z = zeros(N,1);
d = 7;
for n = 0:N-1
    z(n+1) = -(n-(N-1)/2)*d; 
end

p = [zeros(1,N) ; zeros(1,N) ; z'];
phi_vec = 0;
theta_vec = 0:pi/500:pi;
phis = 0;
thetas = pi/2;
weight = 'hanning';
wls = 'none';

f = 90;
c = 1430;

for i = 1:length(f) 
    disp(i)
    lambda = c/f(i);
    [B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas,weight,lambda,wls,f(i));
    [~,loc] = findpeaks(-abs(B));
    angle = -1*(rad2deg(theta_vec)-90);
    half_beamwidth(i) = abs(angle(end)-angle(loc(end)));
    clear loc
end
close all
figure
plot(f,half_beamwidth)
hold on
plot(105*ones(2,1),[20,90]);
plot(80*ones(2,1),[20,90]);
plot(100*ones(2,1),[20,90]);

%% Plot Vert Array Beampattern vs Elevation
clear beampattern_mat
set(0,'DefaultFigureVisible','off')
N = 32;
z = zeros(N,1);
d = 7;
for n = 0:N-1
    z(n+1) = -(n-(N-1)/2)*d; 
end

p = [zeros(1,N) ; zeros(1,N) ; z'];
phi_vec = 0;
theta_vec = 0:pi/700:pi;
phis = 0;
thetas = 0:pi/36:pi;
weight = 'uniform';
wls = 'none';

f = 50;
c = 1430;
lambda = c/f;

for i = 1:length(thetas)
    disp(i)
    [B,vs] = beampattern_plot(p,phi_vec,theta_vec,phis,thetas(i),weight,lambda,wls,f);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    beampattern_mat(i) = getframe(gcf);
end

v_beampattern = VideoWriter('VertArray_Beampattern2.avi','Uncompressed AVI');
v_beampattern.FrameRate = 3;
open(v_beampattern)
writeVideo(v_beampattern,beampattern_mat)
close(v_beampattern)