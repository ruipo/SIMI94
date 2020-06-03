
c = 1435;
FS = 1000;
win_len = 4096;
overlap = 0.5;

data = array_data{1,1}(1:24,3600000*2+1:3600000*3).';
data = data.*(1/(10^(-165/20)));
% data(:,32) = [];
% data(:,29) = [];
% data(:,28) = [];
% data(:,26) = [];
% data(:,20) = [];
% data(:,18) = [];
% 
% x = [15 25.981 30 25.981 15 0 -15 -25.981 -30 -25.981 -15 0 30 60 100.582 51.962 103.923 -36 -108.533 -51.962 -103.923 -197.696 151.901 -144.685 153.432 346.314];
% y = [25.981 15 0 -15 -25.981 -30 -25.981 -15 0 15 25.981 30 51.962 103.923 179.6 -30 -60 -62.354 -188.186 30 60 117.580 -87.7 -5.315 80.531 42.595];

x = [15 25.981 30 25.981 15 0 -15 -25.981 -30 -25.981 -15 0 30 60 100.582 51.962 103.923 207.846 -36 -60 -108.533 -51.962 -103.923 -197.696];
y = [25.981 15 0 -15 -25.981 -30 -25.981 -15 0 15 25.981 30 51.962 103.923 179.6 -30 -60 -120 -62.354 -103.923 -188.186 30 60 117.580];


p = [x; y];

daz = 2*pi/60;
azlist = 0:daz:2*pi-daz;
%%

[time_beamform,t] = time_beamform_2D(data,p,FS,azlist,c,win_len,overlap);

%%
figure
fig = pcolor(t/60,azlist*180./pi,20*log10(time_beamform));
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Time (min)')
ylabel('Azimuth (Degrees)')
title(['Tape 31 Seg.2 XArray Beamform'])
colorbar;
%caxis([120 130]);
colormap 'jet';
%%
figure
for i = 1:45
t_average = 20*log10(time_beamform(:,i));
polarplot([azlist 2*pi],[t_average; t_average(1)]);
set(gca,'Fontsize',30);
rlim([120 130]);
pause
end
%%
figure
t_average = 20*log10(mean(time_beamform,2));
polarplot([azlist 2*pi],[t_average; t_average(1)],'linewidth',3);
set(gca,'Fontsize',30);
rlim([115 125]);
title(['Tape 31 Seg.2 XArray Beamform'])
