%% Element Locations
% Cross Array Only
N = 32;

z = zeros(1,N);

x = [15 25.981 30 25.981 15 0 -15 -25.981 -30 -25.981 -15 0 30 60 100.582 51.962 103.923 207.846 -36 -60 -108.533 -51.962 -103.923 -197.696 151.901 76.845 -144.685 -110.618 -76.125 153.432 346.314 -151.987];
y = [25.981 15 0 -15 -25.981 -30 -25.981 -15 0 15 25.981 30 51.962 103.923 179.6 -30 -60 -120 -62.354 -103.923 -188.186 30 60 117.580 -87.7 -121.796 -5.315 148.248 -131.853 80.531 42.595 87.75];

p = [x; y; z];
p = p';

%% Load Data

data = array_data{1,1}(1:32,1:1000000).';
data = data.*(1/(10^(-175/20)));

%% Beamforming

FS = 1000;
elev = 0:5:180;
az = 0:5:360;
c = 1445;
window = hanning(8192*4);
NFFT = 256;
f_range = [40 60];
overlap = 0.5;
weighting = 'simi_xarray_hanning';

[beamform_output,t,t_end] = beamform_3D_3(data,p,FS,elev,az,c,f_range,NFFT,window,overlap,weighting);

[~,ind] = min(abs(t-t_end));
beamform_output = beamform_output(1:ind-1,:,:,:);
beamform_output_db = squeeze(20*log10(abs(beamform_output)));

%% Plotting
f = linspace(f_range(1),f_range(2),NFFT);
[~,ind1] = min(abs(f - 40));
[~,ind2] = min(abs(f - 60));

[~,ind] = min(abs(t-t_end));
beamform_output_t = squeeze(mean(beamform_output_db(1:ind-1,:,:),1)); 
beamform_80_100hz = squeeze(mean(beamform_output_t(:,ind1:ind2),2));

figure
plot(beamform_80_100hz,elev,'linewidth',2)
set(gca,'Fontsize',30);
xlabel('Power (dB re 1\muPa^2)')
ylabel('Elevation (Degrees)')
ylim([-90 90]);
grid on

%%

beamform_output_f = squeeze(mean(beamform_output_db(:,:,ind1:ind2),3)).'; 

% for n = 1:size(beamform_output_f,1)
%     beamform_output_f(n,1:ind) = (beamform_output_f(n,1:ind)-mean(beamform_output_f(n,1:ind)))/std(beamform_output_f(n,1:ind));
% end

figure
fig = pcolor(t(1:ind-1)/60,elev,beamform_output_f);
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Time (min)')
xlim([0 t(ind-1)/60]);
ylabel('Elevation (Degrees)')
title(['Tape 24 Hour 1; f = 80 - 100Hz'])
colorbar;
colormap 'jet';

%%

figure
fig = pcolor(f_range,elev,beamform_output_t);
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Freqnuency (Hz)')
ylabel('Elevation (Degrees)')
title(['Tape 23 Hour 1; f = 80 - 100Hz'])
colorbar;
colormap 'jet';

%% Plotting
f = linspace(f_range(1),f_range(2),NFFT);

for i = 1;%1:length(elev)
beamform_output_elev = squeeze(beamform_output_db(:,i,:)).';

figure
fig = pcolor(t(1:ind-1)./60,f,beamform_output_elev);
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Time (min)')
xlim([0 t_end/60]);
ylabel('Frequency (Hz)')
title(['Elevation = ',num2str(elev(i))])
colorbar;
colormap 'jet';
end