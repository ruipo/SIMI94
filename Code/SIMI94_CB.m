%% Element Locations
N = 32;
% 0.75m spacing in the middle
z = zeros(N,1);
d = 7;
for n = 0:N-1
    z(n+1) = -(n-(N-1)/2)*d; 
end

p = [zeros(1,N) ; zeros(1,N) ; fliplr(z')];
p = p';

%% Load Data

data = array_data{1,1}(33:64,1:6*3600000).';
data = (data.*(1/(10^(-175/20))));

%% Beamforming

FS = 1000;
elev = -90:1:90;
az = 0;
c = 1445;
window = hanning(4096);
NFFT = 512;
f_range = [0 100];
overlap = 0.5;
weighting = 'hanning';

[beamform_output,t,t_end] = beamform_3D(data,p,FS,elev,az,c,f_range,NFFT,window,overlap,weighting);

[~,ind] = min(abs(t-t_end));
beamform_output = beamform_output(1:ind-1,:,:,:);
beamform_output_db = squeeze(10*log10(abs(beamform_output)./max(max(max(abs(beamform_output))))));

%% Plotting
f = linspace(f_range(1),f_range(2),NFFT);
[~,ind1] = min(abs(f - 80));
[~,ind2] = min(abs(f - 100));

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

for i = 91;%1:length(elev)
beamform_output_elev = squeeze(beamform_output_db(:,i,:)).';

figure
fig = pcolor(t(1:ind-1)./3600,f,beamform_output_elev);
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Time (min)')
xlim([0 t_end/3600]);
ylabel('Frequency (Hz)')
title(['Elevation = ',num2str(elev(i))])
colorbar;
colormap 'jet';
end