

bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',20,'CutoffFrequency2',350,'SampleRate',1000);
data = array_data{1,3}(64,1:end);

data = filtfilt(bandfilt,data);

FS = 1000;
dt = 2048;
overlap = 0.5;

window_start = round(dt-dt*overlap);
num_window = round(length(data)/window_start)-2;

t = zeros(num_window,1);
ts = zeros(dt,num_window);
for l = 1:num_window
    ts(:,l) = data(l*window_start-window_start+1:l*window_start-window_start+dt);
    t(l) = ((l+1)*window_start-window_start+1)/FS;
end

meanlist = mean(ts);
varlist = var(ts);
kurtlist = kurtosis(ts);
skewlist = skewness(ts);


%%

figure
plot(t/60,kurtlist);
hold on
plot(t/60,3*ones(length(t),1),'r','linewidth',2)
xlabel('Time (min)')
ylabel('Kurtosis');
set(gca,'fontsize',25);
grid on

%%

for i = 1:32
    i
FS = 1000; 
window_size = 256;
window = hamming(window_size);
NFFT = 4*window_size;
overlap = 0.5;
df = FS/NFFT;
f = -FS/2:df:FS/2-df;
chn = 1;
data_name = 'Tape 23';

data = array_data{1,3}(65-i,1:2500000).';
data = data.*(1/(10^(-155/20)/1E-6));

dt = 1/FS;
t = 0:dt:length(data)/FS-dt;

[~,psd_med(i,:),freq] = psd(data,window_size,window,overlap,NFFT,FS,chn,data_name);

end

 %%
[~,ind1] = min(abs(freq-17.538));
[~,ind2] = min(abs(freq-22.098));
[~,ind3] = min(abs(freq-27.840));
[~,ind4] = min(abs(freq-35.077));
[~,ind5] = min(abs(freq-44.194));
[~,ind6] = min(abs(freq-55.681));
[~,ind7] = min(abs(freq-70.154));
[~,ind8] = min(abs(freq-88.388));
[~,ind9] = min(abs(freq-111.362));
[~,ind10] = min(abs(freq-140.308));
[~,ind11] = min(abs(freq-176.777));
[~,ind12] = min(abs(freq-222.725));
[~,ind13] = min(abs(freq-280.616));
[~,ind14] = min(abs(freq-350));

third_oband(1,:) = squeeze(mean(psd_med(:,ind1:ind2),2)); %20Hz
third_oband(2,:) = squeeze(mean(psd_med(:,ind2:ind3),2)); %25Hz  
third_oband(3,:) = squeeze(mean(psd_med(:,ind3:ind4),2)); %31Hz
third_oband(4,:) = squeeze(mean(psd_med(:,ind4:ind5),2)); %40Hz
third_oband(5,:) = squeeze(mean(psd_med(:,ind5:ind6),2)); %50Hz
third_oband(6,:) = squeeze(mean(psd_med(:,ind6:ind7),2)); %63Hz
third_oband(7,:) = squeeze(mean(psd_med(:,ind7:ind8),2)); %80Hz
third_oband(8,:) = squeeze(mean(psd_med(:,ind8:ind9),2)); %100Hz
third_oband(9,:) = squeeze(mean(psd_med(:,ind9:ind10),2)); %125Hz
third_oband(10,:) = squeeze(mean(psd_med(:,ind10:ind11),2)); %160Hz
third_oband(11,:) = squeeze(mean(psd_med(:,ind11:ind12),2)); %200Hz
third_oband(12,:) = squeeze(mean(psd_med(:,ind12:ind13),2)); %250Hz
third_oband(13,:) = squeeze(mean(psd_med(:,ind13:ind14),2)); %315Hz


%% Plotting

flist = [20 25 31 40 50 63 80 100 125 160 200 250 315];
y = 63:7:280;

for j = 1:13
figure
plot(10*log10(third_oband(j,:)/(1E-6)^2),y,'linewidth',1.5);
xlabel('Power/Frequency (dB re 1\muPa^2/Hz)');
ylabel('Depth (m)');
ylim([63 280]);
xlim([50 90]);
set (gca,'Ydir','reverse');
set(gca,'Fontsize',30);
title([num2str(flist(j)),' Hz']);
grid on

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf,[pwd ['/SPL_depth/SPL_depth_',num2str(flist(j)),'Hz.fig']]);
end

