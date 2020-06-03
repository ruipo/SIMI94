
%%
FS = 1000; 
window_size = 512;
window = hanning(window_size);
NFFT = window_size;
overlap = 0.5;
df = FS/NFFT;
f = -FS/2:df:FS/2-df;
chn = 1;
data_name = 'Tape 23 Hydrophone 48';

data = array_data{1,1}(40,1:1*3600000).';
data = (data.*(1/((10^(-175/20))/1E-6)));

dt = 1/FS;
t = 0:dt:length(data)/FS-dt;

[psd_all,psd_med,freq] = psd(data,window_size,window,overlap,NFFT,FS,chn,data_name);
%psd_75 = quantile(psd_all,0.75);
%psd_25 = quantile(psd_all,0.25);
%psd_90 = quantile(psd_all,0.90);
%psd_10 = quantile(psd_all,0.10);

figure
plot(freq_simi,10*log10(psd_med_SIMI./1E-12),'b','linewidth',1.5)
 set(gca,'Fontsize',30);
 title([data_name,' Power Spectral Density Estimate']);
 xlabel('Frequency (Hz)');
 xlim([1 350]);
 ylabel('Power/Frequency (dB re 1\muPa^2/Hz)');
 grid on
 %hold on
%plot(freq,10*log10(psd_75./1E-12),'r--','linewidth',1.5);
%plot(freq,10*log10(psd_25./1E-12),'r--','linewidth',1.5);
%plot(freq,10*log10(psd_90./1E-12),'k-.','linewidth',1.5);
%plot(freq,10*log10(psd_10./1E-12),'k-.','linewidth',1.5); 
%legend('Median','75th percentile','25th percentile','90th percentile','10th percentile')

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

third_oband(1,:) = squeeze(mean(psd_all(:,ind1:ind2),2)); %20Hz
third_oband(2,:) = squeeze(mean(psd_all(:,ind2:ind3),2)); %25Hz  
third_oband(3,:) = squeeze(mean(psd_all(:,ind3:ind4),2)); %31Hz
third_oband(4,:) = squeeze(mean(psd_all(:,ind4:ind5),2)); %40Hz
third_oband(5,:) = squeeze(mean(psd_all(:,ind5:ind6),2)); %50Hz
third_oband(6,:) = squeeze(mean(psd_all(:,ind6:ind7),2)); %63Hz
third_oband(7,:) = squeeze(mean(psd_all(:,ind7:ind8),2)); %80Hz
third_oband(8,:) = squeeze(mean(psd_all(:,ind8:ind9),2)); %100Hz
third_oband(9,:) = squeeze(mean(psd_all(:,ind9:ind10),2)); %125Hz
third_oband(10,:) = squeeze(mean(psd_all(:,ind10:ind11),2)); %160Hz
third_oband(11,:) = squeeze(mean(psd_all(:,ind11:ind12),2)); %200Hz
third_oband(12,:) = squeeze(mean(psd_all(:,ind12:ind13),2)); %250Hz
third_oband(13,:) = squeeze(mean(psd_all(:,ind13:ind14),2)); %315Hz


%% Plotting

flist = [20 25 31 40 50 63 80 100 125 160 200 250 315];

for j = 1:13
figure
h1 = histfit(10*log10(third_oband(j,:)/(1E-6)^2),[],'kernel');
x = h1(2).XData;
y = (h1(2).YData/norm(h1(2).YData,1));
plot(x,y,'linewidth',1.5);
delete(h1(1));
delete(h1(2));
xlabel('Power/Frequency (dB re 1\muPa^2/Hz)');
ylabel('Probability Density');
set(gca,'Fontsize',30);
title([num2str(flist(j)),' Hz']);
grid on

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
saveas(gcf,[pwd ['/Tape_31/PSD_PDF_',num2str(flist(j)),'Hz.fig']]);
end

%%
flist = [20 25 31 40 50 63 80 100 125 160 200 250 315];
for i = 1:13
kurtlist(i) = kurtosis(10*log10(third_oband(i,:)/(1E-6)^2));
end

figure
plot(flist,kurtlist,'linewidth',2)
xlabel('Frequency (Hz)');
ylabel('Kurtosis');
set(gca,'Fontsize',30);
title(data_name);
grid on

%%
w = 512;
for k = 1
    figure
     %Spectrogram usage: spectrogram(signal, window, number of overlaps,
     %number of frequency points, sampling frequency)
    spectrogram(1E6*data(:,k),hanning(w),w/2,w,1000,'yaxis');
    %set(gcf,'Position',[1          25        1600         962]);
    title(['Element ' num2str(k)])
    hcolor=colorbar;
    set(hcolor,'FontSize',30);
    caxis([70 120]);
    ylim([0 300]);
   
end