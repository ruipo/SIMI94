
f0 = 'file_name';
f1 = 'band1_slope';
f2 = 'band1_yint';
f3 = 'band1_r2';
f4 = 'band2_slope';
f5 = 'band2_yint';
f6 = 'band2_r2';

psd_fit = struct(f0,[],f1,[],f2,[],f3,[],f4,[],f5,[],f6,[]);

FS = 1000;
window1 = hanning(512);
NFFT1 = 2*length(window1);
window2 = hanning(1024);
window_size = length(window2);
NFFT2 = 2*length(window2);
overlap = 0.5;
NUM_CHANNELS = 16;
chn = 8;
tape_num = 20;

set(0,'DefaultFigureVisible','off')

% Set Path to DATA
prefix = ['/Volumes/icex6/SIMI_Parsed_Data_Poulsen/tape0',num2str(tape_num),'/'];
directory = dir([prefix ['tape0',num2str(tape_num),'_file0*.mat']]);

% bandfilt1 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',57,'CutoffFrequency2',62,'SampleRate',1000);
% bandfilt2 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',116,'CutoffFrequency2',121,'SampleRate',1000);
% bandfilt3 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',66,'CutoffFrequency2',71,'SampleRate',1000);
% bandfilt4 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',76,'CutoffFrequency2',81,'SampleRate',1000);

for i = 1:length(directory)
    
    close all
    disp([num2str(i),' / ', num2str(length(directory))])
    
    filename = [prefix directory(i).name];
    file = importdata(filename);

    timestamp = 725846400 + file.record_headers_concat_file{1,1}.Date(1,2)*24*60*60 + file.record_time_msec_file(1,1)/1000;
    data_name = datestr ((timestamp / 86400) + datenum (1970,1,1), 31);
    
    data = file.array_data_concat_file{1,1}(33:48,:)';
    for j = 1:NUM_CHANNELS
    data(:,j) = data(:,j) - mean(data(:,j));
    end
    clear file
    
    data_fil = data(:,chn);
%     data_fil = filtfilt(bandfilt1,data(:,chn));
%     data_fil = filtfilt(bandfilt2,data_fil);
%     data_fil = filtfilt(bandfilt3,data_fil);
%     data_fil = filtfilt(bandfilt4,data_fil);
    time = (1/(FS))*(0:length(data)-1/FS);
    freq = 0:FS/NFFT2:FS/2;
    
    figure
    spectrogram(160*1E6*data_fil,window1,[],NFFT1,FS,'yaxis');
    title(['Timestamp: ',data_name,'; Channel = ',num2str(chn)]);
    caxis([50 90]);
    set(gca,'Fontsize',20);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,[pwd ['/spectrogram_plot/tape',num2str(tape_num),'_file',num2str(i+2),'.png']]);
    
    
    psdplot = psd(160*data,window_size,window2,overlap,NFFT2,FS,chn,data_name);
    xlim([1 500]);
    ylim([20 120]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %psd_mat(i) = getframe(gcf);
    
    [~,ind0] = min(abs(freq-20));
    [~,ind1] = min(abs(freq-52));
    [~,ind2] = min(abs(freq-62));
    [~,ind3] = min(abs(freq-115));
    [~,ind4] = min(abs(freq-121));
    [~,ind5] = min(abs(freq-174));
    [~,ind6] = min(abs(freq-180));
    [~,ind7] = min(abs(freq-300));
    [~,ind8] = min(abs(freq-65));
    [~,ind9] = min(abs(freq-70));
    [~,ind10] = min(abs(freq-75));
    [~,ind11] = min(abs(freq-80));
    [~,ind12] = min(abs(freq-42));
    [~,ind13] = min(abs(freq-46));
    
    psdplot_fil1 = psdplot([ind0:ind12,ind13:ind1,ind2:ind8,ind9:ind10,ind11:ind3,ind4:ind5,ind6:ind7]);
    freq_fil1 = freq([ind0:ind12,ind13:ind1,ind2:ind8,ind9:ind10,ind11:ind3,ind4:ind5,ind6:ind7]);
    
    psdplot_fil2 = psdplot(ind7:end);
    freq_fil2 = freq(ind7:end);
    
    psd_fit(i).file_name = i+2;
    
    % Fit 20-300Hz
    figure
    plot(freq_fil1,10*log10(psdplot_fil1/1E-12),'o','linewidth',1.5)
    hold on
    set(gca,'Fontsize',20);
    title(['PSD Fit; Time = ',data_name,'; Window size = ', num2str(window_size)]);
    xlabel('Frequency (Hz)');
    xlim([20 300])
    ylabel('Power/Freq (dB/Hz)');
    ylim([40 100]);
    
    [p,s] = polyfit(freq_fil1,10*log10(psdplot_fil1/1E-12),1);
    r2 = 1 - s.normr^2 / norm(10*log10(psdplot_fil1/1E-12)-mean(10*log10(psdplot_fil1/1E-12)))^2;
    plot(freq_fil1,p(1)*freq_fil1+p(2),'r-','linewidth',1.5)
    equation = ['y = ',num2str(p(1)),'x + ',num2str(p(2)),';   r^2 = ',num2str(r2)];
    text(100,90,equation,'FontSize',20);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,[pwd ['/band1_psd_fit/tape',num2str(tape_num),'_file',num2str(i+2),'.png']]);
    
    psd_fit(i).band1_slope = p(1);
    psd_fit(i).band1_yint = p(2);
    psd_fit(i).band1_r2 = r2;
    
    
    
    % Fit 300-500Hz
    figure
    plot(freq_fil2,10*log10(psdplot_fil2/1E-12),'o','linewidth',1.5)
    hold on
    set(gca,'Fontsize',20);
    title(['PSD Fit; Time = ',data_name,'; Window size = ', num2str(window_size)]);
    xlabel('Frequency (Hz)');
    xlim([300 500])
    ylabel('Power/Freq (dB/Hz)');
    ylim([20 90]);
    
    [p,s] = polyfit(freq_fil2,10*log10(psdplot_fil2/1E-12),1);
    r2 = 1 - s.normr^2 / norm(10*log10(psdplot_fil2/1E-12)-mean(10*log10(psdplot_fil2/1E-12)))^2;
    plot(freq_fil2,p(1)*freq_fil2+p(2),'r-')
    equation = ['y = ',num2str(p(1)),'x + ',num2str(p(2)),';   r^2 = ',num2str(r2)];
    text(420,60,equation,'FontSize',20);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    saveas(gcf,[pwd ['/band2_psd_fit/tape',num2str(tape_num),'_file',num2str(i+2),'.png']]);
    
    psd_fit(i).band2_slope = p(1);
    psd_fit(i).band2_yint = p(2);
    psd_fit(i).band2_r2 = r2;
    

end

% v_psd = VideoWriter(['SIMI94_PSD_tape',num2str(tape_num),'_ti819s.avi','Uncompressed AVI']);
% v_psd.FrameRate = 3;
% open(v_psd)
% writeVideo(v_psd,psd_mat)
% close(v_psd)
% 
