

band1filt = designfilt('highpassfir', 'FilterOrder', 500, 'StopbandFrequency', 20, 'PassbandFrequency', 30, 'SampleRate', 1000);
band2filt = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',42,'CutoffFrequency2',46,'SampleRate',1000);
band3filt = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',52,'CutoffFrequency2',62,'SampleRate',1000);
band4filt = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',65,'CutoffFrequency2',70,'SampleRate',1000);
band5filt = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',75,'CutoffFrequency2',80,'SampleRate',1000);
band6filt = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',115,'CutoffFrequency2',121,'SampleRate',1000);
band7filt = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',174,'CutoffFrequency2',180,'SampleRate',1000);

FS = 1000;
NUM_CHANNELS = 16;
chn = 8;
tape_num = 18;

set(0,'DefaultFigureVisible','on')

% Set Path to DATA
prefix = ['/Volumes/icex6/SIMI_Parsed_Data_Poulsen/tape0',num2str(tape_num),'/'];
directory = dir([prefix ['tape0',num2str(tape_num),'_file0*.mat']]);

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
    
    data_fil = filtfilt(band1filt,data(:,chn));
    data_fil = filtfilt(band2filt,data_fil);
    data_fil = filtfilt(band3filt,data_fil);
    data_fil = filtfilt(band4filt,data_fil);
    data_fil = filtfilt(band5filt,data_fil);
    data_fil = filtfilt(band6filt,data_fil);
    data_fil = filtfilt(band7filt,data_fil);
    
    time = (1/(FS))*(0:length(data_fil)-1/FS);
    
    figure
    plot(time,160*data_fil);
    set(gca,'Fontsize',20);
    xlabel('Time (s)');
    ylabel('Amplitude (Pa)');
    title(['Filtered TS (20-500Hz); Timestamp = ',data_name]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    %saveas(gcf,[pwd ['/ts_plot/tape',num2str(tape_num),'_file',num2str(i+2),'.png']]);
    
end