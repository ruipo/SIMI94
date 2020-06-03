FS = 1000;
window = hanning(1024);
NFFT = 2*length(window);
elev = -90:1:90;
z = [-52.5 -45.5 -38.5 -31.5 -24.5 -17.5 -10.5 -3.5 3.5 10.5 17.5 24.5 31.5 38.5 45.5 52.5];
overlap = 0.5;
taper = ones(length(z),1);
f_range = [80 100];
c_0 = 1435;
adaptive = 0;
btr = 1;
lofar = 0;
fraz = 0;

NUM_CHANNELS = 16;
chn = 8;

%set(0,'DefaultFigureVisible','off')

% Set Path to DATA
prefix = '/Volumes/icex6/SIMI_Parsed_Data_Poulsen/tape015/';
directory = dir([prefix 'tape015_file0*.mat']);

bandfilt1 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',57,'CutoffFrequency2',62,'SampleRate',1000);
bandfilt2 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',116,'CutoffFrequency2',121,'SampleRate',1000);
bandfilt3 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',66,'CutoffFrequency2',71,'SampleRate',1000);
bandfilt4 = designfilt('bandstopfir','FilterOrder',500,'CutoffFrequency1',76,'CutoffFrequency2',81,'SampleRate',1000);

for i = 1;
    
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
    
    data_fil = filtfilt(bandfilt1,data);
    data_fil = filtfilt(bandfilt2,data_fil);
    data_fil = filtfilt(bandfilt3,data_fil);
    data_fil = filtfilt(bandfilt4,data_fil);
    time = (1/(FS))*(0:length(data)-1/FS);
    
    [beamform_elev,beamform_elev_allf,f,t] = vert_array_beamform(160*data_fil,elev,z,window,overlap,NFFT,FS,taper,f_range,c_0,adaptive,btr,lofar,fraz,data_name);
end
