
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',22.098,'CutoffFrequency2',27.840,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',27.840,'CutoffFrequency2',35.077,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',35.077,'CutoffFrequency2',44.194,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',44.194,'CutoffFrequency2',55.681,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',55.681,'CutoffFrequency2',70.154,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',70.154,'CutoffFrequency2',88.388,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',88.388,'CutoffFrequency2',111.362,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',111.362,'CutoffFrequency2',140.308,'SampleRate',1000);
%bandfilt =designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',140.308,'CutoffFrequency2',176.777,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',176.777,'CutoffFrequency2',222.725,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',222.725,'CutoffFrequency2',280.616,'SampleRate',1000);
%bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',280.616,'CutoffFrequency2',350,'SampleRate',1000);
bandfilt = designfilt('bandpassfir','FilterOrder',400,'CutoffFrequency1',160,'CutoffFrequency2',320,'SampleRate',1000);

data = array_data{1,1}(33:64,3600000*7:end).';
data = filtfilt(bandfilt,data);
data_fil = data;

%%

FS = 1000;         
NUM_CHANNELS = 32;
view_chn = 32;
array_spacing = 7.5;
c_0 = 1440;
w = 10;
r = 0.01;
pthres = 0.00001;
tbe_list = [];
data_name = 'Tape 23 Hour 1';

f1 = 'event_number';
f2 = 'start_time';
f3 = 'end_time';
f4 = 'event_dur';
f5 = 'num_arrivals';
f6 = 'peak_amp';
event_disp = struct(f1,[],f2,[],f3,[],f4,[],f5,[],f6,[]);

[peak_event_mat,loc_event_mat,num_event] = np_eSelect(data,FS,w,r,pthres,view_chn,array_spacing,c_0,data_name);

%%
% Sort Event location
loc_event = [];
peak_event = [];
for chn = 1:32
    loc_event_temp = loc_event_mat(chn,:)';
    loc_event_temp = loc_event_temp(~isnan(loc_event_temp));
    loc_event = [loc_event;loc_event_temp];
    peak_event_temp = peak_event_mat(chn,:)';
    peak_event_temp = peak_event_temp(~isnan(peak_event_temp));
    peak_event = [peak_event;peak_event_temp];
end

loc_peak_mat = zeros(length(loc_event),2);
loc_peak_mat(:,1) = loc_event;
loc_peak_mat(:,2) = peak_event;
loc_peak_mat = sortrows(loc_peak_mat);

loc_event = loc_peak_mat(:,1);
peak_event = loc_peak_mat(:,2);

time_betw_events = zeros(length(loc_event)-1,1);

count = 0;
if ~isempty(loc_event);
    t_thres = 0.05;
    t_start = loc_event(1);
    ind = 1;
    for t = 2:length(loc_event)
        time_betw_events(t-1) = loc_event(t)-loc_event(t-1);

        if loc_event(t)-loc_event(t-1) > t_thres
            count = count+1;
            event_disp(count).event_number = count;
            event_disp(count).start_time = t_start;
            event_disp(count).end_time = loc_event(t-1);
            event_disp(count).event_dur = loc_event(t-1)-t_start;
            event_disp(count).num_arrivals = ((t-1)-ind+1)/32;
            event_disp(count).peak_amp = max(peak_event(ind:t-1));
            %data(floor(t_start*FS):ceil(loc_event(t-1)*FS),:) = NaN;
            data_fil(floor(t_start*FS)+1:ceil(loc_event(t-1)*FS),:) = NaN;
            t_start = loc_event(t);
            ind = t;
        end     

    end

    count = count+1; 
    event_disp(count).event_number = count;
    event_disp(count).start_time = t_start;
    event_disp(count).end_time = loc_event(end);
    event_disp(count).event_dur = loc_event(end)-t_start;
    event_disp(count).num_arrivals = (length(loc_event)-ind+1)/32;
    event_disp(count).peak_amp = max(peak_event(ind:end));
    %data_fil(floor(t_start*FS):ceil(loc_event(end)*FS),:) = NaN;
    data_fil(floor(t_start*FS)+1:ceil(loc_event(end)*FS),:) = NaN;
end 


%%
data_temp = data_fil(:,32);
data_temp = data_temp(~isnan(data_temp));
FS = 1000; 
window_size = 256;
window = hamming(window_size);
NFFT = 4*window_size;
overlap = 0.5;
df = FS/NFFT;
f = -FS/2:df:FS/2-df;
chn = 1;
data_name = 'Tape 23';

[psd_all,psd_med,freq] = psd(data_temp,window_size,window,overlap,NFFT,FS,chn,data_name);

plot(freq,10*log10(psd_med./1E-12),'linewidth',1.5)
 set(gca,'Fontsize',30);
 title([data_name,' Power Spectral Density Estimate']);
 xlabel('Frequency (Hz)');
 xlim([1 350]);
 ylabel('Power/Frequency (dB re 1\muPa^2/Hz)');
 grid on
 hold on
 
 %%
 time = (1/(12000))*(0:length(data)-1/12000);
 
 figure
 for i = 1:32
     subplot(32,1,i)
     plot(time(113.69*12000:113.73*12000),data(113.69*12000:113.73*12000,i))
 end
