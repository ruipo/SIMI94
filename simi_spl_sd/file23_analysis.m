% Element Locations
N = 32;
% 0.75m spacing in the middle
z = zeros(N,1);
d = 7;
for n = 0:N-1
    z(n+1) = -(n-(N-1)/2)*d; 
end

p = [zeros(1,N) ; zeros(1,N) ; fliplr(z')];
p = p';

% Load Data
data = array_data{1,1}(33:64,:).';
data = (data.*(1/(10^(-175/20))));

%%
starttime_str = '1994-04-18 21:39:14';
starttime = datetime(starttime_str);
starttime_utc = posixtime(starttime);

window_size = 512;
window = hamming(window_size);
NFFT = 4*window_size;
overlap = 0.5;
elev = -90:1:90;
az = 0;
c = 1445;
windowb = hanning(2048);
NFFTb = 128;
f_range = [20 100];
weighting = 'icex_hanning';
FS = 1000;
samp_len = 30*FS;
data_len = size(data,1);
channels = 1:1:32;

fr = linspace(20,100,NFFTb);
f1 = 80;
f2 = 90;
[~,ind1] = min(abs(fr-f1));
[~,ind2] = min(abs(fr-f2));

aco_in = zeros(samp_len, 32);
time_vec = zeros(1,round(data_len/samp_len)-1);
psd_mat = zeros(round(data_len/samp_len)-1,overlap*NFFT+1,size(data,2));
beamform_mat = zeros(round(data_len/samp_len)-1,length(elev));
beamform_f_mat = zeros(length(elev),NFFTb,round(data_len/samp_len)-1);

for index = 0:round(data_len/samp_len)-2
    
    disp([num2str(index),' / ', num2str(round(data_len/samp_len)-1)])

    % Nomalized to zero mean;
    aco_in = data(index*samp_len+1:(index+1)*samp_len,:);
    for j = 1:32       
        data_in(:,j) = ((aco_in(:,j)-mean(aco_in(:,j))))./10^6;
    end
    
    time_vec(index+1) = (samp_len/FS/2) + index*(samp_len/FS);
    
    % PSD
    [psd_chn,~,~,freq] = psd(data_in,window_size,window,overlap,NFFT,FS,channels,'test');
    
    psd_mat(index+1,:,:) = psd_chn;
    
    % BEAMFORMING

    [beamform_output,t,t_end] = beamform_3D(data_in,p,FS,elev,az,c,f_range,NFFTb,windowb,overlap,weighting);
    beamform_f = squeeze(mean(beamform_output,1));
    beamform = squeeze(mean(mean(beamform_output(:,:,:,ind1:ind2,1)),4));
    
    beamform_mat(index+1,:) = beamform;
    beamform_f_mat(:,:,index+1) = beamform_f;

end

%%

timevec_utc = time_vec+starttime_utc;
datetimevec = datetime(timevec_utc, 'convertfrom','posixtime');

beamform_mat(end,:) = [];
beamform_f_mat(:,:,end) = [];

%% Plotting - time vs freq @ bf angle

figpath = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/SIMI94/simi_spl_sd/ang_figures/';

fr = linspace(20,100,NFFTb);
f1 = 80;
f2 = 90;
[~,ind1] = min(abs(fr-f1));
[~,ind2] = min(abs(fr-f2));

bf_all = squeeze(mean(median(beamform_f_mat(:,ind1:ind2,:),3),2));

figure('units','normalized','outerposition',[0 0 1 1]);
s1 = subplot(3,3,[1,2,4,5]);
ylabel('Frequency (Hz)');
xlabel('Time (04-18/19-1994 UTC)')
datetick('x', 'HH:MM','keeplimits');
ylim([20 100])
xlim([datenum(datetimevec(1)) datenum(datetimevec(end))])
%Xticks([])
%set(gca,'Xticklabel',[]);
%set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
caxis([80,125])
colorbar
a = colorbar;
a.Label.String = 'dB';
a.Location = 'northoutside';
a.Direction = 'reverse';
colormap jet
hold on

s2 = subplot(3,3,[7,8]);
ylim([80 125])
xlim([20 100])
xlabel('Frequency (Hz)')
ylabel('NL (dB re 1\muPa^2/Hz)')
set(gca,'fontsize',15)
grid on
hold on

s3 = subplot(3,3,[3,6,9]);
plot(10*log10(bf_all/(1E-6)^2),elev,'b','linewidth',2)
ylim([-90 90])
xlim([80 125])
yticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Yticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]) 
grid on
ylabel('Elevation (Degrees)');
xlabel('NL (dB re 1\muPa^2/Hz)')
set(gca,'fontsize',15)
title('Frequency = 850 Hz')
hold on


for ee = 1:length(elev)
    ee
    figure(1)
    axes(s1);
    h = pcolor(datenum(datetimevec),fr,10*log10(squeeze(beamform_f_mat(ee,:,:))/(1E-6)^2));
    set(h,'Edgecolor','None')
    datetick('x', 'HH:MM','keeplimits');
    ylim([20 100])
    xlim([datenum(datetimevec(1)) datenum(datetimevec(end))])
    title(['Elevation = ' num2str(elev(ee)) ' Degrees'])

    axes(s2)
    f1 = plot(fr,10*log10(squeeze(mean(beamform_f_mat(ee,:,:),3))/(1E-6)^2),'k','linewidth',2);
    
    axes(s3);
    f2 = plot([70 130],ones(2,1)*elev(ee), 'k','linewidth',1.5);
    
    saveas(gcf,[figpath num2str(ee) '.png']);
    pause(0.1)
    
    delete(h);
    delete(f1);
    delete(f2);
    
end

%% make video - ANG

testf = dir([figpath '*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([figpath 'ang.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1 : length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);

%% plotting - Depth vs dB @ frequency
figure
depths = -(-p(:,3)-max(p(:,3))-30);

flist = [20 40 60 80 100];
color = ['k','r','b','m','g'];
%shape = ['o','*','+','x'];
for ff = 1:length(flist)
    f = flist(ff);
    [~,ind] = min(abs(f-freq));

    plot(10*log10(squeeze(mean(psd_mat(:,ind,:),1))/(1E-6)^2),depths,[color(ff)],'linewidth',1)
    hold on
    title('Noise Level vs. Depth');
    xlabel('NL (dB re 1\muPa^2/Hz)')
    xlim([100 115])
    ylabel('Depth (m)')
    yticks([0 25 50 75 100 125 150 175 200 225 250 275])
    set(gca,'Yticklabel',[0 25 50 75 100 125 150 175 200 225 250 275]);
    set(gca,'Ydir','reverse')
    grid on
    set(gca,'fontsize',20);
end

legend('20Hz', '40Hz','60Hz','80Hz','100Hz');
%% Plotting - NL

figpath = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/SIMI94/simi_spl_sd/figures/';

figure('units','normalized','outerposition',[0 0 1 1]);

s2 = subplot(3,2,[1,2,3,4]);
h = pcolor(freq,datenum(datetimevec),10*log10(psd_mat(:,:,16)/(1E-6)^2));
xlabel('Frequency (Hz)');
ylabel('Time (04-18/19-1994 UTC)')
datetick('y', 'HH:MM','keeplimits');
xlim([20 100])
set(h,'Edgecolor','None')
%set(gca,'XScale','log')
set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
%xticks([50 100 200 400 800 1600 3200 6000])
%set(gca,'Xticklabel',[50 100 200 400 800 1600 3200 6000]) 
caxis([80,125])
colorbar
a = colorbar;
a.Label.String = 'dB';
a.Location = 'northoutside';
a.Direction = 'reverse';
colormap jet
title('Depth = 135m')
hold on

s3 = subplot(3,2,[5,6]);
xlim([20 100])
ylim([80 125])
%xticks([50 100 200 400 800 1600 3200 6000])
%set(gca,'Xticklabel',[50 100 200 400 800 1600 3200 6000]) 
grid on
%set(gca,'XScale','log')
xlabel('Frequency (Hz)');
ylabel('NL (dB re 1\muPa^2/Hz)')
set(gca,'fontsize',15)
hold on


for tt = 1:length(datetimevec)
    tt
    figure(1)

    axes(s2);
    f1 = plot([20 100],ones(2,1)*datenum(datetimevec(tt)),'k','linewidth',2);

    axes(s3);
    f2 = plot(freq,10*log10(psd_mat(tt,:,16)/(1E-6)^2),'b','linewidth',2);
    
    saveas(gcf,[figpath num2str(tt) '.png']);
    pause(0.1)

    delete(f1);
    delete(f2);
end


%% make video - NL

testf = dir([figpath '*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([figpath 'nl.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1 : length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);
  
%% Plotting - depth vs bf max
figure
hold on
title('Max NL Elevation Angle vs. Depth (Frequency = 85Hz)');
xlabel('Elevation (Degrees)')
ylabel('Time (UTC)')
xlim([-90 90])
xticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Xticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]);
set(gca,'Ydir','reverse')
grid on
set(gca,'fontsize',20);

for ii = 1:size(beamform_mat,1)
    [~, inds] = max(10*log10(abs(beamform_mat(ii,:))/1E-12));
    
    plot(elev(inds),datetimevec(ii),'bo','linewidth',1.5);
end
    
%% Plotting - depth vs bf

figpath = '/Users/Rui/Documents/Graduate/Research/ICEX:SIMI/SIMI94/simi_spl_sd/bf_figures/';

figure('units','normalized','outerposition',[0 0 1 1]);

s2 = subplot(3,2,[1,2,3,4]);
h = pcolor(elev,datenum(datetimevec),10*log10(abs(beamform_mat)/1E-12));
xlabel('Elevation (Degrees)');
ylabel('Time (04-18/19-1994 UTC)')
datetick('y', 'HH:MM','keeplimits');
xlim([-90 90])
xticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Xticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]);
set(h,'Edgecolor','None')
set(gca,'Ydir','reverse')
set(gca,'fontsize',15)
caxis([80,125])
colorbar
a = colorbar;
a.Label.String = 'dB';
a.Location = 'northoutside';
a.Direction = 'reverse';
colormap jet
hold on
title('Frequency = 85Hz')

s3 = subplot(3,2,[5,6]);
xlim([-90 90])
xticks([-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90])
set(gca,'Xticklabel',[-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90]);
ylim([80,125])
grid on
xlabel('Elevation (Degrees)');
ylabel('NL (dB re 1\muPa^2/Hz)')
set(gca,'fontsize',15)
hold on


for tt = 1:length(datetimevec)
    tt
    figure(1)

    axes(s2);
    f1 = plot([-90 90],ones(2,1)*datenum(datetimevec(tt)),'k','linewidth',2);

    axes(s3);
    f2 = plot(elev,10*log10(abs(beamform_mat(tt,:))/1E-12),'b','linewidth',2);
    
    saveas(gcf,[figpath num2str(tt) '.png']);
    pause(0.1)

    delete(f1);
    delete(f2);
end

%% make video - BF

testf = dir([figpath '*.png']);
[~, reindex] = sort( str2double( regexp( {testf.name}, '\d+', 'match', 'once' )));
testf = testf(reindex);

% Create a VideoWriter object to write the video out to a new, different file.
  writerObj = VideoWriter([figpath 'bf.avi']);
  writerObj.FrameRate = 8;
  open(writerObj);
  
  for frame = 1 : length(testf)
    disp([num2str(frame) '/' num2str(length(testf))])
    % Construct an output image file name.
    outputBaseFileName = sprintf(testf(frame).name);
    outputFullFileName = fullfile(testf(frame).folder, outputBaseFileName);
    % Read the image in from disk.
    thisFrame = imread(outputFullFileName);
    % Convert the image into a "movie frame" structure.
    recalledMovie(frame) = im2frame(thisFrame);
    % Write this frame out to a new video file.
    writeVideo(writerObj, thisFrame);
  end
  close(writerObj);



