
data = array_data{1,1}(33:64,1:6*3600000);
data = (data.*(1/((10^(-175/20))/1E-6)));
data = data.';
window = hanning(512);
overlap = 0.5;
NFFT = length(window);
FS = 1000;
c = 1435;
N = 32;
d = 7;
z = zeros(N,1);

for n = 0:N-1
    z(n+1) = -(n-(N-1)/2)*d; 
end
z = flipud(z);

%%
[P,freq,k,cov_mat] = f_k_spec_hr(data,window,overlap,NFFT,FS,c,z);

%%
figure
fig = pcolor(k,freq,10*log10(abs(P)/1E-12));
set(fig,'Edgecolor', 'none');
set(gca,'Fontsize',30);
xlabel('Wavenumber')
ylabel('Frequency (Hz)')
colorbar;
colormap 'jet';

hold on
f1list = fliplr(freq);
f1list(end) = [];
flist = [f1list,freq];
plot(k,flist,'k','linewidth',3)

% alpha = 1/2200;
% k1 = freq*alpha;
% k2 = -fliplr(k1);
% k2(end) = [];
% k_1 = [k2 k1];
% plot(k_1,flist,'k','linewidth',3)


%%
finterest = 1:2:80;

for i = 1:length(finterest)
    
flist = freq - finterest(i);
[~,loc] = min(abs(flist));
kinterest = freq(loc)/c;
dk = k(2)-k(1);

klist = k-(-0.07);
[~,loc1] = min(abs(klist));

klist = k-(-kinterest);
[~,loc2] = min(abs(klist));

klist = k-kinterest;
[~,loc3] = min(abs(klist));

klist = k-(0.07);
[~,loc4] = min(abs(klist));


P_outside(i) = 10*log10(mean(mean(abs(P(loc,loc1:loc2))) + mean(abs(P(loc,loc3:loc4))))/1E-12);
P_cone(i) = 10*log10(mean(abs(P(loc,loc2:loc3)))/1E-12);
end

figure
hold on
plot(finterest,P_cone)
hold on
plot(finterest,P_outside)
hold on
plot(finterest,P_cone-P_outside)