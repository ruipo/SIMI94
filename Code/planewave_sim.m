
FS = 1000;
c = 1435;
win_len = 512;
overlap = 0.5;
daz = 2*pi/500;
azlist = 0:daz:2*pi-daz;

%sUCA = phased.UCA('NumElements',12,'Radius',30.0);
sUCA = phased.ULA('NumElements',32,'ElementSpacing',7,'ArrayAxis','z');
t = (0:0.001:100)';

p = getElementPosition(sUCA);
p(3,:) = [];

delay = phased.ElementDelay('SensorArray',sUCA,'PropagationSpeed',c);
tau =  step(delay,[0;30]);
tau = tau-tau(1);

amp1 = rand(100001,1);
amp2 = rand(100001,1);
amp3 = rand(100001,1);
for i = 1:length(tau)
data(:,i) = amp1.*sin(2*pi*50*(t-tau(i))) + amp2.*sin(2*pi*85.450*(t-tau(i))) + amp3.*sin(2*pi*199.30*(t-tau(i)));
end


[time_beamform,t] = time_beamform_2D(data,p,FS,azlist,c,win_len,overlap);

%%
figure
t_average = mean(time_beamform,2);
polarplot([azlist 2*pi],[t_average; t_average(1)]);

%%
figure
for i = 1:4
t_average = time_beamform(:,i);
polarplot([azlist 2*pi],[t_average; t_average(1)]);
set(gca,'Fontsize',30);
pause
end