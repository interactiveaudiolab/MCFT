%% Example(1):
 % Generate and plot upward and downward STRFs with a scale of
 % 1 cycle per octave and a rate of 4 Hz
 
S=1; R=4;
 
S_params=struct('hslen',128,'ripple_freq',24,'type','bandpass');
R_params=struct('time_const',3.5,'hrlen',500,'frame_per_sec',250,'type','bandpass');

[h_up,~]=gen_hsr(S,R,S_params,R_params,'up');
[h_down,~]=gen_hsr(S,R,S_params,R_params,'down');

h_up=fftshift(h_up,1);
h_down=fftshift(h_down,1);

ytl={'\bf.25 f0','\bf.5 f0','\bf1 f0','\bf2 f0','\bf4 f0','\bf8 f0'};
xtl={'\bf0.2','\bf0.4','\bf0.6','\bf0.8','\bf1.0'};

figure;
subplot(211)
imagesc(abs(h_up))
set(gca,'ydir','normal','yticklabel',ytl,'xticklabel',xtl)
colormap('jet')
ylabel('\bf Log Frequency (Hz)')
xlabel('Time (sec)')
title('Upward Ripple')


subplot(212)
colormap('jet')
imagesc(abs(h_down))
set(gca,'ydir','normal','yticklabel',ytl,'xticklabel',xtl)
ylabel('scale (cyc/oct)')
ylabel('scale (cyc/oct)')
ylabel('\bf{Log Frequency}')
xlabel('\bf{Time (sec)}')
title('Downward Ripple')
