function plotDataSpectrum(src,event)
    sampleRate = getGlobalFs;
    nfft = getGlobalNFFT;
    psd = getGlobalPSD;
    window = 1; % plot window
    
    %% resample data
%     [lpf_b, lpf_a] = butter(6,500/(sampleRate/2),'low');
%     lpData = filtfilt(lpf_b,lpf_a,x);
%     [bs_b, bs_a] = butter(2,[58 62]/(sampleRate/2),'stop');
%     f_Data = filtfilt(bs_b,bs_a,lpData);
%     M = 100;
%     f_Data = resample(x, 1, M);
%     f_Data = x; % don't resample
%     M = 10; % resampling ratio
%     f_Data = resample(f_Data, 1, M);
%     sampleRate = sampleRate / M;

    %%% Alan added this 2018-02-09
    t = event.TimeStamps;
    x = event.Data;
%     n = length(x);
%     T = n/sampleRate;
    resampleRate = 5000; % samp/s
    M = sampleRate/resampleRate; % resampling factor
    P = 500; % padding factor?
    x = [x(end-(P*M-1):end); x; x(1:(P*M))]; % pad
    x_r = resample(x, 1, M); % resample
%     t_r = t(1:M:end); % resampled timestamps
%     x_r_2 = interp1(t_r, x_r, linspace(t_r(1), t_r(end), T*resampleRate)); % fix integer nonsense
%     f0 = 60; Q = 35; BW = (f0/(resampleRate/2))/Q; % make 60 Hz comb filter
%     [b, a] = iircomb(ceil(resampleRate/f0), BW, 'notch'); % Note type flag 'notch'
    f0 = 60; Q = 1; BW = (f0/(resampleRate/2))/Q; % make 60 Hz notch filter
    [b, a] = iirnotch((f0/resampleRate/2), BW);
    x_r_n = filter(b, a, x_r); % notch filter
    x_r_n = x_r_n(1+P:end-P); % remove zero-pad
    f_Data = x_r_n;
%     pause(0);
%%%

    %% plot raw data
    set(gcf, 'units', 'normalized', 'outerposition', [0.1 0.1 0.9 0.9]); % fullscreen
    subplot(3,1,1)
    h1 = animatedline('MaximumNumPoints',100000);
    axis_x2 = max(t);
    axis_x1 = axis_x2-window;
    axis_y1 = min(event.Data);
    axis_y2 = max(event.Data);

    addpoints(h1,t,event.Data);
    [~,ydata] = getpoints(h1);
    axis_y1 = min(ydata);
    axis_y2 = max(ydata);
    if axis_x1 >= 0
        axis([axis_x1,axis_x2,axis_y1,axis_y2]);
    else
        axis([0,window,axis_y1,axis_y2]);
    end
%     set(gca,'xlim',[t(1), t(end)]);
    ax = axis;
%     drawnow
    title('Raw data');
    xlabel('Time (s)','FontSize',12);
    ylabel('V','FontSize',12);
    
     %% plot filtered data
    subplot(3,1,2)
    h2 = animatedline('MaximumNumPoints',ceil(1000000/M));
    addpoints(h2, linspace(t(1), t(end), length(f_Data)), f_Data);
%     axis_y1 = min(ydata);
%     axis_y2 = max(ydata);
%     if axis_x1 >= 0
%         axis([axis_x1,axis_x2,axis_y1,axis_y2]);
%     else
%         axis([0,window,axis_y1,axis_y2]);
%     end
    axis(ax);
    drawnow
    title(sprintf('Resampled to f_{max} = %d Hz and 60 Hz notch filtered', resampleRate/2));
    xlabel('Time (s)','FontSize',12);
    ylabel('V','FontSize',12);
    
    %% plot spectrum
    subplot(3,1,3)
%     [psd_new, freq] = pwelch(f_Data, [], [], nfft, sampleRate/M);
    sig = event.Data;
    freq = (0:nfft/2) * sampleRate/nfft;
    psd_new = fft( hanning(length(sig)) .* (sig - mean(sig)), nfft);
    psd_new = abs(psd_new(1:nfft/2+1)).^2 / (sampleRate/nfft);
    lambda = 0.5;
    psd = (1-lambda)*psd_new + lambda*psd;
    psd_dB = 10*log10(psd);
    setGlobalPSD(psd);
%     fprintf('%f\n', psd_dB(1));
%     h2 = animatedline('MaximumNumPoints', length(freq));
%     clearpoints(h2);
%     addpoints(h2, freq, psd_dB);
    plot(freq, psd_dB);
    set(gca, 'xscale', 'log');
%     [~,ydata] = getpoints(h1);
%     axis_y1 = min(psd_dB);
%     axis_y2 = max(psd_dB);
%     axis_x1 = freq(1);
%     axis_x2 = freq(end);
%     axis([axis_x1,axis_x2,axis_y1,axis_y2]);
%     set(gca,'xscale','log');
    grid on;
    grid minor;
    axis tight;
    drawnow
    title('Power Spectral Density')
    xlabel('Frequency (Hz)','FontSize',12);
    ylabel('PSD (dB/Hz)','FontSize',12);
end