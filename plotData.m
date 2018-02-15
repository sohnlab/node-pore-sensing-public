function plotData(src,event)
    sampleRate = getGlobalFs;
    window = 1; % plot window
    padding = 0; % proportion of signal used for y-axis padding
%     downsampleFactor = 1;
%     tempData = downsample(event.Data,downsampleFactor); % downsample by 10
    
    %% filter data
%     tempDataExp = [tempData(end-1000:end);tempData;tempData(1:1000)]; % add elements to beginning and end of tempData for filtering to eliminate ringing edge effects
%     [lpf_b, lpf_a] = butter(6,500/(sampleRate/(2*downsampleFactor)),'low');
%     lpData = filtfilt(lpf_b,lpf_a,tempDataExp);
%     [bs_b, bs_a] = butter(2,[55 65]/(sampleRate/(2*downsampleFactor)),'stop');
%     f_Data = filtfilt(bs_b,bs_a,lpData);
%     f_Data = f_Data(1001:end-1001); % remove elements we artificially added
%     f_Data = upsample(f_Data,downsampleFactor); % upsample data by factor for plotting
    
%%% Alan added this 2018-02-08
    t = event.TimeStamps;
    x = event.Data;
%     n = length(x);
%     T = n/sampleRate;
    resampleRate = 1000; % samp/s
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
    
%     [lpf_b, lpf_a] = butter(6,500/(sampleRate/2),'low');
%     lpData = filtfilt(lpf_b,lpf_a,event.Data);
%     [bs_b, bs_a] = butter(2,[55 65]/(sampleRate/2),'stop');
%     f_Data = filtfilt(bs_b,bs_a,lpData);

    %% plot raw data
    subplot(2,1,1)
    h1 = animatedline('MaximumNumPoints',1000000);
    minax = min(event.Data)*(1-padding);
    maxax = max(event.Data)*(1+padding);
    axis_x2 = max(event.TimeStamps);
    axis_x1 = axis_x2-window;
    axis_y1 = minax;
    axis_y2 = maxax;

    addpoints(h1,event.TimeStamps,event.Data);
    if axis_x1 >= 0
        axis([axis_x1,axis_x2,axis_y1,axis_y2]);
    else
        axis([0,window,axis_y1,axis_y2]);
    end
    ax = axis;
%     drawnow
    title('Raw data');
    xlabel('time (s)','FontSize',12);
    ylabel('V','FontSize',12);
    
    %% plot filtered data
    subplot(2,1,2)
    h2 = animatedline('MaximumNumPoints',ceil(1000000/M));
    addpoints(h2, linspace(event.TimeStamps(1), event.TimeStamps(end), length(f_Data)), f_Data);
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
    xlabel('time(s)','FontSize',12);
    ylabel('V','FontSize',12);
end