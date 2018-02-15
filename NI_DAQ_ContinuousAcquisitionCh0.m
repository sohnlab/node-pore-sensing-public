filename = 'testdata.mat';
sampleRate = 100000; % sample rate in Hz
updateFrequency = 10; % update frequency of the plot in Hz

s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');
s.Rate = sampleRate; % rate in Hz

fid1 = fopen('log.bin','w'); % write acquired data to binary file

btn = uicontrol('Style', 'pushbutton', 'String', 'Stop',...
        'Position', [10 10 40 20],...
        'Callback', 's.stop');%This push button is suppose to stop
                                %the analog input object
                                
lh1 = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps, event.Data)); % plot data in real time
lh2 = addlistener(s,'DataAvailable', @(src,event) logData(src,event,fid1)); % log data to file

s.IsContinuous = true;
s.startBackground();

while s.IsRunning
    pause(1)
    fprintf('While loop: Scans acquired = %d\n', s.ScansAcquired)
end

%after completion, delete listeners reponsible for logging data and
%updating plot
delete(lh1);
delete(lh2);
fclose(fid1);

fid2 = fopen('log.bin','r');
[data,count] = fread(fid2,[2,inf],'double'); % load data into MATLAB
fclose(fid2)

save(filename,'data'); %save .m containing the data + timestamps to file