%% NPS data collection
% Thomas Carey, 1/17/2018

% This script begins collecting data at the specified sampling frequency
% when it is run. When the stop button is clicked, it will write all
% acquired data to a .mat file with the specified file name.

% If you close the figure or stop the script (ctrl+C) instead of clicking
% the 'stop' button, data will !not be saved! and the DAQ will continue the
% DAQ will continue to acquire data. To stop the DAQ, type 's.stop' at the
% command line.

% To save to a different folder, first set the current directory to that
% folder. 

% Note: the script will overwrite data if a file with the specified
% filename already exists; make sure to change it. 

%% set the filename and click "run"
% to stop acquisition, press "stop" on the bottom-left corner of the plot
filename = 'device4_200mbar_1mhz.mat';

sampleRate = 1000000; % sampling frequency, Hz

setGlobalFs(sampleRate);
updateFrequency = 1; % number of plot updates per second

% close figure 
f1 = figure(1);
if ishandle(f1)
    close(f1)
end

%% data acquisition
s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');
s.Rate = sampleRate; % rate in Hz
s.NotifyWhenDataAvailableExceeds = sampleRate/updateFrequency;

fid1 = fopen('log.bin','w'); % write acquired data to binary file

btn = uicontrol('Style', 'pushbutton', 'String', 'Stop',...
        'Position', [10 10 40 20],...
        'Callback', 's.stop');%This push button is suppose to stop
                                %the analog input object

figure(1)

% lh1 = addlistener(s,'DataAvailable', @(src,event) plot(event.TimeStamps, event.Data)); % plot data in real time
lh1 = addlistener(s,'DataAvailable', @plotData); % plot data in real time
lh2 = addlistener(s,'DataAvailable', @(src,event) logData(src,event,fid1)); % log data to file

s.IsContinuous = true;
s.startBackground();

while s.IsRunning
    pause(5);
    fprintf('Samples acquired = %d\n', s.ScansAcquired);
end

%after completion, delete listeners reponsible for logging data and
%updating plot
delete(lh1);
delete(lh2);
fclose(fid1);
fprintf('Stopping ... \n'),

%% save data

fid2 = fopen('log.bin','r');
[data,count] = fread(fid2,[2,inf],'double'); % load data into MATLAB
fclose(fid2);


fprintf('Saving, please wait... \n'),
save(filename,'data'); %save .m containing the data + timestamps to file

return

