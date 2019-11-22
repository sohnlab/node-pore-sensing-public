function [ output_matrix ] = mNPS_procKim( filepath, thresholds )
% [ output_matrix ] = sNPS( start, number_of_files, thresholds )
%   Reads all sNPS data, analyzes data and returns final output matrix
%   Needs a vector of 2 thresholds for initial thresholding. Afterwards,
%   QC, thresholding, and analysis for individual pulses will need to
%   proceed with user input
    
    %% apply default thresholds
    
    load(filepath,'data');
    data = data(2,:);
    sampleRate = 60000;
    
    if nargin < 2
        thresholds = [1e-4, 1e-3];
        fprintf('Auto thresholds set to %3.2e, %3.2e\n',thresholds);
    end
    
    %% read all, search for pulses
    all_out = mNPS_readKim(data, sampleRate, thresholds, false, false);
    
    %% remove duplicate files to read
        
    [uni_win, output_matrix] = mNPS_cleanKim(all_out);
    
    clear all_out,
    %% analyze all time windows, prompt user for input to see if the pulse looks good
    
    good_index = 1; % index through output of good pulses
    i = 1;
    new_th = thresholds; % reset threshold values
    warning('off', 'curvefit:cfit:subsasgn:coeffsClearingConfBounds'); % annoying warning when fit fails to converge
    
    while (i < size(uni_win,1))
        
        % checking for working pulses
        searchflag = true; % true to search, false to stop
        retryflag = true; % true for first attempt, false if trying again with default thresholds
        
        while (searchflag)
            try
                fprintf('\nNow reading index: %d\n',uni_win(i));
                fprintf('Progress: %2.1f %%\n',i/length(uni_win)*100);
                
                iterdata = data(20*(uni_win(i)-200):20*(uni_win(i)+800));
                % measure a new pulse
                [iter_out, emptyflag] = mNPS_readKim(iterdata, sampleRate, new_th, false, false); 
                % iter_out: output of one iteration
                % emptyflag: skip pulse if TRUE
                % auto: values for computing auto-threshold value
                % 
                % function arguments: read iterdata from uni_win entries
                % plot and fit off for speed
                
                
                if emptyflag
                    fprintf('No pulse! Skipping...\n');
                    i = i+1;
                    searchflag = true;
                else
                    [iter_out, ~, auto] = mNPS_readKim(iterdata, sampleRate, new_th, true, true);
                    searchflag = false;
                end
                
            catch ME
                
                %% print errors to command line for debug purposes
                fprintf('-----\n%s\n',ME.identifier),
                for errorstack_i = 1:length(ME.stack)

                    fprintf('Line: %d --- %s\n',ME.stack(errorstack_i).line,ME.stack(errorstack_i).name),
                end
                if retryflag % retry once with default thresholds
                    new_th = thresholds;
                    searchflag = true;
                    fprintf('Error occured, retrying with default thresholds...\n');
                    retryflag = false;
                else % fail on second try -> skip entirely
                    new_th = thresholds;
                    i = i + 1; % skip this pulse
                    fprintf('Error occured, skipping this file...\n');
                    retryflag = true;
                end
            end
            
            if i > size(uni_win,1) % return if finished
                fprintf('Done reading, check output!\n');
                close all,
                return
            end
        end
        
        % prompt for input to determine next operation
%         fprintf(['Is anything wrong? Enter P to save this data\n\n' ... 
%                     'ENTER will skip this pulse,\n' ... 
%                     'T will adjust top threshold only\n' ...
%                     'Y will adjust bottom threshold only\n' ...
%                     '+ or = will auto-adjust according to top threshold' ...
%                     '\nOr, just enter a new top threshold to\nautomatically set both\n' ...
%                     ]);
        fprintf('\n%d cell(s) saved\n',good_index-1);
        fprintf('\nCurrent thresholds: %3.2e, %3.2e\n',new_th);
        OK = input('---\n','s');
        
        %% input case structures
        switch OK
            case [] % empty
            
            fprintf('Skipping this pulse...\n');
            i = i+1;
            new_th = thresholds;
            
        % save data, auto-choose window
            case {'P','p', '.', '/'}
            
%             fprintf('Ok, displaying windows and file numbers...\n'),
            
            % make and display table of pulses and indices
            indices = iter_out(:,1); winds = 1:length(indices);
            table_data = [winds', indices];
            
            % clean up empty table entries
            cci = 1;
            stopc = size(table_data,1);
            while(cci <= stopc)
                if table_data(cci,2) == 0
                    table_data(cci,:) = [];
                    stopc = stopc - 1;
                else
                    cci = cci + 1;
                end
            end
                        
            
            if ~isempty(table_data) % make sure WindowTable is NOT empty                        
    %             iter_out_index = input('Select window to save:\n---\n');
     
                WindowTable = array2table(table_data,'VariableNames',{'Window','StartIndex'}),
                
                iter_out_index = 1;
                fprintf('Ok, saving data...\n');

                if (0 < iter_out_index) && (iter_out_index <= max(winds))

                    output_matrix(good_index,:) = iter_out(iter_out_index,:); % save to output matrix
                    output_matrix(good_index,1) = output_matrix(good_index,1) + uni_win(i) - 200;
                    good_index = good_index + 1;
                    i = i+1;
                    new_th = thresholds; % reset thresholds
                else
                    fprintf('Unrecognized input.\n');
                    beep,
                end      
                
            else
                fprintf('Pulse table is empty; force retry\n');
                new_th = thresholds; % reset thresholds
            end
            
            % save data, user picks window
            case {'//'}
            
%             fprintf('Ok, displaying windows and file numbers...\n'),
            
            % make and display table of pulses and indices
            indices = iter_out(:,1); winds = 1:length(indices);
            table_data = [winds', indices];
            
            % clean up empty table entries
            cci = 1;
            stopc = size(table_data,1);
            while(cci <= stopc)
                if table_data(cci,2) == 0
                    table_data(cci,:) = [];
                    stopc = stopc - 1;
                else
                    cci = cci + 1;
                end
            end
            
            if ~isempty(table_data) % make sure WindowTable is NOT empty
                        
                WindowTable = array2table(table_data,'VariableNames',{'Window','StartIndex'}),                        
                iter_out_index = input('Select window to save:\n---\n');
                fprintf('Ok, saving data...\n');

                if (0 < iter_out_index) && (iter_out_index <= max(winds))

                    output_matrix(good_index,:) = iter_out(iter_out_index,:); % save to output matrix
                    output_matrix(good_index,1) = output_matrix(good_index,1) + uni_win(i) - 200;
                    good_index = good_index + 1;
                    i = i+1;
                    new_th = thresholds; % reset thresholds
                else
                    fprintf('Unrecognized input.\n');
                    beep,
                end      
                
            else
                fprintf('Pulse table is empty; force retry\n');
                new_th = thresholds; % reset thresholds
            end    
            
        % in case of bad threshold, ask for new thresholds and retry   
            case {'T','t'}
            
            th_input = input('Please input new top threshold\n---\n');
            
                if isempty(th_input)
                    fprintf('Unrecognized input: stopping and dumping data.\n');
                    break

                else
                    new_th(2) = th_input;
                    fprintf('New top threshold: %3.2e\n',new_th(2));
                end
            
            case {'Y','y'}
            
            th_input = input('Please input new bottom threshold\n---\n');
            
                if isempty(th_input)
                    fprintf('Unrecognized input: stopping and dumping data.\n');
                    break

                else
                    new_th(1) = th_input;
                    fprintf('New bottom threshold: %3.2e\n',new_th(1));
                end
            
                % automatically compute threshold values based on max vals
                % from difference vector
            case {'+', '=', ''''}
                new_th(2) = auto*0.85;
                new_th(1) = new_th(2)*0.08;
                if new_th(1) < thresholds(1)
                    new_th(1) = thresholds(1);
                end
                fprintf('--------------\nSet new bottom threshold to %3.2e\n',new_th(1));
                fprintf('\nSet new upper threshold to %3.2e\n',new_th(2));
                
            otherwise
                if length(OK) > 1
                    th_input = str2double(OK);
                    new_th(1) = th_input/10;
                    new_th(2) = th_input;
                    fprintf('New bottom threshold: %3.2e\n',new_th(1));
                    fprintf('New top threshold: %3.2e\n',new_th(2));
            
                else
                    fprintf('Unrecognized input: stopping and dumping data.\n');
                    break
                end
        end
    end
    
    t = linspace(0,1,2^16); % time-samples
    Fs = 2^16; % sampling frequency
    y = 0.3*exp(-4*t).*sin(t*2*pi*440); % sinusoid
    sound(y,Fs),
    
    fprintf('Done reading, check output!\n');

end

