function [ OUT, empty, auto_thresh_value ] = mNPS_readKim( data_vector, sampleRate, ...
    thresholds, plotflag, fitflag )
% [ OUT, empty, auto_thresh_value ] = measureKimNPS( data_vector, sampleRate,
%     thresholds, plotflag, fitflag )
%   Reads mNPS data and returns OUT matrix

    %% SECTION 1: load data and perform basic signal conditioning
    
    Fs = sampleRate/1000; % convert to kHz
    
    N = 20; % downsample factor
    
    % flip if negative
    if (mean(data_vector) < 0 )
        data_vector = -data_vector;
    end
    
    y_smoothed = mNPS_fastsmooth(data_vector',200,1,1); % perform rectangular smoothing
    
    ym = downsample(y_smoothed,N); % downsample by period N
    
    if size(ym,1) < size(ym,2)
        ym = ym'; % transpose if vector is of the wrong dimension
    end
    
    y_detrend = ym - ASLS(ym,1e9,1e-4,10); % remove trend
    
    %% SECTION 2: threshold signal by differences
    % take the difference of ym, threshold by lower value
    % thresholds for differences, user provided
    
    ym_diff = diff(y_detrend); % compute difference
    ym_diff(abs(ym_diff) < thresholds(1)) = 0; % threshold values below thresholds(1)
    ym_diff(1:10) = 0; % zero out the first few values
    ym_diff(end-10:end) = 0; % zero out the last few values

    % ensure all values in squeeze channel are zero
    squeeze_begin = find(ym_diff <= -thresholds(2),1); % find where squeeze starts

    % advance by 1 until large difference ends
    while (ym_diff(squeeze_begin) <= -thresholds(2))
        squeeze_begin = squeeze_begin + 1;
    end

    squeeze_end = find(ym_diff(squeeze_begin+20:end) >= thresholds(2),1); % find end of squeeze
    ym_diff(squeeze_begin+1:squeeze_begin+squeeze_end+20) = 0; % zero out the squeeze channel

    %% SECTION 3: identify nonzero differences, nz_mat is the matrix of nonzero differences

    nz_mat = ones(2,length(nonzeros(ym_diff))); % preallocation
    k = 1; % index for A
    for i = 1:length(ym_diff)
        if (ym_diff(i) ~= 0) % look for nonzeros
            nz_mat(1,k) = i; % array index of nonzero
            nz_mat(2,k) = ym_diff(i); % nonzero value
            k = k+1;
        end
    end
    

    %% SECTION 4: remove error from A
    
    k=1;
    while (k < ceil(log(length(nz_mat))))
        i=1;
        while (i < length(nz_mat))

            % Case 1: current and next both positive && next > current
            if nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                    nz_mat(2,i+1) > nz_mat(2,i)
                % move next into current
                nz_mat(:,i) = nz_mat(:,i+1);

            % Case 2: current and next both positive && current > next
            elseif nz_mat(2,i) > 0 && nz_mat(2,i+1) > 0 && ...
                    nz_mat(2,i) > nz_mat(2,i+1)
                % move current into next
                nz_mat(:,i+1) = nz_mat(:,i);

            % Case 3: current and next both negative && current > next
            elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                    nz_mat(2,i) > nz_mat(2,i+1)
                % move next into current
                nz_mat(:,i) = nz_mat(:,i+1);

            % Case 4: current and next both negative && next > current
            elseif nz_mat(2,i) < 0 && nz_mat(2,i+1) < 0 && ...
                    nz_mat(2,i+1) > nz_mat(2,i)
                % move current into next
                nz_mat(:,i+1) = nz_mat(:,i);
            end

            i=i+1;
        end
        k=k+1;
    end
    
    %% SECTION 5: remove repeats in A
    
    unique_xs = unique(nz_mat(1,:)); % unique values
    unique_is = ones(1,length(unique_xs)); % unique indices
    for i = 1:length(unique_is)
        unique_is(i) = find(nz_mat(1,:) == unique_xs(i), 1);
    end
    unique_ys = nz_mat(2,unique_is);
    nz_mat = [unique_xs; unique_ys];
   

    %% SECTION 6: rectangularize pulses
    
    ym_rect = y_detrend;
    k = 1;
    while (k <= 50)
        i = 1;
        while (i < length(nz_mat))

            if nz_mat(2,i) < 0 && nz_mat(2,i+1) > 0 % look for sign change in differences
                ym_rect(nz_mat(1,i):nz_mat(1,i+1)) = ...
                    mean(y_detrend(nz_mat(1,i):nz_mat(1,i+1))); 
                % replace all values in between with mean
            end

            i=i+1;
        end
        k=k+1;
    end

    %% SECTION 7: Plot figures if flag is true
    
    if max(ym_diff) < thresholds(1) % waste of time, no pulse
        plotflag = false;
        empty = true;
    else
        empty = false;
    end

    if plotflag == true
        Pix_SS = get(0,'screensize');
        figh = figure(42); 
        figsize = [0.1 0.1 0.45 0.75]*Pix_SS(4);
        set(figh,'units','pixels','pos',figsize);
        
        % take top and bottom 3 values
        nsorted_d = sort(ym_diff);
        min_vals = nsorted_d(1:3);
        
        psorted_d = sort(ym_diff,'descend');
        max_vals = psorted_d(1:3);
        
        % set auto-thresholds
        if abs(min_vals(3)) < abs(max_vals(3))
            auto_thresh_value = abs(min_vals(3));
        else
            auto_thresh_value = abs(max_vals(3));
        end
        
        % difference plot
        ax1 = subplot(3,1,1);
        figwin_tighten();
        difp = plot(ym_diff,'k-'); difp.LineWidth = 1;
        title('y_{diff}'), set(gca,'FontSize',10),
        grid(ax1,'on'), ax1.XMinorGrid = 'on';
        axis([0, length(ym_diff), 1.1*mean(min_vals), 1.1*mean(max_vals)]),
        for i = 1:length(max_vals)
            label_str = sprintf('%3.3e',min_vals(i));
            text(i*600,1.35*max_vals(1),label_str,'FontSize',12);
            
            label_str = sprintf('%3.3e',max_vals(i));
            text(i*600,1.85*max_vals(1),label_str,'FontSize',12);
        end
        
        % plot thresholds
        hold on,
        linel = length(ym_diff);
        bthlin_l = line([0 linel], [-thresholds(1), -thresholds(1)]);
        bthlin_u = line([0 linel], [thresholds(1), thresholds(1)]);
        
        tthlin_l = line([0 linel], [-thresholds(2), -thresholds(2)]);
        tthlin_u = line([0 linel], [thresholds(2), thresholds(2)]);
        
        bthlin_l.Color = [1 0 0]; bthlin_u.Color = [1 0 0];
        tthlin_l.Color = [0 0 1]; tthlin_u.Color = [0 0 1];
        hold off,
        
        % rectangularized
        ax2 = subplot(3,1,2);
        figwin_tighten();
        plot(ym_rect), title('y_{rect}'), set(gca,'FontSize',10), 
        grid(ax2,'on'), ax2.XMinorGrid = 'on';
        axis([0, length(ym_rect), 1.1*min(ym_rect), 0.01]),
        
        % smoothed
        subplot(3,1,3),
        figwin_tighten();
        plot(y_smoothed,'k-'), title('y_{LP}'), set(gca,'FontSize',10),
        axis([0, length(y_smoothed), 0.999*min(y_smoothed), 1.001*max(y_smoothed)]),
                
    end
    
    %% SECTION 8: Detect NPS pulses
    %       pulse_series is a matrix with the indices and parameters for
    %       rectangular pulses
    
    i=1;
    k = 0;
    backset = 10;
    pulse_series = ones(length(nz_mat),4);
    while (i < length(nz_mat))
        if nz_mat(2,i) < 0 && nz_mat(2,i+1) > 0 % starts negative and flips sign
            k = k + 1;
            pulse_series(k,1) = nz_mat(1,i);  % Start index
            pulse_series(k,2) = nz_mat(1,i+1); % End index
            pulse_series(k,3) = mean(ym((nz_mat(1,i)-backset):(nz_mat(1,i)-backset+10))); % normalized baseline current
            pulse_series(k,4) = mean(y_detrend(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % avg current drop between pulses
            pulse_series(k,5) = std(y_detrend(nz_mat(1,i)+1:nz_mat(1,i+1)-1)); % std dev of current drop
        end
        i=i+1;
    end
    
    % remove empty entries in P
    cci = 1;
    stopc = size(pulse_series,1);
    while(cci <= stopc)
        if (pulse_series(cci,1) == 1 && pulse_series(cci,2) == 1)
            pulse_series(cci,:) = [];
            stopc = stopc - 1;
        else
            cci = cci + 1;
        end
    end
    
    %% SECTION 9: Extract Data for mNPS-R
    
    total_segs = 6;
    out = ones(length(pulse_series) - (total_segs - 1),11);
    for k = 1:length(pulse_series)+1 - total_segs
        start_index = pulse_series(k,1); % starting index
        I = pulse_series(k,3); % baseline current
        dI = -mean([pulse_series(k,4), pulse_series(k+1,4)]); % node-pore current drop
        dIC = -pulse_series(k+2,4); % squeeze current drop
        dICstd = pulse_series(k+2,5); % squeeze current drop std. dev
        dT = mean( [    (pulse_series(k,2) - pulse_series(k,1))/Fs*N, ...
                  (pulse_series(k+1,2) - pulse_series(k+1,1))/Fs*N]); % NP transit time (ms)
        dTsq = (pulse_series(k+2,2) - pulse_series(k+2,1))/Fs*N; % squeeze transit time (ms)
        
        % post-squeeze NP current drops
        dI4 = -pulse_series(k+3,4);
        dI5 = -pulse_series(k+4,4);
        dI6 = -pulse_series(k+5,4);
        
        % recovery time determined when post-squeeze NP current drop reaches
        % pre-squeeze NP current drop
        
        if abs(dI6-dI)/dI < 0.08
            Tr = (pulse_series(k+4,2)-pulse_series(k+2,1))/Fs*N;
            
        elseif abs(dI5-dI)/dI < 0.08
            Tr = (pulse_series(k+3,2)-pulse_series(k+2,1))/Fs*N;
            
        elseif abs(dI4-dI)/dI < 0.08
            Tr = 0;
            
        else
            Tr = Inf;
        end
        
        if fitflag
            
            % populate vector rT with time-points of recovery pulses
            rT = [0, ...
                (pulse_series(k+3,2)-pulse_series(k+2,1))/Fs*N, ...
                (pulse_series(k+4,2)-pulse_series(k+2,1))/Fs*N];
            
            % populate vector rdI with current drop-amplitudes of recovery
            % pulses. Anticipate that rdI should increase
            %
            % NOTE: cannot approximate size for ellipsoid particle
            rdI = [abs(dI4), ...
                abs(dI5), ...
                abs(dI6)];
            
            % use MATLAB fit() to fit a linear polynomial to recovery data.
            % fo has fields pertaining to fit parameters:
            %   fo.p1 is the slope of the line, not sign-bounded
            %   fo.p2 is the y-intercept, not sign-bounded
            % gof has fields pertaining to goodness of fit, but only
            % gof.rsquare will be used
            [fo, gof] = fit(rT',rdI','poly1');
            
        else
            fo.p1 = 0;
            fo.p2 = 0;
            gof.rsquare = 0;
        end
            
        out(k,:) = [start_index, I, dI, dIC, dICstd, dT, dTsq, Tr, fo.p1, fo.p2, gof.rsquare];
        
    end
    
    if fitflag
        
        Pix_SS = get(0,'screensize');
        
        figf = figure(43);
        title('Curve Fit'),
        ffunc = @(x) fo.p2+fo.p1*x;
        fplot(ffunc,[rT(1), rT(end)]),
        hold on, scatter(rT,rdI), hold off,
        figsize = [0.70 0.6 1/3 2/9]*Pix_SS(4);
        set(figf,'units','pixels','pos',figsize);
        
    end


    %% Section 9b: export data
    
    L = 6750;
    Dnp = 27.65;
    npL = 700;
    sqL = 3000;
    H = 23;
    wC = 10.5;
    wNP = 22;
    DC = Dnp*(wC/wNP)^(0.5);
    
    calculated = zeros(size(out,1),8);
    
    % diameter
    calculated(:,1) = ((out(:,3)./out(:,2)*Dnp^2*L)./ ...
        (1+0.8*L/Dnp*out(:,3)./out(:,2))).^(1/3);
    
    % strain
    calculated(:,2) = (calculated(:,1)-wC)./calculated(:,1);
    
    % np velocity
    calculated(:,3) = npL./out(:,6);
    
    % sq velocity
    calculated(:,4) = sqL./out(:,7);
    
    % deformed diameter
    calculated(:,5) = 0.01*(pi/4*wC)*(((out(:,4)./out(:,2)*DC^2*L)./ ...
        (1+0.8*L/DC*out(:,4)./out(:,2))).^(1/3)).^2;
    
    % wCDI
    calculated(:,6) = calculated(:,4)./calculated(:,3) .* calculated(:,1)./H;
    
    % recovery time
    calculated(:,7) = out(:,8);
    
    % recovery rate
    calculated(:,8) = fo.p1;
    
    
    out = [out, calculated];
    
    OUT = out;

end

