function ProcessIntan(intan_fn, varargin)

opts.DataOverwrite = 0;
opts.SampleFrequency = 25000;
opts.Interactive = 1;
opts.SaveRawData = 1;
%Notch filter options
opts.NotchFrequencies = [];
opts.NotchOrder = 16;
opts.NotchPassbandRipple = 0.1;
opts.NotchStopbandAttenuation = 60;
opts.NotchQuality = 100; %Q = F0/BW where BW is the bandwidth of the notch
%LFP Filter options
opts.doLFPFilt = 0;
opts.LFPFilt_FStop1 = 0.5; opts.LFPFilt_FPass1 = 1; opts.LFPFilt_AStop1 = 40;
opts.LFPFilt_FStop2 = 150; opts.LFPFilt_FPass2 = 125; opts.LFPFilt_AStop2 = 40;
opts.LFPFilt_APass = 1;
%Spike Filter options
opts.doSpikeFilt = 1;
opts.SpikeFilt_FStop1 = 150; opts.SpikeFilt_FPass1 = 250; opts.SpikeFilt_AStop1 = 40;
opts.SpikeFilt_FStop2 = 8000; opts.SpikeFilt_FPass2 = 6000; opts.SpikeFilt_AStop2 = 40;
opts.SpikeFilt_APass = 1;
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs.'); end
for i = 1:2:length(varargin),
    if ~isfield(opts, varargin{i}),
        error(sprintf('Key %s not valid.', varargin{i}));
    end
    opts.(varargin{i}) = varargin{i+1};
end

%Create LFP filter
match  = 'stopband';  % Band to match exactly
% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(opts.LFPFilt_FStop1, opts.LFPFilt_FPass1, opts.LFPFilt_FPass2, opts.LFPFilt_FStop2, ...
    opts.LFPFilt_AStop1, opts.LFPFilt_APass, opts.LFPFilt_AStop2, opts.SampleFrequency);
LFP_Hd = design(h, 'cheby2', 'MatchExactly', match);

%Create Spike filter
match  = 'passband';  % Band to match exactly
% Construct an FDESIGN object and call its CHEBY2 method.
h  = fdesign.bandpass(opts.SpikeFilt_FStop1, opts.SpikeFilt_FPass1, opts.SpikeFilt_FPass2, opts.SpikeFilt_FStop2, ...
    opts.SpikeFilt_AStop1, opts.SpikeFilt_APass, opts.SpikeFilt_AStop2, opts.SampleFrequency);
Spike_Hd = design(h, 'cheby2', 'MatchExactly', match);

if ~iscell(intan_fn),
    if ~isstr(intan_fn),
        error('Must pass a filename (or matching string) or cell array of filenames.');
    end
    if ~exist(intan_fn, 'file'),
        temp_files = dir(intan_fn);
        if isempty(temp_files),
            error('Couldn''t find any files matching: %s.\n', intan_fn);
        end
        intan_fn_list = struct2cell(temp_files);
        intan_fn_list = intan_fn_list(1, :);
    else
        intan_fn_list = {intan_fn};
    end
else
    intan_fn_list = intan_fn;
end

orig_NotchFrequencies = opts.NotchFrequencies;
for file_ind = 1:length(intan_fn_list),
    intan_fn = intan_fn_list{file_ind};
    
    %Save file name
    [pathstr, name, ext] = fileparts(intan_fn);
    if isempty(pathstr), pathstr = pwd; end
    save_fn = strcat(pathstr, '\', name);
    
    %Read header information from intan file
    [num_chans, chan_list] = read_intan_header(intan_fn);
    
    %Save encode information
    save_filename = sprintf('%s_Encodes.mat', save_fn);
    if ~exist(save_filename, 'file') || opts.DataOverwrite,
        fprintf('Loading encodes...\n');
        [~, ~, aux, t] = read_intan_data(intan_fn, 'ReadChannels', chan_list(1));
        fprintf('Saving to MAT...');
        save(save_filename, 't', 'aux', 'opts', '-v7.3');
        fprintf('done.\n');
    else
        fprintf('File %s, Encodes already processed. Skipping.\n', save_fn);
    end
    
    %Interactive filtering?
    opts.NotchFrequencies = orig_NotchFrequencies;
    if opts.Interactive,
        fprintf('Loading data for interactive filtering of noise...\n');
        
        %Loop through channels, filter LFP and Spike
        Pxx = []; clear Hpsd;
        for cur_chan = 1:length(chan_list),
            save_filename = sprintf('%s_Ch%d.mat', save_fn, chan_list(cur_chan));
            if exist(save_filename, 'file') && ~opts.DataOverwrite,
                fprintf('\tFile %s, Channel %d already processed. Skipping.\n', save_fn, chan_list(cur_chan));
                continue;
            end
            
            %Load channel data
            [~, chan_data, ~, t] = read_intan_data(intan_fn, 'ReadChannels', chan_list(cur_chan));
            chan_data = double(chan_data);
            
            nfft = 2^min(21, nextpow2(length(chan_data)));
            Pxx(:, cur_chan) = abs(fft(chan_data, nfft)).^2/length(chan_data)/opts.SampleFrequency;
            Hpsd(cur_chan) = dspdata.psd(Pxx(1:size(Pxx, 1)/2, cur_chan),'Fs',opts.SampleFrequency);
        end
        if isempty(Pxx), continue; end
        avgHpsd = dspdata.psd(nanmean(Pxx(1:size(Pxx, 1)/2, :), 2),'Fs',opts.SampleFrequency);
        
        freq_ind = (Hpsd(1).Frequencies < 400);
        stop_filtering = 0;
        while ~stop_filtering,
            figure;
            for i = 1:length(Hpsd),
                plot(Hpsd(i).Frequencies(freq_ind), 10*log10(Hpsd(i).Data(freq_ind))); hold all;
            end
            plot(avgHpsd.Frequencies(freq_ind), 10*log10(avgHpsd.Data(freq_ind)), 'k-', 'LineWidth', 1.5);
            xlabel('Frequency (Hz'); ylabel('Power (dB)');
            title('Select peaks to notch.  Click outside figure to quit.');
            for i = 1:length(opts.NotchFrequencies),
                [~, mi] = min(abs(avgHpsd.Frequencies - opts.NotchFrequencies(i)));
                plot(avgHpsd.Frequencies(mi), 10*log10(avgHpsd.Data(mi)), 'r*', 'MarkerSize', 12);
            end
            v = axis; orig_v = v;
            button = 3; zoom_count = 0; max_zoom_count = 5;
            while (button ~= 1),
                [x,y,button] = ginput(1);
                if (button == 3),
                    v = axis;
                    v(1:2) = x + [-0.25 0.25]*diff(v(1:2));
                    axis(v);
                    zoom_count = zoom_count + 1;
                    if zoom_count > max_zoom_count,
                        zoom_count = 0;
                        axis(orig_v);
                    end
                elseif (button == 2),
                    v = axis;
                    v(1:2) = x + [2 2]*diff(v(1:2));
                    axis(v);
                    zoom_count = max(0, zoom_count - 1);
                end
                drawnow;
            end
            close(gcf);
            if ~((x >= v(1)) && (x <= v(2)) && (y >= v(3)) && (y <= v(4))), stop_filtering = 1; continue; end
            
            %Find peak that was clicked
            [~,mi] = min(abs(avgHpsd.Frequencies - x));
            mi = mi + [-20:20]; mi = mi(mi > 0); mi = mi(mi <= length(avgHpsd.Frequencies));
            [pks, locs] = findpeaks(10*log10(avgHpsd.Data(mi)), 'MINPEAKDISTANCE', 4, 'SORTSTR', 'descend');
            good_ind = (abs(pks - y)/(v(4) - v(3)) < 0.05);
            locs = locs(good_ind);
            if isempty(locs), fprintf('ERROR: Couldn''t find a peak near the position clicked.  Try again.\n'); continue; end
            opts.NotchFrequencies = cat(1, opts.NotchFrequencies, avgHpsd.Frequencies(mi(locs(1))));
            opts.NotchFrequencies = unique(opts.NotchFrequencies);
        end
        drawnow;
    end
    
    %Loop through channels, filter LFP and Spike
    for cur_chan = 1:length(chan_list),
        
        save_filename = sprintf('%s_Ch%d.mat', save_fn, chan_list(cur_chan));
        if exist(save_filename, 'file') && ~opts.DataOverwrite,
            fprintf('File %s, Channel %d already processed. Skipping.\n', save_fn, chan_list(cur_chan));
            continue;
        end
        
        %Load channel data
        fprintf('Reading in channel %d...\n', chan_list(cur_chan));
        [~, chan_data, ~, ~] = read_intan_data(intan_fn, 'ReadChannels', chan_list(cur_chan));
        chan_data = double(chan_data);
        
        %Do notch filtering
        fprintf('Notch filtering at ['); fprintf('%4.1f ', opts.NotchFrequencies); fprintf('\b] Hz...');
        for cur_notch_f = 1:length(opts.NotchFrequencies),
            d = fdesign.notch('N,F0,Q,Ap,Ast', opts.NotchOrder, opts.NotchFrequencies(cur_notch_f), ...
                opts.NotchQuality, opts.NotchPassbandRipple, opts.NotchStopbandAttenuation, opts.SampleFrequency);
            notchHd(cur_notch_f) = design(d);
            chan_data = filtfilt(notchHd(cur_notch_f).sosMatrix, notchHd(cur_notch_f).ScaleValues, chan_data);
        end %notch freq loop
        fprintf('done.\n');
        
        %Save initial information to file
        save(save_filename, 'num_chans', 'cur_chan', 'chan_list', 'notchHd', 'opts', '-v7.3');
        
        %Filter channel data
        if opts.doSpikeFilt,
            fprintf('Filtering spikes...');
            spk_data = single(filtfilt(Spike_Hd.sosMatrix, Spike_Hd.ScaleValues, chan_data));
            fprintf('saving...');
            save(save_filename, 'spk_data', '-append', '-v7.3');
            fprintf('done.\n');
            clear('spk_data');
        end
        if opts.doLFPFilt,
            fprintf('Filtering LFPs...');
            lfp_data = single(filtfilt(LFP_Hd.sosMatrix, LFP_Hd.ScaleValues, chan_data));
            fprintf('saving...');
            save(save_filename, 'lfp_data', '-append', '-v7.3');
            fprintf('done.\n');
            clear('lfp_data');
        end
        if opts.SaveRawData,
            chan_data = single(chan_data);
            fprintf('Saving raw data...');
            save(save_filename, 'chan_data', '-append', '-v7.3');
            fprintf('done.\n');
        end
        
    end %channel loop
    
end %file loop