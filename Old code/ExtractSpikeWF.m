function ExtractSpikeWF(spk_data, save_fn, varargin),

%Options
opts.Overwrite = 0;
opts.SampleFrequency = 25000;
opts.FilterSettleTime = 10; %seconds, how long to ignore at start/end of file
opts.SamplePoints = 32; %how many points to downsample the waveform to
opts.SpikeInterpFreq = 10^5; %in Hz, also used for spike alignment so keep high for good alignment
opts.DoAtOnce = 1;

%Spike threshold options
opts.SpikeExtract_SlidingWindow = 0;
opts.SpikeExtract_SlidingWindowTimeWidth = 15*60; %in seconds
opts.SpikeExtract_SlidingWindowTimeStep = 5*60; %in seconds
opts.SpikeExtract_STDSpikeThresh = -4;
opts.SpikeExtract_STDArtifactThresh = -7;
opts.SpikeExtract_WFRange = [-0.4 1.2]; %in ms

%Spike re-alignment options
opts.SpikeAlign = 1; %0 - no alignment; 1 - global min/max alignment; 2 - local min/max alignment
opts.SpikeAlignRange = [-0.1 0.4]; %in ms, what range to search for min/max

%Parse variable inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs.'); end
for i = 1:2:length(varargin),
    if ~isfield(opts, varargin{i}),
        error(sprintf('Key %s not valid.', varargin{i}));
    end
    opts.(varargin{i}) = varargin{i+1};
end

%Save file name
[pathstr, name, ext] = fileparts(save_fn);
if isempty(pathstr), pathstr = pwd; end
if isempty(ext), ext = '.wf'; end
save_fn = strcat(pathstr, '\', name, ext);
if exist(save_fn, 'file') && ~opts.Overwrite,
    fprintf('File %s already exists. Skipping.\n', save_fn);
    return;
end

%Is spk_sig a filename?
if ischar(spk_data),
    if ~exist(spk_data, 'file'),
        error(sprintf('Can''t find specified file: %s.', spk_data));
    end
    fprintf('Loading spike data from file...');
    load(spk_data, 'spk_data');
    fprintf('done.\n');
end

max_spk_ind = length(spk_data);
fprintf('Processing %d samples.\n', max_spk_ind);

%Ratio of interpolation to sampling frequency and rescaling points
interp_rescale_ratio = opts.SampleFrequency/opts.SpikeInterpFreq;
if ~isempty(opts.SamplePoints),
    sample_rescale_ratio = floor(diff(opts.SpikeExtract_WFRange)*opts.SampleFrequency/1000 + 1)/opts.SamplePoints;
else
    sample_rescale_ratio = 1;
end

%Initialize indexing variables
opts.SpikeExtract_WFRange = sort(opts.SpikeExtract_WFRange);
wf_range = round(opts.SpikeExtract_WFRange/1000*opts.SampleFrequency);
wf_range_ind = [wf_range(1):wf_range(2)];
%Number of points in waveform
if isempty(opts.SamplePoints),
    Nwf = length(wf_range_ind);
else
    Nwf = opts.SamplePoints;
end

%Was an alignment range specified?
opts.SpikeAlignRange = sort(opts.SpikeAlignRange);
if isempty(opts.SpikeAlignRange), opts.SpikeAlignRange = opts.SpikeExtract_WFRange; end

%Blank start/end
spk_data(1:opts.FilterSettleTime*opts.SampleFrequency) = NaN;
spk_data((end - opts.FilterSettleTime*opts.SampleFrequency):end) = NaN;

%Threshold spiking activity
fprintf('Thresholding spiking data:\n');
if ~opts.SpikeExtract_SlidingWindow,
    fprintf('\tUsing grand average.\n');
    spk_thresh = nanmean(spk_data) + opts.SpikeExtract_STDSpikeThresh*nanstd(spk_data);
    art_thresh = nanmean(spk_data) + opts.SpikeExtract_STDArtifactThresh*nanstd(spk_data);
else
    opts.SpikeExtract_SlidingWindowTimeStep = round(opts.SpikeExtract_SlidingWindowTimeStep*opts.SampleFrequency);
    opts.SpikeExtract_SlidingWindowTimeWidth = round(opts.SpikeExtract_SlidingWindowTimeWidth*opts.SampleFrequency);
    fprintf('\tUsing sliding window: %4.0f/%4.0f', 0, ceil(max_spk_ind./opts.SpikeExtract_SlidingWindowTimeStep));
    
    cur_t = 1; cur_trange_ind = [1 1]; count = 0;
    spk_thresh = zeros(size(spk_data)); art_thresh = spk_thresh;
    while (cur_t < max_spk_ind),
        cur_trange = cur_t + [-0.5 0.5]*opts.SpikeExtract_SlidingWindowTimeWidth;
        cur_trange(1) = max(1, cur_trange(1)); cur_trange(2) = min(cur_trange(2), max_spk_ind);
        spk_thresh(cur_trange(1):cur_trange(2)) = nanmean(spk_data(cur_trange(1):cur_trange(2))) + ...
            opts.SpikeExtract_STDSpikeThresh*nanstd(spk_data(cur_trange(1):cur_trange(2)));
        art_thresh(cur_trange(1):cur_trange(2)) = nanmean(spk_data(cur_trange(1):cur_trange(2))) + ...
            opts.SpikeExtract_STDArtifactThresh*nanstd(spk_data(cur_trange(1):cur_trange(2)));
        
        cur_t = cur_t + opts.SpikeExtract_SlidingWindowTimeStep;
        count = count + 1;
        fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', count, ceil(max_spk_ind./opts.SpikeExtract_SlidingWindowTimeStep));
    end
    fprintf('\n');
end

%Are we searching for above threshold?  Then invert everything
if (opts.SpikeExtract_STDSpikeThresh > 0),
    spk_data = -spk_data;
    spk_thresh = -spk_thresh;
    art_thresh = -art_thresh;
end

%Threshold our data
fprintf('\tThresholding.\n');
spk_thresh = (spk_data <= spk_thresh) & (spk_data >= art_thresh);

%Find where it crosses (and doesn't go back)
fprintf('\tFinding spikes...');
spk_ind = find(spk_thresh & (diff([0; spk_thresh]) == 1));
clear('spk_thresh');
fprintf('found %d putative spikes.\n', length(spk_ind));

%Remove 'duplicates' where they are within a spike window of each other
fprintf('\tRemoving spikes within an extract window of each other...');
bad_spk_ind = false(length(spk_ind), 1); cur_good = 1; min_dist = (opts.SpikeExtract_WFRange(2)/1000*opts.SampleFrequency);
for i = 2:length(spk_ind),
    if (spk_ind(i) - spk_ind(cur_good)) <= min_dist,
        bad_spk_ind(i) = 1;
    else
        cur_good = i;
    end
end
fprintf('removed %d spikes.\n', sum(bad_spk_ind));
spk_ind = spk_ind(~bad_spk_ind);
fprintf('\tFound %d spikes (%4.2f Hz).\n', length(spk_ind), length(spk_ind)./(max_spk_ind./opts.SampleFrequency));

%Open output file, write header information
fid = fopen(save_fn, 'w');
fwrite(fid, 'EXTRACTEDSPIKEWF', 'char*1');
fwrite(fid, length(spk_ind), 'int64');
fwrite(fid, opts.SampleFrequency, 'double');
fwrite(fid, opts.FilterSettleTime, 'double');
fwrite(fid, opts.SamplePoints, 'double');
fwrite(fid, opts.SpikeInterpFreq, 'double');
fwrite(fid, opts.SpikeExtract_SlidingWindow, 'uint8');
fwrite(fid, opts.SpikeExtract_SlidingWindowTimeWidth, 'double');
fwrite(fid, opts.SpikeExtract_SlidingWindowTimeStep, 'double');
fwrite(fid, opts.SpikeExtract_STDSpikeThresh, 'double');
fwrite(fid, opts.SpikeExtract_STDArtifactThresh, 'double');
fwrite(fid, opts.SpikeExtract_WFRange, 'double');
fwrite(fid, opts.SpikeAlign, 'uint8');
fwrite(fid, opts.SpikeAlignRange, 'double');
fwrite(fid, Nwf, 'int32');
    
%Do we do everything at once?  Takes more memory but faster and does better alignment to zero
if opts.DoAtOnce,
    
    if (opts.SpikeAlign > 0),
        %Add to front/back for alignment shifts
        orig_wf_range = wf_range;
        wf_range(1) = opts.SpikeExtract_WFRange(1) + min(0, opts.SpikeAlignRange(1));
        wf_range(2) = opts.SpikeExtract_WFRange(2) + max(0, opts.SpikeAlignRange(2));
        wf_range = round(wf_range/1000*opts.SampleFrequency);
        wf_range_ind = [wf_range(1):wf_range(2)];
    end
        
    %Do we need to center interpolation index on 0?
    if (wf_range(1) < 0) && (wf_range(2) > 0),
        %Center on zero
        pre_wf_interp_ind = -[0:interp_rescale_ratio:(-wf_range(1))];
        post_wf_interp_ind = [0:interp_rescale_ratio:wf_range(2)];
        wf_interp_ind = cat(2, pre_wf_interp_ind(end:-1:2), post_wf_interp_ind);
    else
        %Just fill the range from the left side
        wf_interp_ind = [wf_range(1):interp_rescale_ratio:wf_range(2)];
    end

    spk_wf = zeros(length(wf_range_ind), length(spk_ind));
    for cur_spk = 1:length(spk_ind),
        % Grab the current waveform (plus more for alignment purposes if needed)
        cur_wf_ind = wf_range_ind + spk_ind(cur_spk);
        valid_cur_wf_ind = ((cur_wf_ind > 0) & (cur_wf_ind <= max_spk_ind));
        spk_wf(valid_cur_wf_ind, cur_spk) = spk_data(cur_wf_ind(valid_cur_wf_ind));
    end
    %Do spline fit on our local pieces
    spline_wf = spline(wf_range_ind, spk_wf', wf_interp_ind)';
    
    if (opts.SpikeAlign > 0),
        %Align spikes on a limited range (if specified)
        fprintf('\tRe-aligning to peak...');
        if ~isempty(opts.SpikeAlignRange),
            align_range_ind = find((wf_interp_ind >= (opts.SpikeAlignRange(1)/1000*opts.SampleFrequency)) & (wf_interp_ind <= (opts.SpikeAlignRange(2)/1000*opts.SampleFrequency)));
        else
            align_range_ind = find((wf_interp_ind >= (opts.SpikeExtract_WFRange(1)/1000*opts.SampleFrequency)) & (wf_interp_ind <= (opts.SpikeExtract_WFRange(2)/1000*opts.SampleFrequency)));
        end
                    
        %Find minima
        mi = align_range_ind(1)*ones(size(spline_wf, 2), 1);
        if opts.SpikeAlign == 1,
            %Find global minima
            fprintf('finding global minima...');
            [~, mi] = min(spline_wf(align_range_ind, :), [], 1);
        elseif opts.SpikeAling == 2,
            %Find first local minima
            fprintf('finding first local minima...');
            for cur_spk = 1:length(spk_ind),
                mi(cur_spk) = find(diff(diff(spline_wf(align_range_ind, cur_spk), 1, 1) >= 0, 1, 1) == 1, 1, 'first') + 1;
            end
        end
        spk_ind = spk_ind + wf_interp_ind(align_range_ind(mi))';
        fprintf('done.\n');
        
        %Shift by our alignment
        if (opts.SpikeExtract_WFRange(1) < 0) && (opts.SpikeExtract_WFRange(2) > 0),
            %Do we need to center interpolation index on 0?
            %Center on zero
            pre_orig_wf_interp_ind = -[0:interp_rescale_ratio:(-orig_wf_range(1))];
            post_orig_wf_interp_ind = [0:interp_rescale_ratio:orig_wf_range(2)];
            orig_wf_interp_ind = cat(2, pre_orig_wf_interp_ind(end:-1:2), post_orig_wf_interp_ind);
            
            %Rescale index (into interpolated data)
            pre_rescale_ind = length(pre_orig_wf_interp_ind) - [0:(sample_rescale_ratio/interp_rescale_ratio):length(pre_orig_wf_interp_ind)];
            post_rescale_ind = (length(pre_orig_wf_interp_ind)-1) + [0:(sample_rescale_ratio/interp_rescale_ratio):length(post_orig_wf_interp_ind)];
            rescale_ind = cat(2, pre_rescale_ind(end:-1:2), post_rescale_ind);
            rescale_ind = round(rescale_ind+1);
        else
            %Just fill the range from the left side
            orig_wf_interp_ind = [orig_wf_range(1):interp_rescale_ratio:orig_wf_range(2)];
            rescale_ind = round([0:(Nwf-1)]/interp_rescale_ratio/sample_rescale_ratio)+1;
        end
        
        for cur_spk = 1:length(spk_ind),
            spline_wf(1:length(orig_wf_interp_ind), cur_spk) = spline_wf(align_range_ind(mi(cur_spk)) + orig_wf_interp_ind/interp_rescale_ratio, cur_spk);
        end
        spline_wf = spline_wf(1:length(orig_wf_interp_ind), :);
    end
        
    %Rescale
    spline_wf = spline_wf(rescale_ind, :);
    if (opts.SpikeExtract_STDSpikeThresh > 0),
        spline_wf = -spline_wf; %Invert back for finding peaks
    end
    
    fprintf('\tWriting to file...');
    %For each spike index grab its time, write to file
    fwrite(fid, (spk_ind - 1)/opts.SampleFrequency, 'double');
    %Write to disk (with zeros for cluster marker)
    fwrite(fid, cat(1, spline_wf, zeros(1, size(spline_wf, 2))), 'double');
else
    %Do we need to center interpolation index on 0?
    if (opts.SpikeExtract_WFRange(1) < 0) && (opts.SpikeExtract_WFRange(2) > 0),
        %Center on zero
        pre_wf_interp_ind = -[0:interp_rescale_ratio:(-wf_range(1))];
        post_wf_interp_ind = [0:interp_rescale_ratio:wf_range(2)];
        wf_interp_ind = cat(2, pre_wf_interp_ind(end:-1:2), post_wf_interp_ind);
        
        %Rescale index (into interpolated data)
        pre_rescale_ind = length(pre_wf_interp_ind) - [0:(sample_rescale_ratio/interp_rescale_ratio):length(pre_wf_interp_ind)];
        post_rescale_ind = length(pre_wf_interp_ind) + [0:(sample_rescale_ratio/interp_rescale_ratio):length(post_wf_interp_ind)];
        rescale_ind = cat(2, pre_rescale_ind(end:-1:2), post_rescale_ind);
        rescale_ind = round(rescale_ind + 1);
    else
        %Just fill the range from the left side
        wf_interp_ind = [wf_range(1):interp_rescale_ratio:wf_range(2)];
        rescale_ind = round([0:(Nwf-1)]/interp_rescale_ratio/sample_rescale_ratio)+1;
    end

    if (opts.SpikeAlign > 0),
        %If aligning, need indexing variables for alignment purposes
        align_range_ind = round(opts.SpikeAlignRange/1000*opts.SampleFrequency);
        if (opts.SpikeAlignRange(1) < 0) && (opts.SpikeAlignRange(2) > 0),
            align_interp_ind = -[0:interp_rescale_ratio:(-align_range_ind(1))];
            align_interp_ind = cat(2, align_interp_ind(end:-1:2), [0:interp_rescale_ratio:align_range_ind(2)]);
        else
            align_interp_ind = [align_range_ind(1):interp_rescale_ratio:align_range_ind(2)];
        end
        align_range_ind = [align_range_ind(1):align_range_ind(2)];
        
        fprintf('\tRe-aligning to peak...');
        align_wf = zeros(length(align_range_ind), length(spk_ind));
        for cur_spk = 1:length(spk_ind),
            %Grab the portion around the threshold to search for min/max
            cur_wf_ind = align_range_ind + spk_ind(cur_spk);
            valid_cur_wf_ind = ((cur_wf_ind > 0) & (cur_wf_ind <= max_spk_ind));
            align_wf(valid_cur_wf_ind, cur_spk) = spk_data(cur_wf_ind(valid_cur_wf_ind));
        end
        %Do spline fit on our local pieces
        spline_wf = spline(align_range_ind, align_wf', align_interp_ind)';
        
        %Find minima
        if opts.SpikeAlign == 1,
            %Find global minima
            fprintf('finding global minima...');
            [~, mi] = min(spline_wf, [], 1);
        elseif opts.SpikeAling == 2,
            %Find first local minima
            fprintf('finding first local minima...');
            for cur_spk = 1:length(spk_ind),
                mi(cur_spk) = find(diff(diff(spline_wf(:, cur_spk), 1, 1) >= 0, 1, 1) == 1, 1, 'first') + 1;
            end
        end
                
        fprintf('done.\n');
    end
    
    %For each spike index grab its time, write to file
    fwrite(fid, (round(spk_ind) - 1)/opts.SampleFrequency, 'double');
    
    fprintf('\tWriting to file...');
    %For each spike, grab a waveform, write to file
    for cur_spk = 1:length(spk_ind),
        
        % Grab the current waveform
        cur_wf_ind = round(wf_range_ind + spk_ind(cur_spk));
        cur_interp_ind = wf_interp_ind + spk_ind(cur_spk);
        valid_cur_wf_ind = ((cur_wf_ind > 0) & (cur_wf_ind <= max_spk_ind));
        
        %Do spline fit on our local piece
        spline_wf = spline(cur_wf_ind(valid_cur_wf_ind), spk_data(cur_wf_ind(valid_cur_wf_ind)), cur_interp_ind);
        spline_wf((cur_interp_ind <= 0) | (cur_interp_ind > max_spk_ind)) = 0;
        
        %Use rescaling to change #samples and write to disk
        if (opts.SpikeExtract_STDSpikeThresh > 0),
            fwrite(fid, -spline_wf(rescale_ind), 'double'); %waveform data
        else
            fwrite(fid, spline_wf(rescale_ind), 'double'); %waveform data
        end
        fwrite(fid, 0, 'double'); %cluster marker
    end
end
fprintf('done.\n');

%Close file
fclose(fid);

