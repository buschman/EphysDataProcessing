function ExtractTrodeSpikeWF(spk_data, save_fn, varargin),

%Options
opts.Overwrite = 0;
opts.SampleFrequency = 25000;
opts.FilterSettleTime = 10; %seconds, how long to ignore at start/end of file
opts.SamplePoints = 32; %how many points to downsample the waveform to
opts.SpikeInterpFreq = 10^5; %in Hz, also used for spike alignment so keep high for good alignment

%Spike threshold options
opts.SpikeExtract_ThreshType = 0; %0 - STD, 1 - raw
opts.SpikeExtract_SlidingWindow = 0;
opts.SpikeExtract_SlidingWindowTimeWidth = 15*60; %in seconds
opts.SpikeExtract_SlidingWindowTimeStep = 5*60; %in seconds
opts.SpikeExtract_SpikeThresh = -4;
opts.SpikeExtract_ArtifactThresh = -7;
opts.SpikeExtract_WFRange = [-0.4 1.2]; %in ms

%Spike re-alignment options
opts.SpikeAlign = 1; %0 - no alignment; 1 - global min/max alignment; 2 - local min/max alignment
opts.SpikeAlignRange = [-0.1 0.4]; %in ms, what range to search for min/max

%Parse variable inputs
if (length(varargin) == 1) && iscell(varargin{1}),
    varargin = varargin{1}; %probably inherited from another function
end
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

max_spk_ind = size(spk_data, 1);
num_chans = size(spk_data, 2);
fprintf('Processing %d samples from %d channels.\n', max_spk_ind, num_chans);

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
    Nwf = length(wf_range_ind)*num_chans;
else
    Nwf = opts.SamplePoints*num_chans;
end

%Was an alignment range specified?
opts.SpikeAlignRange = sort(opts.SpikeAlignRange);
if isempty(opts.SpikeAlignRange), opts.SpikeAlignRange = opts.SpikeExtract_WFRange; end

%Blank start/end
spk_data(1:opts.FilterSettleTime*opts.SampleFrequency, :) = NaN;
spk_data((end - opts.FilterSettleTime*opts.SampleFrequency):end, :) = NaN;

%Threshold spiking activity
fprintf('Thresholding spiking data:\n');
if opts.SpikeExtract_ThreshType == 1,
    fprintf('\tUsing raw values.\n');
    %Use thresholds as raw value
    spk_thresh = opts.SpikeExtract_SpikeThresh;
    art_thresh = opts.SpikeExtract_ArtifactThresh;
elseif opts.SpikeExtract_ThreshType == 0,
    if ~opts.SpikeExtract_SlidingWindow,
        fprintf('\tUsing grand average.\n');
        temp_mean = nanmean(spk_data);
        temp_std = nanstd(spk_data);
        spk_thresh = repmat(temp_mean + opts.SpikeExtract_SpikeThresh*temp_std, [max_spk_ind 1]);
        art_thresh = repmat(temp_mean + opts.SpikeExtract_ArtifactThresh*temp_std, [max_spk_ind 1]);
    else
        opts.SpikeExtract_SlidingWindowTimeStep = round(opts.SpikeExtract_SlidingWindowTimeStep*opts.SampleFrequency);
        opts.SpikeExtract_SlidingWindowTimeWidth = round(opts.SpikeExtract_SlidingWindowTimeWidth*opts.SampleFrequency);
        fprintf('\tUsing sliding window: %4.0f/%4.0f', 0, ceil(max_spk_ind./opts.SpikeExtract_SlidingWindowTimeStep));
        
        cur_t = 1; cur_trange_ind = [1 1]; count = 0;
        spk_thresh = zeros(size(spk_data)); art_thresh = spk_thresh;
        while (cur_t < max_spk_ind),
            cur_trange = cur_t + [-0.5 0.5]*opts.SpikeExtract_SlidingWindowTimeWidth;
            cur_trange(1) = max(1, cur_trange(1)); cur_trange(2) = min(cur_trange(2), max_spk_ind);
            temp_mean = nanmean(spk_data(cur_trange(1):cur_trange(2), :));
            temp_std = nanstd(spk_data(cur_trange(1):cur_trange(2), :));
            
            spk_thresh(cur_trange(1):cur_trange(2), :) = repmat(temp_mean + opts.SpikeExtract_SpikeThresh*temp_std, [diff(cur_trange)+1 1]);
            art_thresh(cur_trange(1):cur_trange(2), :) = repmat(temp_mean + opts.SpikeExtract_ArtifactThresh*temp_std, [diff(cur_trange)+1 1]);
            
            cur_t = cur_t + opts.SpikeExtract_SlidingWindowTimeStep;
            count = count + 1;
            fprintf('\b\b\b\b\b\b\b\b\b%4.0f/%4.0f', count, ceil(max_spk_ind./opts.SpikeExtract_SlidingWindowTimeStep));
        end
        fprintf('\n');
    end
    clear('temp_mean', 'temp_std');
end

%Are we searching for above threshold?  Then invert everything
if (opts.SpikeExtract_SpikeThresh > 0),
    spk_data = -spk_data;
    spk_thresh = -spk_thresh;
    art_thresh = -art_thresh;
end

%Threshold our data
fprintf('\tThresholding.\n');
spk_thresh = (spk_data <= spk_thresh);
art_thresh = (spk_data <= art_thresh);

%Find where it crosses (and doesn't go back)
fprintf('\tFinding spikes...\n');
spk_ind = find(any(spk_thresh & (diff([zeros(1, num_chans); spk_thresh], [], 1) == 1), 2));
art_ind = find(any(art_thresh & (diff([zeros(1, num_chans); art_thresh], [], 1) == 1), 2));
clear('spk_thresh', 'art_thresh');
good_ind = true(length(spk_ind), 1);
for i = 1:(length(spk_ind)-1),
    if any((art_ind >= spk_ind(i)) & (art_ind <= spk_ind(i+1))),
        good_ind(i) = 0;
    end
end
if any(art_ind >= spk_ind(end)), good_ind(end) = 0; end
spk_ind = spk_ind(good_ind);
fprintf('\t\tfound %d putative spikes (%d removed as artifacts)...\n', length(spk_ind), sum(~good_ind));
if isempty(spk_ind), fprintf('Nothing to write.\n'); return; end

%Remove 'duplicates' where they are within a spike window of each other
bad_spk_ind = false(length(spk_ind), 1); cur_good = 1; min_dist = (opts.SpikeExtract_WFRange(2)/1000*opts.SampleFrequency);
for i = 2:length(spk_ind),
    if (spk_ind(i) - spk_ind(cur_good)) <= min_dist,
        bad_spk_ind(i) = 1;
    else
        cur_good = i;
    end
end
fprintf('\t\tremoved %d spikes within a window...\n', sum(bad_spk_ind));
fprintf('\t\tleaving %d spikes (%4.2f Hz).\n', sum(~bad_spk_ind), sum(~bad_spk_ind)./(max_spk_ind./opts.SampleFrequency));
clear('spk_thresh', 'art_thresh');

Nspk = length(spk_ind);

%Open output file, write header information
fid = fopen(save_fn, 'w');
fwrite(fid, 'EXTRACTEDSPIKEWF', 'char*1');
fwrite(fid, Nspk, 'int64');
fwrite(fid, opts.SampleFrequency, 'double');
fwrite(fid, opts.FilterSettleTime, 'double');
fwrite(fid, opts.SamplePoints, 'double');
fwrite(fid, opts.SpikeInterpFreq, 'double');
fwrite(fid, opts.SpikeExtract_ThreshType, 'uint8');
fwrite(fid, opts.SpikeExtract_SlidingWindow, 'uint8');
fwrite(fid, opts.SpikeExtract_SlidingWindowTimeWidth, 'double');
fwrite(fid, opts.SpikeExtract_SlidingWindowTimeStep, 'double');
fwrite(fid, opts.SpikeExtract_SpikeThresh, 'double');
fwrite(fid, opts.SpikeExtract_ArtifactThresh, 'double');
fwrite(fid, opts.SpikeExtract_WFRange, 'double');
fwrite(fid, opts.SpikeAlign, 'uint8');
fwrite(fid, opts.SpikeAlignRange, 'double');
fwrite(fid, Nwf, 'int32');

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

%Extract spike waveforms
spk_wf = zeros(length(wf_range_ind), num_chans, Nspk);
spline_wf = zeros(length(wf_interp_ind), num_chans, Nspk);
for cur_spk = 1:Nspk,
    % Grab the current waveform (plus more for alignment purposes if needed)
    cur_wf_ind = wf_range_ind + spk_ind(cur_spk);
    valid_cur_wf_ind = ((cur_wf_ind > 0) & (cur_wf_ind <= max_spk_ind));
    spk_wf(valid_cur_wf_ind, :, cur_spk) = spk_data(cur_wf_ind(valid_cur_wf_ind), :);
end
%Do spline fit on our local pieces
for cur_chan = 1:num_chans,
    spline_wf(:, cur_chan, :) = reshape(spline(wf_range_ind, squeeze(spk_wf(:, cur_chan, :))', wf_interp_ind)', [length(wf_interp_ind) 1 Nspk]);
end

if (opts.SpikeAlign > 0),
    %Align spikes on a limited range (if specified)
    fprintf('\tRe-aligning to peak on minimum channel...');
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
        [~, mi] = min(min(spline_wf(align_range_ind, :, :), [], 2), [], 1);
        mi = squeeze(mi);
    elseif opts.SpikeAlign == 2,
        %Find first local minima
        fprintf('finding first local minima...');
        mi = NaN*ones(Nspk, 1); tmi = Inf*ones(num_chans, 1); mv = Inf*ones(num_chans, 1);
        for cur_spk = 1:Nspk,
            for cur_chan = 1:num_chans,
                cur_mi = find(diff(diff(spline_wf(align_range_ind, :, cur_spk), 1, 1) >= 0, 1, 1) == 1, 1, 'first') + 1;
                if isempty(cur_mi), cur_mi = length(align_range_ind); end
                tmi(cur_chan) = cur_mi;
                mv(cur_chan) = spline_wf(align_range_ind(cur_mi), cur_chan, cur_spk);
            end %channel loop
            [~,ttmi] = min(mv);
            mi(cur_spk) = tmi(ttmi);
        end %spike loop
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
    
    for cur_spk = 1:Nspk,
        spline_wf(1:length(orig_wf_interp_ind), :, cur_spk) = spline_wf(align_range_ind(mi(cur_spk)) + orig_wf_interp_ind/interp_rescale_ratio, :, cur_spk);
    end
    spline_wf = spline_wf(1:length(orig_wf_interp_ind), :, :);
end

%Rescale
spline_wf = reshape(spline_wf(rescale_ind, :, :), [(length(rescale_ind)*num_chans) Nspk]);
if (opts.SpikeExtract_SpikeThresh > 0),
    spline_wf = -spline_wf; %Invert back for finding peaks
end

fprintf('\tWriting to file...');
%For each spike index grab its time, write to file
fwrite(fid, (spk_ind - 1)/opts.SampleFrequency, 'double');
%Write to disk (with zeros for cluster marker)
fwrite(fid, cat(1, spline_wf, zeros(1, Nspk)), 'double');

%Close file
fclose(fid);
fprintf('done.\n');