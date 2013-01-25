function AddSimpleClustToWF(varargin),

opts.UpdateWaitbarStep = 10000; %every this # of ts; 0 turns it off

%GUI to select files?
if isempty(varargin),
    %Use dialog
    [FileName,PathName] = uigetfile('*.wf', 'Select waveform files to update...', 'MultiSelect', 'on');
    wf_file = strcat(PathName, FileName);
else
    %Change to cell if a single file passed
    if ischar(varargin{1}),
        wf_file = varargin(1);
    else
        wf_file = varargin{1};
    end
end

if opts.UpdateWaitbarStep, wh = waitbar(0, 'Adding cluster information to WF files...'); end

%Loop through files
for cur_file = 1:length(wf_file), 
    if opts.UpdateWaitbarStep, waitbar((cur_file - 1)./length(wf_file), wh); end
    
    fn = wf_file{cur_file}; 
    %Check if file exists
    if ~exist(fn, 'file'),
        fprintf('File %s doesn''t exist. Skipping.\n', wf_file{cur_file});
    end
    
    %Simple clust file
    [pathstr, name, ext] = fileparts(fn);
    chstr = regexpi(name, '(\w+)_Ch(\d)', 'tokens');
    clust_fn = sprintf('%s\\%s_Ch%s_clustered.mat', pathstr, chstr{1}{1}, chstr{1}{2});
    %Check if file exists
    if ~exist(clust_fn, 'file'),
        fprintf('File %s doesn''t have a simple_clust match. Skipping.\n', wf_file{cur_file});
    end
       
    
    %Check to make sure the file is valid
    fid = fopen(fn, 'r+');
    c = char(fread(fid, length('EXTRACTEDSPIKEWF'), 'char*1'));
    if ~all(c(:)' == char('EXTRACTEDSPIKEWF')), error('Passed filename does not appear to be valid waveform file.'); end
    
    Nspk = fread(fid, 1, 'int64');
    
    %Load other variables (could skip most of these in the future?)
    opts.SampleFrequency = fread(fid, 1, 'double'); %64 bits, 8 bytes
    opts.FilterSettleTime = fread(fid, 1, 'double'); %64 bits, 8 bytes
    opts.SamplePoints = fread(fid, 1, 'double'); %64 bits, 8 bytes
    opts.SpikeInterpFreq = fread(fid, 1, 'double'); %64 bits, 8 bytes
    opts.SpikeExtract_ThreshType = fread(fid, 1, 'uint8'); %8 bits
    opts.SpikeExtract_SlidingWindow = fread(fid, 1, 'uint8'); %8 bits
    opts.SpikeExtract_SlidingWindowTimeWidth = fread(fid, 1, 'double'); %64 bits
    opts.SpikeExtract_SlidingWindowTimeStep = fread(fid, 1, 'double'); %64 bits
    opts.SpikeExtract_SpikeThresh = fread(fid, 1, 'double'); %64 bits
    opts.SpikeExtract_ArtifactThresh = fread(fid, 1, 'double'); %64 bits
    opts.SpikeExtract_WFRange = fread(fid, 2, 'double'); %64 bits
    opts.SpikeAlign = fread(fid, 1, 'uint8');
    opts.SpikeAlignRange = fread(fid, 2, 'double'); %64 bits
    Nwf = fread(fid, 1, 'int32'); %32 bits
    opts.NumPointsInWF = Nwf;

    %Read in timestamps
    ts = fread(fid, Nspk, 'double');

    %Load the cluster information
    spikes = load(clust_fn, 'spikes');
    spikes = spikes.spikes;
    
    %Match the timestamps
    [ts_match, ts_clust_ind] = ismember(ts, spikes.ts);
    
    %Convert each category to negative (noise), positive (for mua and single units)
    [uniq_clust, ~, uniq_clust_ind] = unique(spikes.cluster_is(ts_clust_ind));
    num_noise = -1; num_mua = 0; num_unit = 100;
    clust_map = zeros(length(uniq_clust), 1);
    for i = 1:length(uniq_clust),
        if ~isempty(regexpi(spikes.labelcategories{spikes.clusterlabels(uniq_clust(i))}, '(unit)')),
            clust_map(i) = num_unit;
            num_unit = num_unit + 1;
        elseif ~isempty(regexpi(spikes.labelcategories{spikes.clusterlabels(uniq_clust(i))}, '(mua)')) || isempty(spikes.labelcategories{spikes.clusterlabels(uniq_clust(i))}),
            clust_map(i) = num_mua;
            num_mua = num_mua + 1;
        elseif ~isempty(regexpi(spikes.labelcategories{spikes.clusterlabels(uniq_clust(i))}, '(noise|artefact)')),
            clust_map(i) = num_noise;
            num_noise = num_noise - 1;
        end
    end
    
    %Move through waveforms, updating as we go
    for i = 1:length(ts),
        fseek(fid, Nwf*8, 'cof');
        if ~ts_match(i),
            fseek(fid, 8, 'cof');
        else
            fwrite(fid, clust_map(uniq_clust_ind(ts_clust_ind(i))), 'double');
        end
        if opts.UpdateWaitbarStep && (mod(i, opts.UpdateWaitbarStep) == 0), waitbar((cur_file - 1)./length(wf_file) + i./length(ts)./length(wf_file), wh); end
    end
    fclose(fid);
    clear spikes;
end
close(wh);