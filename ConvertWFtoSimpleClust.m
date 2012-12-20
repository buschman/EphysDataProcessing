function ConvertWFtoSimpleClust(wf_file, out_file, varargin),

%Options
opts.Overwrite = [];
opts.SampleFrequency = 25000;
opts.SourceChannel = [];
%Check optional inputs
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs.'); end
for i = 1:2:length(varargin),
    if ~isfield(opts, varargin{i}),
        error(sprintf('Key %s not valid.', varargin{i}));
    end
    opts.(varargin{i}) = varargin{i+1};
end

%Check if output file exists
if exist(out_file, 'file'),
    if isempty(opts.Overwrite),
        dlgans = questdlg(sprintf('File %s exists. Overwrite?', out_file), 'Overwrite file?', 'Yes', 'No', 'Yes');
        if strcmpi(dlgans, 'yes'), opts.Overwrite = 1; else, opts.Overwrite = 0; end
    end
    if ~opts.Overwrite, fprintf('File %s already exists. Skipping.\n', out_file); end
end

%Load waveforms
n = LoadSpikeWF(wf_file, [], 5);
[t, wv] = LoadSpikeWF(wf_file, [1 n], 4);

%Create MUA variable
mua.ts = t; 
mua.Nspikes = n;
if ~isempty(opts.SourceChannel),
    mua.sourcechannel = opts.SourceChannel;
else
    chstr = regexpi(wf_file, '_Ch(\d)', 'tokens');
    if isempty(chstr),
        mua.sourcechannel = 0;
    else
        %We'll take the first channel as our channel
        mua.sourcechannel = str2num(chstr{1}{1});
    end
end

%Copy waveforms over
mua.waveforms = zeros(size(wv, 3), size(wv, 2), size(wv, 1));
for i = 1:size(wv, 1),
    mua.waveforms(:, :, i) = squeeze(wv(i, :, :))';
end

%Save to file
save(out_file, 'mua');