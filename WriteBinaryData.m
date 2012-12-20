function WriteBinaryData(mat_file, out_file, varargin),

%Options
opts.Overwrite = [];
opts.VarToWrite = 'spk_data';
opts.SampleFrequency = 25000;
opts.MaxScale = 1; %1 for scaling as best as possible, 0 for using scale factor
opts.ScaleFactor = 1000;
opts.BitDepth = 16;
opts.FilterSettleTime = 5*60; %in seconds, this amount of data will be zero'd
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
out_fid = fopen(out_file, 'w');

%Read data
if ~exist(mat_file, 'file'), error('File %s doesn''t exist.', mat_file); end
data = load(mat_file, opts.VarToWrite);
data = data.(opts.VarToWrite);
data(1:(opts.SampleFrequency*opts.FilterSettleTime), :) = 0;
data((end - opts.SampleFrequency*opts.FilterSettleTime + 1):end, :) = 0;

%Scale data to fit bit depth maximally?
if opts.MaxScale,
    opts.ScaleFactor = max(abs(data), [], 1);
    opts.ScaleFactor = max(opts.ScaleFactor); %keep everything equal
end

%Write data
switch opts.BitDepth,
    case 8,
        fwrite(out_fid, int8(data'./opts.ScaleFactor.*(2^(opts.BitDepth-1)-1)), 'int8');
    case 16,
        fwrite(out_fid, int16(data'./opts.ScaleFactor.*(2^(opts.BitDepth-1)-1)), 'int16');
    case 32,
        fwrite(out_fid, int32(data'./opts.ScaleFactor.*(2^(opts.BitDepth-1)-1)), 'int32');
    case 64,
        fwrite(out_fid, int64(data'./opts.ScaleFactor.*(2^(opts.BitDepth-1)-1)), 'int64');
    otherwise,
        error('Unsupported bit-depth.  Should be a 8, 16, 32, or 64 bits.  Or add code to support intermediates.');
end

fclose(out_fid);
