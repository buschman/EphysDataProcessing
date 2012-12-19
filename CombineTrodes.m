function CombineTrodes(base_fn, ch_list, varargin),

% CombineTrodes(base_fn, ch_list)
%
% Combine multiple electrodes for procesing as trodal.
%
% Inputs:
%   base_fn         Base file name.  Use wildcard * to indicate where
%                   Ch + channel number should go.  So, for Test_Ch5.mat
%                   you should pass 'Test_*.mat'.
%   ch_list         Cell array of channels to combine.  Each element in
%                   the array should be a list of channels to combine into
%                   a single file.
%
% Optional Inputs:
%   var_name        String of variable to load.  Default is 'spk_data'.
%   overwrite       Whether to forcefully overwrite file.
%
% Outputs:          None.
%
% 12/17/12  Written TJB

if (length(varargin) >= 1) && ~isempty(varargin{1}),
    var_name = varargin{1};
else
    var_name = 'spk_data';
end

if (length(varargin) >= 2) && ~isempty(varargin{2}),
    overwrite = varargin{2};
else
    overwrite = 0;
end

%Loop through new trodes
for cur_trode = 1:length(ch_list),
    cur_ch_list = ch_list{cur_trode};
    if isempty(cur_ch_list), continue; end
    
    %Initialize our Matlab file for writing
    channel_string = sprintf('Ch%d', cur_ch_list(:));
    save_fn = strrep(base_fn, '*', channel_string);
    if exist(save_fn, 'file') & ~overwrite,
        fprintf('File %s already exists. Skipping.\n', save_fn);
    end
    if exist(save_fn, 'file'), delete(save_fn); end
    %Load options from first file (all others will be ignored)
    opts = load(strrep(base_fn, '*', sprintf('Ch%d', cur_ch_list(1))), 'opts');
    opts = opts.opts;
    save(save_fn, 'ch_list', 'opts', '-v7.3');
    saveMatObj = matfile(save_fn, 'Writable', true);
    
    %Loop through channels
    for cur_chan = 1:length(cur_ch_list),
        %Open file for reading
        load_fn = sprintf('Ch%d', cur_ch_list(cur_chan));
        load_fn = strrep(base_fn, '*', load_fn);
        loadMatObj = matfile(load_fn, 'Writable', false);
        [Nt, ~] = size(loadMatObj, var_name);
        
        %Write to file
        saveMatObj.(var_name)(1:Nt, cur_chan) = loadMatObj.(var_name)(1:Nt, 1);
        
        clear('loadMatObj');
    end %channel loop
    clear('saveMatObj');
end %trode loop