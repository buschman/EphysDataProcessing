function ReReferenceSignals(file_list, ref_list, varargin)

% Re-references neural signals.
%
% Inputs:
%   file_list       List of file to process.
%   ref_list        List of which files to use as reference (averaged if more than 1)
%
% Optional:
%   FileAppendStr       String to append to re-referenced files (default is 'ReRef')
%   RescaleRef          Flag for rescaling reference before subtracting (default = 1)
%   ReRefVar            Variable to re-reference (either a string or cell array
%                       of strings). Default is 'spk_data' and 'lfp_data' if they exist.
%   FilterSettleTime    Amount of time to cut from beginning/end due to the
%                       filter settling.  Default is 10s @ 25 kHz (i.e. 250,000 points)
%
% 12/11/12, TJB, first version

%Options
fopts.Overwrite = 0;
fopts.FileAppendStr = 'ReRef';
fopts.RescaleRef = 1;
fopts.ReRefVar = {'spk_data', 'lfp_data'};
fopts.FilterSettleTime = 10*25000; %how many indexes to ignore at start/end of file (default is 10s @ 25kHz)

if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs.'); end
for i = 1:2:length(varargin),
    if ~isfield(fopts, varargin{i}),
        error(sprintf('Key %s not valid.', varargin{i}));
    end
    fopts.(varargin{i}) = varargin{i+1};
end

%Expand file_list if just a string
if isstr(file_list),
    orig_file_list = file_list;
    filedir = fileparts(orig_file_list);
    file_list = dir(file_list);
    file_list = struct2cell(file_list);
    file_list = file_list(1, :);
    good_ind = [];
    for i = 1:length(file_list),
        if isempty(strfind(file_list{i}, 'ReRef')), %don't use files already re-referenced
            good_ind = cat(1, good_ind, i);
        end
    end
    file_list = file_list(good_ind);
    file_list = strcat(filedir, '\', file_list(1, :));
    fprintf('Processing %d files that match ''%s''...\n', length(file_list), orig_file_list);
    clear('filedir', 'orig_file_list');
end

%Check inputs
if length(file_list) < 2, error('Must specify at least 2 files.'); end
if isempty(ref_list), 
    fprintf('Using all channels as a reference.\n');
    ref_list = [1:length(file_list)];
end
if ~all(ismember(ref_list, [1:length(file_list)])), error('Reference list must index into the current file list.'); end
if isempty(fopts.ReRefVar), return; end

% Make sure all files exist, check what variables are in each file
var_exists = false(length(file_list), length(fopts.ReRefVar));
for cur_file = 1:length(file_list), 
    if ~exist(file_list{cur_file}, 'file'), error(sprintf('Can''t find %s.', file_list{cur_file})); end
    
    %Which variables exist
    vars = struct2cell(whos('-file', file_list{cur_file}));
    vars = vars(1, :);
    for cur_var = 1:length(fopts.ReRefVar),
        var_exists(cur_file, cur_var) = any(strcmpi(fopts.ReRefVar(cur_var), vars));
    end    
end
% Make sure re-reference variables are in all of the specified files
for cur_var = 1:length(fopts.ReRefVar),
    if ~all(var_exists(:, cur_var)),
        fprintf('Couldn''t find variable %s in all of the reference files.\n', fopts.ReRefVar{cur_var});
    end
end
fopts.ReRefVar = fopts.ReRefVar(all(var_exists, 1));
if isempty(fopts.ReRefVar), 
    fprintf('No variables to process.\n');
end
       
%% Process each variable in turn
for cur_var = 1:length(fopts.ReRefVar),
    
    fprintf('Re-referencing variable %s...\n', fopts.ReRefVar{cur_var});
    
    %Create grand-averaged version
    fprintf('Calculating overall average variable: [%4.0f/%4.0f]', 1, length(ref_list));
    temp = load(file_list{ref_list(1)}, fopts.ReRefVar{cur_var});
    ovr_avg = double(temp.(fopts.ReRefVar{cur_var}));
    for cur_ref = 2:length(ref_list),
        fprintf('\b\b\b\b\b\b\b\b\b\b\b[%4.0f/%4.0f]', cur_ref, length(ref_list));
        temp = load(file_list{ref_list(cur_ref)}, fopts.ReRefVar{cur_var});
        ovr_avg = ovr_avg + double(temp.(fopts.ReRefVar{cur_var}));
    end
    fprintf('\n');
    ovr_avg = ovr_avg./length(ref_list);
    
    %What indexes are 'good'?
    good_ind = [fopts.FilterSettleTime:(length(ovr_avg) - fopts.FilterSettleTime)];
        
    %Remove from each file in turn
    for cur_file = 1:length(file_list),
        
        fprintf('Re-referencing file %d: %s...\n', cur_file, file_list{cur_file});
        %Load the current value
        temp = load(file_list{cur_file}, fopts.ReRefVar{cur_var}, 'opts');
        opts = temp.opts;
        temp = double(temp.(fopts.ReRefVar{cur_var}));
        
        %Re-reference
        if ~fopts.RescaleRef,
            fprintf('\tSubtracting raw reference...\n');
            temp(good_ind) = temp(good_ind) - ovr_avg(good_ind);
        else
            fprintf('\tSubtracting scaled reference...');
            %A = ovr_avg(good_ind)\temp(good_ind);
            A = pinv(ovr_avg(good_ind))*temp(good_ind);
            fprintf('scaled by %5.3e.\n', A);
            temp(good_ind) = temp(good_ind) - A*ovr_avg(good_ind);
        end
        temp = single(temp);
                    
        %Save to file
        [pathstr, name, ext] = fileparts(file_list{cur_file});
        save_fn = strcat(pathstr, name, fopts.FileAppendStr, ext);
        eval(sprintf('%s = temp;', fopts.ReRefVar{cur_var}));
        clear('temp');
        if exist(save_fn, 'file'),
            save(save_fn, fopts.ReRefVar{cur_var}, 'opts', '-append', '-v7.3');
        else
            save(save_fn, fopts.ReRefVar{cur_var}, 'opts', '-v7.3');
        end
        clear('opts');
    end %file loop
    
end