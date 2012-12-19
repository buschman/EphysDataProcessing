function ExtractAllSpikeWF(fn, varargin)


if ischar(fn),
    fn = {fn};
end

for cur_fn = 1:length(fn),
       
    filedir = fileparts(fn{cur_fn});
    file_list = dir(fn{cur_fn});
    file_list = struct2cell(file_list);
    files = strcat(filedir, '\', file_list(1, :));
    fprintf('Processing %d files that match ''%s''...\n', length(files), fn{cur_fn});
    
    for i = 1:length(files),
        [pathstr, name, ext] = fileparts(files{i});
        fprintf('Analyzing %s...\n', name);
        load(files{i}, 'spk_data');
        if ~isempty(pathstr), save_fn = sprintf('%s\\%s.wf', pathstr, name);
        else save_fn = sprintf('%s.wf', name); end
        ExtractTrodeSpikeWF(spk_data, save_fn, varargin);
        clear('spk_data');
    end
end