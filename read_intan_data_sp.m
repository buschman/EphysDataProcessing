function [amps,data,aux,varargout] = read_intan_data(varargin)

% [amps,data,aux, (t)] = read_intan_data
%
% Opens file selection GUI to select and then read data from an Intan
% amplifier data file (*.int).
%
% t = time vector (in seconds)
% amps = vector listing active amplifier channels
% data = matrix of electrode-referred amplifier signals (in microvolts)
% aux = matrix of six auxiliary TTL input signals
%
% Example usage:
%  >> [t,amps,data,aux] = read_intan_data;
%  >> plot(t,data(:,1));
%
% Version 1.1, June 26, 2010
% (c) 2010, Intan Technologies, LLC
% For more information, see http://www.intantech.com
% For updates and latest version, see http://www.intantech.com/software.html
%
% 06-22-10 Added GUI file selection and optimized: Craig Patten, Plexon, Inc.
% 09-06-12 Add ability to specify filename, select channels, and some options, TJB 9/6/12

% Filename specified?
if isempty(varargin) || isempty(varargin{1}),
    % use MATLAB predefined gui uigetfile to select the file(s) to analyze
    [file, path, filterindex] = uigetfile('*.int','Select a .int file','MultiSelect', 'off');
    filename = [path,file];
else
    filename = varargin{1};
end

%Parse optional variables
varargin = varargin(2:end);
opts.ReadChannels = [];
opts.TimeChunkSize = 5*10^6; %in time samples
opts.MaxFastTrackSize = 2*10^9; %in bytes
opts.UpdateChunkNumber = 5; %in chunks
opts.Verbose = 0;
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs.'); end
for i = 1:(length(varargin)/2),
    opts.(varargin{2*i-1}) = varargin{2*i};
end

%Open file
fid = fopen(filename, 'r');

% Read first three header bytes encoding file version
for i=1:3
    header(i) = fread(fid, 1, 'uint8');
end

if (header(1) ~= 128)
    error('Improper data file format.');
end

if (header(2) ~= 1 || header(3) ~= 1)
    warning('Data file version may not be compatible with this m-file.');
end

% Now see which amplifier channels are saved in this file.
for i=1:64
    amp_on(i) = fread(fid, 1, 'uint8');
end
num_amps = sum(amp_on);

% Create a list of amplifier channels in this file.
amps = find(amp_on == 1);
if isempty(opts.ReadChannels),
    opts.ReadChannels = amps;
end

% Now search for the end of the file to find out the length of the data.
% t_count = 0;
% while (~feof(fid))
%    fread(fid, 1+4*num_amps, 'uint8');
%    t_count = t_count + 1;
% end
% t_count = t_count - 1;
% t_max = t_count/25000;

%-----------------------------------
% replace above code with a more efficient method CDP 06-24-10
s = dir(filename);
filesize = s.bytes;
t_count = (filesize - 67)/(num_amps*4 + 1);
if rem(t_count, 1) ~= 0,
    warning('File appears to have been truncated.  Will continue to process with available timepoints.');
    t_count = floor(t_count);
end
t_max = t_count/25000;
%-----------------------------------

% print channel (singular) when there is only one channel! CDP 06-24-10
if opts.Verbose,
if num_amps == 1;
    fprintf(1, '\nData file contains %0.2f seconds of data (%d points) from %d amplifier channel.\n', t_max, t_count, num_amps);
    fprintf(1, 'Channel: ');
else
    fprintf(1, '\nData file contains %0.2f seconds of data (%d points) from %d amplifier channels.\n', t_max, t_count, num_amps);
    fprintf(1, 'Channels: ');
end
for i=1:num_amps
    fprintf(1, '%d ', amps(i));
end
fprintf(1, '\n');
end

%Write time vector if requested
if nargout >= 4,
    varargout{1} = [0:1:(t_count-1)]'/25000;
end    

%--------------------------------------
% Replace code code below with much faster code CDP 06-24-10
% Go back to the beginning of the file...
frewind(fid);

% ...skip the header this time...
fread(fid, 3+64, 'uint8');

% allocate space to read the entire file
slow_road = 1;
if (filesize-67) > opts.MaxFastTrackSize,
    if opts.Verbose, fprintf('Can''t allocate an array of that size. '); end
    slow_road = 1;
end

if ~slow_road,   
    % read the entire file
    %data2 = fread(fid,(filesize-67),'uint8=>uint8');
    data2 = fread(fid, (num_amps*4+1)*t_count, 'uint8=>uint8');
    
    % extract the digital data
    aux_data = data2((num_amps*4)+1:num_amps*4+1:filesize-67);
    
    % extract individual bits
    aux = [bitget(aux_data,6),bitget(aux_data,5),bitget(aux_data,4),bitget(aux_data,3),bitget(aux_data,2),bitget(aux_data,1)];
    clear aux_data;
    
    % delete the digital data
    %data2((num_amps*4)+1:num_amps*4+1:filesize-67) = [];
        data2((num_amps*4)+1:num_amps*4+1:end) = [];

    
    % convert the remaining data from bytes to single
    data2 = typecast(data2,'single');
    
    data = zeros(t_count,num_amps);
    % de-mux the channels
    for ind = 1:num_amps
        data(:,ind) = data2(ind:num_amps:end);
    end
    %Extract our channels to read
    data = data(:, ismember(amps, opts.ReadChannels));

elseif ~isempty(opts.ReadChannels),
    % Go back to the beginning of the file...
    frewind(fid);
    
    % ...skip the header this time...
    fread(fid, 3+64, 'uint8');
    data_start_pos = ftell(fid);
    
    %Move to beginning of auxilary data
    fseek(fid, num_amps*4, 'cof');
    % read the auxilary data
    aux_data = fread(fid, t_count, 'uint8', 4*num_amps);    
    % extract individual bits
    aux = [bitget(aux_data,6),bitget(aux_data,5),bitget(aux_data,4),bitget(aux_data,3),bitget(aux_data,2),bitget(aux_data,1)];
    clear aux_data;
    
    %Pre-allocate data channel
    data = single(NaN*ones(t_count, length(opts.ReadChannels)));
    %Loop through channels
    for chan_ind = 1:length(opts.ReadChannels),
        %What is current channel
        cur_chan = opts.ReadChannels(chan_ind);
        amp_ind = find(amps == cur_chan, 1, 'first');
        
        %Move to beginning of data section
        fseek(fid, data_start_pos, 'bof');
        
        %Move forward channel-1 positions
        fseek(fid, (amp_ind-1)*4, 'cof');
        
        %Read data, skipping intervening channels (and auxilary data)
        data(:, chan_ind) = fread(fid, t_count, 'single', 4*(num_amps-1)+1);
        
    end %channel loop
else
    % Go back to the beginning of the file...
    frewind(fid);
    
    % ...skip the header this time...
    fread(fid, 3+64, 'uint8');
    data_start_pos = ftell(fid);
    
    % ...and read all the data.
    if opts.Verbose, fprintf('Reading data in chunks. (This may take a while.)\n'); end
    
    %What channels to read
    [~, amp_ind] = ismember(opts.ReadChannels, amps);
    
    %Create matrices
    wh = waitbar(0, 'Creating matrices...'); drawnow;
    data = zeros(t_count, length(opts.ReadChannels), 'single');
    aux = zeros(t_count, 6, 'uint8');
    
    %Read data
    waitbar(0, wh, 'Reading data in chunks...'); start_time = now; drawnow;
    chunk_count = 0;
    for i=1:opts.TimeChunkSize:t_count,
        chunk_count = chunk_count + 1;
        
        time_ind = i + [0:(opts.TimeChunkSize-1)];
        time_ind  = time_ind((time_ind <= t_count));
        time_data = fread(fid, (num_amps*4+1)*length(time_ind), 'uint8=>uint8');

        %Extract out aux inputs
        aux_data = time_data((num_amps*4)+1:num_amps*4+1:end);
        % extract individual bits
        aux(time_ind, :) = [bitget(aux_data,6),bitget(aux_data,5),bitget(aux_data,4),bitget(aux_data,3),bitget(aux_data,2),bitget(aux_data,1)];
        clear aux_data;
        
        % delete the aux data
        time_data((num_amps*4)+1:num_amps*4+1:end) = [];
        
        % convert the remaining data from bytes to single
        time_data = typecast(time_data,'single');
        time_data = reshape(time_data, [num_amps length(time_ind)]);
        
        %Save to our variable
        data(time_ind, :) = time_data(ismember(amps, opts.ReadChannels), :)';
        %         for j=1:num_amps,
        %             %temp_val = double(fread(fid, 1, 'float32'));
        %             temp_val = fread(fid, 1, 'single'); %this should be equivalent to float32?
        %             data(i, ismember(opts.ReadChannels, j)) = temp_val;
        %         end %channel loop
        
        %         aux_byte = fread(fid, 1, 'uint8');
        %
        %         % Decode auxiliary TTL inputs
        %         aux(i, :) = bitget(aux_byte, [6:-1:1]);
        %         if aux_byte >= 32
        %             aux(i,6) = 1;
        %             aux_byte = aux_byte - 32;
        %         end
        %         if aux_byte >= 16
        %             aux(i,5) = 1;
        %             aux_byte = aux_byte - 16;
        %         end
        %         if aux_byte >= 8
        %             aux(i,4) = 1;
        %             aux_byte = aux_byte - 8;
        %         end
        %         if aux_byte >= 4
        %             aux(i,3) = 1;
        %             aux_byte = aux_byte - 4;
        %         end
        %         if aux_byte >= 2
        %             aux(i,2) = 1;
        %             aux_byte = aux_byte - 2;
        %         end
        %         if aux_byte >= 1
        %             aux(i,1) = 1;
        %             aux_byte = aux_byte - 1;
        %         end
        if mod(chunk_count, opts.UpdateChunkNumber) == 1,
            waitbar(i./t_count, wh, sprintf('Reading data...[%s to complete]', datestr((now-start_time)/i*(t_count-i), 'HH:MM:SS'))); drawnow;
        end
    end %time loop
    close(wh); drawnow;
    
    %Correct channel listing
    amps = opts.ReadChannels;
end
% Close file, and we're done.
fclose(fid);
