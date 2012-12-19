function [t, varargout] = LoadSpikeWF(fn, varargin)

% Loading engine requirements
% A loading engine must be callable as a Matlab function.  This means it must be either a .m function-file or a
% compiled mex function-file.  It must take as input one or three inputs and provide one or two outputs.
% MClust-3.0 should work with any loading engine that supplies this functionality.
%
% INPUTS
% •	fn = file name string
% •	records_to_get = a range of values
% •	record_units = a flag taking one of 5 cases (1,2,3,4 or 5)
% 1.	implies that records_to_get is a timestamp list.
% 2.	implies that records_to_get  is a record number list
% 3.	implies that records_to_get  is range of timestamps (a vector with 2 elements: a start and an end timestamp)
% 4.	implies that records_to_get  is a range of records (a vector with 2 elements: a start and an end record number)
% 5.	asks to return the count of spikes (records_to_get should be [] in this case)
% In addition, if only fn is passed in then the entire file should be read.
% OUTPUT
% •	[t, wv]
% •	t = n x 1: timestamps of each spike in file
% •	wv = n x 4 x 32 waveforms
% EXAMPLES
% •	[t,wv] = myLoadingEngine(‘myfile.dat’, 1:10, 2) should return the time and waveforms for the first 10 spikes in the file.
% •	t = myLoadingEngine(‘myfile.dat’) should return all the timestamps from the file.
% •	n = myLoadingEngine(‘myfile.dat’, [], 5) should return the number of spikes in the file.

if ~isempty(varargin),
    records_to_get = varargin{1};
    record_units = varargin{2};
else
    record_units = [];
    records_to_get = [];
end

%Check to make sure the file is valid
fid = fopen(fn, 'r');
c = char(fread(fid, length('EXTRACTEDSPIKEWF'), 'char*1'));
if ~all(c(:)' == char('EXTRACTEDSPIKEWF')), error('Passed filename does not appear to be valid waveform file.'); end

Nspk = fread(fid, 1, 'int64');
%Just asking for the # of spikes?
if record_units == 5,
    t = Nspk;
    fclose(fid);
    return;
end

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

%Add an option to get information about extraction
if record_units == 6,
    t = opts;
    fclose(fid);
    return;
end

if isempty(record_units),
    %Return timestamps
    t = fread(fid, Nspk, 'double');
    fclose(fid);
    return;
end
if (record_units == 3) || (record_units == 4),
    if (record_units == 3),
        %Read all of the time points in
        t = fread(fid, Nspk, 'double');
        %Convert timestamp range to record range
        records_to_get(1) = find(t >= records_to_get(1), 1, 'first');
        records_to_get(2) = find(t <= records_to_get(2), 1, 'last');
        t = t(records_to_get(1):records_to_get(2));
        num_records = diff(records_to_get) + 1;
    elseif (record_units == 4),
        num_records = diff(records_to_get) + 1;
        %Seek to start of range in the time array
        fseek(fid, (records_to_get(1) - 1)*8, 'cof');
        %Read data
        t = fread(fid, num_records, 'double');
        %Seek to start of wf data
        fseek(fid, (Nspk - records_to_get(2))*8, 'cof');
    end
    
    %Seek to start of range in wf data
    fseek(fid, (Nwf+1)*(records_to_get(1) - 1)*8, 'cof');
    %Read data
    wf = fread(fid, [(Nwf+1) num_records], 'double');
    wf = wf(1:Nwf, :);
elseif (record_units == 1) || (record_units == 2),
    %Read all of the time points in
    t = fread(fid, Nspk, 'double');
    
    if (record_units == 1),
        %Convert specific times to specific records
        [valid_ts, records_to_get] = ismember(records_to_get, t);
        records_to_get = records_to_get(valid_ts);
    end
    
    %We need to extract specific records -- probably faster to read all in rather than looping through each
    wf = fread(fid, [(Nwf+1) Nspk], 'double');
    wf = wf(1:Nwf, records_to_get);
    t = t(records_to_get);
else
    error('Unexpected record_units passed.');
end

%Do we need to pad with zeros?
if (size(wf, 1) < 128),
    wf = cat(1, wf, zeros(128-size(wf, 1), size(wf, 2)));
elseif (size(wf, 1) > 128),
    wf = wf(:, 1:128);
end
    
%Convert our Nwf x Nrecords to Nrecords x Nwf for now
wf = wf';

%Reshuffle column indexes (going 1-4 in turn)
reshape_ind = repmat([0:3]'*32, [1 32]) + repmat([1:32], [4 1]);
wf = wf(:, reshape_ind(:));

%Reshape
wf = reshape(wf, [size(wf, 1) 4 32]);
varargout{1} = wf;