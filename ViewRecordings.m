function varargout = ViewRecordings(varargin)
% VIEWRECORDINGS MATLAB code for ViewRecordings.fig
%      VIEWRECORDINGS, by itself, creates a new VIEWRECORDINGS or raises the existing
%      singleton*.
%
%      H = VIEWRECORDINGS returns the handle to a new VIEWRECORDINGS or the handle to
%      the existing singleton*.
%
%      VIEWRECORDINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEWRECORDINGS.M with the given input arguments.
%
%      VIEWRECORDINGS('Property','Value',...) creates a new VIEWRECORDINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ViewRecordings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ViewRecordings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ViewRecordings

% Last Modified by GUIDE v2.5 10-Dec-2012 15:09:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @ViewRecordings_OpeningFcn, ...
    'gui_OutputFcn',  @ViewRecordings_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ViewRecordings is made visible.
function ViewRecordings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ViewRecordings (see VARARGIN)

% Choose default command line output for ViewRecordings
handles.output = hObject;

%Initialize GUI
set([handles.ContinuousAxes handles.HistAxes handles.RawSpikeAxes handles.AvgSpikeAxes], 'YTick', [], 'XTick', []);
set(handles.TimeSlider, 'Enable', 'off');

handles.DataVarNameList = {'spk_data', 'chan_data'};
handles.Filename = '';

handles.ContTimeViewWidth = 10; %in seconds
handles.ContAxes = [];
handles.ZoomStep = [1 0.2];

%How many waveforms to draw as image compared to lines
handles.WaveformImageCutoff = 1;
handles.NumImageSteps = 40;
handles.ImageTimeInterpScale = 5;
handles.ImageContrastGamma = 1;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ViewRecordings wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ViewRecordings_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function FileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to FileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FileEdit as text
%        str2double(get(hObject,'String')) returns contents of FileEdit as a double

%Turn off everything (except a few items)
set(cat(1, findobj('Style', 'radiobutton'), findobj('Style', 'edit'), findobj('Style', 'slider'), ...
    findobj('Style', 'pushbutton')), 'Enable', 'off');
axes(handles.ContinuousAxes); cla; axis off;
axes(handles.HistAxes); cla; axis off;
axes(handles.RawSpikeAxes); cla; axis off;
axes(handles.AvgSpikeAxes); cla; axis off; legend off;
set([hObject handles.BrowseButton], 'Enable', 'on');

handles.Filename = get(hObject,'String');
if ~exist(handles.Filename, 'file'),
    error('Couldn''t find specified filename.');
    handles.Filename = '';
    handles.WFFilename = '';
end

%Set waveform filename
[pathstr, name, ext] = fileparts(handles.Filename);
handles.WFFilename = sprintf('%s\\%s.wf', pathstr, name);
if ~exist(handles.WFFilename, 'file'), handles.WFFilename = ''; end

%Create matlab file object
handles.DataMatObj = matfile(handles.Filename);
opts = load(handles.Filename, 'opts');
handles.SampleFrequency = opts.SampleFrequency;

%Update variable selector
file_vars = struct2cell(whos('-file', handles.Filename));
is_var = ismember(file_vars(1, :), handles.DataVarNameList);
set(handles.VarPopup, 'String', file_vars(1, is_var), 'Enable', 'on');
set(handles.VariableText, 'Enable', 'on');

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function FileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function TimeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to TimeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

handles = PlotContinuousData(handles);
handles = PlotRawWFData(handles);



% --- Executes during object creation, after setting all properties.
function TimeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TimeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in BrowseButton.
function BrowseButton_Callback(hObject, eventdata, handles)
% hObject    handle to BrowseButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile({'*.mat',  'MAT-files (*.mat)'; ...
    '*.*',  'All Files (*.*)'}, ...
    'Load data file...');

%Turn off everything (except a few items)
set(cat(1, findobj('Style', 'radiobutton'), findobj('Style', 'edit'), findobj('Style', 'slider'), ...
    findobj('Style', 'pushbutton')), 'Enable', 'off');
axes(handles.ContinuousAxes); cla; axis off;
axes(handles.HistAxes); cla; axis off;
axes(handles.RawSpikeAxes); cla; axis off; 
axes(handles.AvgSpikeAxes); cla; axis off; legend off;
set(hObject, 'Enable', 'on');

%Set filename
handles.Filename = sprintf('%s%s', pathname, filename);
set(handles.FileEdit, 'Enable', 'on', 'String', handles.Filename);

%Set waveform filename
[pathstr, name, ext] = fileparts(handles.Filename);
handles.WFFilename = sprintf('%s\\%s.wf', pathstr, name); handles.WFFiletype = 0;
if ~exist(handles.WFFilename, 'file'), handles.WFFilename = ''; handles.WFFiletype = -1; end

%Create matlab file object
handles.DataMatObj = matfile(handles.Filename);
load(handles.Filename, 'opts');
handles.SampleFrequency = opts.SampleFrequency;

%Update variable selector
file_vars = struct2cell(whos('-file', handles.Filename));
is_var = ismember(file_vars(1, :), handles.DataVarNameList);
set(handles.VarPopup, 'String', file_vars(1, is_var), 'Enable', 'on');
set(handles.VariableText, 'Enable', 'on');

guidata(hObject, handles);

% --- Executes on selection change in VarPopup.
function VarPopup_Callback(hObject, eventdata, handles)
% hObject    handle to VarPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns VarPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from VarPopup

contents = cellstr(get(hObject,'String'));
cur_var = contents{get(hObject,'Value')};

LoadData(handles, cur_var);

% --- Executes during object creation, after setting all properties.
function VarPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to VarPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in YZoomInButton.
function YZoomInButton_Callback(hObject, eventdata, handles)
% hObject    handle to YZoomInButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.ContAxes),
    handles.ContAxes(3:4) = handles.ContAxes(3:4)*2^(-handles.ZoomStep(2));
end
handles = PlotContinuousData(handles);
handles = PlotRawWFData(handles);
guidata(hObject, handles);


% --- Executes on button press in YZoomOutButton.
function YZoomOutButton_Callback(hObject, eventdata, handles)
% hObject    handle to YZoomOutButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.ContAxes),
    handles.ContAxes(3:4) = handles.ContAxes(3:4)*2^(handles.ZoomStep(2));
end
handles = PlotContinuousData(handles);
handles = PlotRawWFData(handles);
guidata(hObject, handles);


function LoadData(handles, cur_var)

%Display wait message
axes(handles.ContinuousAxes); cla; axis([0 100 0 50]);
text(50, 25, 'Loading data...', 'HorizontalAlignment', 'center');
drawnow;

%Initialize the channel selector (if needed)
t_len = size(handles.DataMatObj, cur_var);
handles.NumChannels = t_len(2);
handles.CurChannel = 1;

if handles.NumChannels > 1,
    set([handles.ChannelText handles.ChannelPopup], 'Enable', 'on');
    chan_str = {};
    for i = 1:handles.NumChannels,
        chan_str{i} = sprintf('Channel %d', i);
    end
    set(handles.ChannelPopup, 'String', chan_str, 'Value', 1);
else
    set([handles.ChannelText handles.ChannelPopup], 'Enable', 'off');
end

%Set variables for viewing continuous time
handles.ContTime = [0 (t_len(1)-1)]/handles.SampleFrequency; %Full time (in seconds)
handles.ContTimeStep = 1/handles.SampleFrequency; %time step (in seconds)
handles.ContTimeViewWidth = 20; %time, in seconds

%Initialize time (x-axis) slider
set(handles.TimeSlider, 'Value', mean(handles.ContTime), ...
    'Min', handles.ContTime(1) + handles.ContTimeViewWidth/2, 'Max', handles.ContTime(2) - handles.ContTimeViewWidth/2, ...
    'SliderStep', [handles.ContTimeViewWidth./diff(handles.ContTime)/10 handles.ContTimeViewWidth./diff(handles.ContTime)], 'Enable', 'on');

%Turn on zoom buttons
set([handles.YZoomInButton handles.YZoomOutButton handles.XZoomInButton handles.XZoomOutButton], 'Enable', 'on');

%Do we have a file for extracting waveforms?
if ~isempty(handles.WFFilename),
    opts = LoadSpikeWF(handles.WFFilename, [], 6);
    handles.SpikeExtract_ThreshType = opts.SpikeExtract_ThreshType;
    handles.SpikeExtract_SpikeThresh = opts.SpikeExtract_SpikeThresh;
    handles.SpikeExtract_ArtifactThresh = opts.SpikeExtract_ArtifactThresh;
    handles.SpikeExtract_WFRange = opts.SpikeExtract_WFRange;
else
    %Use defaults
    handles.SpikeExtract_ThreshType = 0; %std
    handles.SpikeExtract_SpikeThresh = -6;
    handles.SpikeExtract_ArtifactThresh = -25;
    handles.SpikeExtract_WFRange = [-0.4 1.2]; %in ms
end

%Set radio buttons for threshold type
if handles.SpikeExtract_ThreshType == 0,
    set(handles.ThreshTypeRB_V, 'Value', 0);
    set(handles.ThreshTypeRB_STD, 'Value', 1);
elseif handles.SpikeExtract_ThreshType == 1,
    set(handles.ThreshTypeRB_V, 'Value', 1);
    set(handles.ThreshTypeRB_STD, 'Value', 0);
end
set([handles.ThreshTypeRB_V handles.ThreshTypeRB_STD handles.SpikeThreshEdit handles.ArtifactThreshEdit ...
    handles.ReExtractWFButton handles.SpikeThresholdText handles.ArtifactThresholdText], 'Enable', 'on');
set(handles.SpikeThreshEdit, 'String', sprintf('%4.3e', handles.SpikeExtract_SpikeThresh));
set(handles.ArtifactThreshEdit, 'String', sprintf('%4.3e', handles.SpikeExtract_ArtifactThresh));

%Determine the mean/std of data, avoid first 60 seconds (just to make sure
%filter settles)
% for cur_chan = 1:handles.NumChannels,
%     temp_data = handles.DataMatObj.(cur_var)(round(60/handles.ContTimeStep):(t_len(1) - round(60/handles.ContTimeStep)), cur_chan);
%     handles.ContDataMean(cur_chan) = nanmean(temp_data, 1);
%     handles.ContDataSTD(cur_chan) = nanstd(temp_data, [], 1);
% end
temp_data = handles.DataMatObj.(cur_var)(round(60/handles.ContTimeStep):(t_len(1) - round(60/handles.ContTimeStep)), 1:handles.NumChannels);
handles.ContDataMean = nanmean(temp_data, 1);
handles.ContDataSTD = nanstd(temp_data, [], 1);
clear('temp_data');

%Display continuous data
handles = PlotContinuousData(handles);

%Display spiking data
handles = PlotRawWFData(handles);
handles = PlotSummaryWFData(handles);

%Update
guidata(handles.TimeSlider, handles);


function handles = PlotContinuousData(handles)

%What variable is currently selected?
contents = cellstr(get(handles.VarPopup,'String'));
cur_var = contents{get(handles.VarPopup,'Value')};

%What time/indices of data to plot?
time_to_plot = get(handles.TimeSlider, 'Value') + handles.ContTimeViewWidth*[-1 1]/2;
time_to_plot(1) = max(time_to_plot(1), 0); time_to_plot(2) = min(time_to_plot(2), handles.ContTime(2));
plot_ind = round((time_to_plot - handles.ContTime(1))./handles.ContTimeStep) + 1;

handles.ContData = handles.DataMatObj.(cur_var)(plot_ind(1):plot_ind(2), handles.CurChannel);
cur_time = time_to_plot(1) + [0:(length(handles.ContData)-1)]*handles.ContTimeStep;

%Plot Continuous data
axes(handles.ContinuousAxes); axis on; cla;
plot(cur_time, handles.ContData); hold on;
plot(cur_time([1 end]), [0 0], 'k-');
if handles.SpikeExtract_ThreshType == 0,
    plot(cur_time([1 end]), handles.ContDataMean(handles.CurChannel) + handles.SpikeExtract_SpikeThresh*handles.ContDataSTD(handles.CurChannel)*[1 1], 'b-');
    plot(cur_time([1 end]), handles.ContDataMean(handles.CurChannel) + handles.SpikeExtract_ArtifactThresh*handles.ContDataSTD(handles.CurChannel)*[1 1], 'r-');
elseif handles.SpikeExtract_ThreshType == 1,
    plot(cur_time([1 end]), handles.SpikeExtract_SpikeThresh*[1 1], 'b-');
    plot(cur_time([1 end]), handles.SpikeExtract_ArtifactThresh*[1 1], 'r-');
end
if ~isempty(handles.ContAxes),
    set(handles.ContinuousAxes, 'XLim', time_to_plot, 'YLim', handles.ContAxes(3:4));
else
    axis tight;
    handles.ContAxes = axis;
end

%Plot histogram of continuous data
axes(handles.HistAxes); axis on; cla;
[n,x] = hist(handles.ContData, handles.ContAxes(3):(diff(handles.ContAxes(3:4))/100):handles.ContAxes(4));
barh(x, log10(n)); hold on;
set(handles.HistAxes, 'YLim', handles.ContAxes(3:4), 'YTickLabel', []);
v = axis;
if handles.SpikeExtract_ThreshType == 0,
    plot(v(1:2), handles.ContDataMean(handles.CurChannel) + handles.SpikeExtract_SpikeThresh*handles.ContDataSTD(handles.CurChannel)*[1 1], 'b-');
    plot(v(1:2), handles.ContDataMean(handles.CurChannel) + handles.SpikeExtract_ArtifactThresh*handles.ContDataSTD(handles.CurChannel)*[1 1], 'r-');
elseif handles.SpikeExtract_ThreshType == 1,
    plot(v(1:2), handles.SpikeExtract_SpikeThresh*[1 1], 'b-');
    plot(v(1:2), handles.SpikeExtract_ArtifactThresh*[1 1], 'r-');
end

function handles = PlotSummaryWFData(handles)

axes(handles.AvgSpikeAxes); cla; axis on; legend(handles.AvgSpikeAxes, 'off');
%Do we have spikes?
if isempty(handles.WFFilename),
    axis([0 40 0 10]);
    text(20, 5, sprintf('No waveforms extracted yet.'), 'HorizontalAlignment', 'center');
    axis off;
    set(handles.AvgSpikeTitle, 'Enable', 'off');
    return;
end
set(handles.AvgSpikeTitle, 'Enable', 'on');

%Extraction options
Nspk = LoadSpikeWF(handles.WFFilename, [], 5);
opts = LoadSpikeWF(handles.WFFilename, [], 6);
%Number of channels
num_chans = round(opts.NumPointsInWF/opts.SamplePoints);

%Get all waveforms
[~, wf, wf_class] = LoadSpikeWF(handles.WFFilename, [1 Nspk], 4);
uniq_class = unique(wf_class);

%Plot each classified type as a seperate waveform (add NaNs in between for
%spacers)
colors = {[0 0 1], [0 1 0], [1 0 0], [0 1 1], [1 0 1], [0.1 0.1 0.1], [1 0.8 0]};
leg_str = {}; lh = [];
for cur_class = 1:length(uniq_class),
    uniq_class(cur_class)
    cur_ind = (wf_class == uniq_class(cur_class));
    
    %Average waveforms
    mean_wf = []; std_wf = []; spk_t = [];
    for cur_chan = 1:num_chans,
        mean_wf = cat(1, mean_wf, squeeze(nanmean(wf(cur_ind, cur_chan, :), 1)), NaN*ones(4, 1));
        std_wf = cat(1, std_wf, squeeze(nanstd(wf(cur_ind, cur_chan, :), [], 1)), NaN*ones(4, 1));
        spk_t = cat(1, spk_t, (cur_chan - 1)*opts.SamplePoints + [1:opts.SamplePoints]', NaN*ones(4, 1));
    end
    th = plot(spk_t, mean_wf, '-', 'Color', colors{mod(cur_class - 1, length(colors)) + 1});
    lh = cat(1, lh, th);
    hold on;
    th = plot(spk_t, mean_wf + std_wf, ':', 'Color', colors{mod(cur_class - 1, length(colors)) + 1});
    lh = cat(1, lh, th);
    th = plot(spk_t, mean_wf - std_wf, ':', 'Color', colors{mod(cur_class - 1, length(colors)) + 1});
    lh = cat(1, lh, th);
    if uniq_class(cur_class) < 0,
        leg_str{cur_class} = sprintf('Noise %d', uniq_class(cur_class));
    elseif uniq_class(cur_class) < 100,
        leg_str{cur_class} = sprintf('MUA %d', uniq_class(cur_class));
    else
        leg_str{cur_class} = sprintf('Unit %d', uniq_class(cur_class));
    end
end
axis tight;
v = axis;
lh = cat(1, lh(2:3:end), lh(1:3:end), lh(3:3:end));
set(handles.AvgSpikeAxes, 'Children', lh);
legend(leg_str);
plot(v(1:2), [0 0], 'k-');
%axis off;


function handles = PlotRawWFData(handles)

axes(handles.RawSpikeAxes); cla; axis on;
%Do we have spikes?
if isempty(handles.WFFilename),
    axis([0 40 0 10]);
    text(20, 5, sprintf('No waveforms extracted yet.'), 'HorizontalAlignment', 'center');
    axis off;
    set(handles.RawSpikeTitle, 'Enable', 'off');
    return;
end
set(handles.RawSpikeTitle, 'Enable', 'on');

%Extraction options
opts = LoadSpikeWF(handles.WFFilename, [], 6);
%Number of channels
num_chans = round(opts.NumPointsInWF/opts.SamplePoints);

%Get waveforms from the current window of time
time_to_plot = get(handles.TimeSlider, 'Value') + handles.ContTimeViewWidth*[-1 1]/2;
time_to_plot(1) = max(time_to_plot(1), 0); time_to_plot(2) = min(time_to_plot(2), handles.ContTime(2));
[~, wf, wf_class] = LoadSpikeWF(handles.WFFilename, time_to_plot, 3);
uniq_class = unique(wf_class);

if isempty(wf),
    axis([0 40 0 10]);
    text(20, 5, sprintf('No waveforms for section.'), 'HorizontalAlignment', 'center');
    axis off;
    set(handles.RawSpikeTitle, 'Enable', 'off');
    return;
end

colors = {[0 0 1], [0 1 0], [1 0 0], [0 1 1], [1 0 1], [0.1 0.1 0.1], [1 0.8 0]};
%Plot waveforms
if size(wf, 1) > handles.WaveformImageCutoff,
    %Plotting as an image
    
    %Concatenate waveforms across channels
    cat_wf = []; spk_t = [];
    for cur_chan = 1:num_chans,
        cat_wf = cat(2, cat_wf, reshape(wf(:, cur_chan, :), [size(wf, 1) size(wf, 3)]), NaN*ones(size(wf, 1), 4));
        spk_t = cat(1, spk_t, (cur_chan - 1)*(opts.SamplePoints+4) + [1:(opts.SamplePoints+4)]');
    end
    cat_wf = cat_wf(:, 1:(end-4), :);
    spk_t = spk_t(1:(end-4));
    
    wf_max = prctile(cat_wf(:), 99);
    wf_min = prctile(cat_wf(:), 1);
    y_scale = linspace(wf_min, wf_max, handles.NumImageSteps);
    x_scale = linspace(spk_t(1), spk_t(end), handles.ImageTimeInterpScale*length(spk_t));
    cur_img = zeros(handles.NumImageSteps, handles.ImageTimeInterpScale*length(spk_t), 3);
    mean_wf = NaN*ones(length(uniq_class), length(spk_t));
    for cur_class = 1:length(uniq_class),
        cur_ind = (wf_class == uniq_class(cur_class));
        if sum(cur_ind) == 1,
            wf_img = hist(repmat(cat_wf(cur_ind, :), [2 1]), y_scale);
            wf_img = wf_img./2;
        else
            wf_img = hist(interp1(spk_t, cat_wf(cur_ind, :)', x_scale, 'cubic')', y_scale);
        end
        mean_wf(cur_class, :) = nanmedian(cat_wf(cur_ind, :), 1);
        cur_img = cur_img + (repmat(wf_img, [1 1 3]).*repmat(reshape(1 - colors{cur_class}, [1 1 3]), [handles.NumImageSteps length(x_scale) 1]));
    end
    cur_img = cur_img./size(cat_wf, 1);
    cur_img = 1 - cur_img;
    cur_img(repmat(all(cur_img == 0, 3), [1 1 3])) = 1;
    cur_img = imadjust(cur_img, stretchlim(cur_img), [0 1], handles.ImageContrastGamma);
    image(x_scale, y_scale, cur_img);
    hold on;
    for cur_class = 1:length(uniq_class), plot(spk_t, mean_wf(cur_class, :), '-', 'Color', colors{cur_class}, 'LineWidth', 3); end
    set(gca, 'YDir', 'Normal', 'YLim', [wf_min wf_max]);
    axis tight; hold on;
    plot(spk_t, zeros(length(spk_t), 1), 'k-');
else
    %Plot individual waveforms
    
    %Concatenate waveforms across channels
    cat_wf = []; spk_t = [];
    for cur_chan = 1:num_chans,
        cat_wf = cat(2, cat_wf, reshape(wf(:, cur_chan, :), [size(wf, 1) size(wf, 3)]), NaN*ones(size(wf, 1), 4));
        spk_t = cat(1, spk_t, (cur_chan - 1)*opts.SamplePoints + [1:opts.SamplePoints]', NaN*ones(4, 1));
    end
    
    leg_str = {}; lh = [];
    for cur_class = 1:length(uniq_class),
        cur_ind = (wf_class == uniq_class(cur_class));
        
        plot(spk_t, cat_wf(cur_ind, :), '-', 'Color', colors{mod(cur_class - 1, length(colors)) + 1});
        hold on;
    end
    axis tight;
    v = axis;
    plot(v(1:2), [0 0], 'k-');
    %axis off;
end


% --- Executes on button press in XZoomInButton.
function XZoomInButton_Callback(hObject, eventdata, handles)
% hObject    handle to XZoomInButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ContTimeViewWidth = handles.ContTimeViewWidth*2^(-handles.ZoomStep(1));
handles = PlotContinuousData(handles);
handles = PlotRawWFData(handles);
guidata(hObject, handles);

% --- Executes on button press in XZoomOutButton.
function XZoomOutButton_Callback(hObject, eventdata, handles)
% hObject    handle to XZoomOutButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.ContTimeViewWidth = handles.ContTimeViewWidth*2^(handles.ZoomStep(1));
handles = PlotContinuousData(handles);
handles = PlotRawWFData(handles);
guidata(hObject, handles);


% --- Executes on selection change in ChannelPopup.
function ChannelPopup_Callback(hObject, eventdata, handles)
% hObject    handle to ChannelPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChannelPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChannelPopup

handles.CurChannel = get(hObject, 'Value');

%Display continuous data
handles = PlotContinuousData(handles);

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ChannelPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChannelPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SpikeThreshEdit_Callback(hObject, eventdata, handles)
% hObject    handle to SpikeThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SpikeThreshEdit as text
%        str2double(get(hObject,'String')) returns contents of SpikeThreshEdit as a double

handles.SpikeExtract_SpikeThresh = str2double(get(hObject, 'String'));
handles = PlotContinuousData(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SpikeThreshEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SpikeThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ThreshTypeRB_STD.
function ThreshTypeRB_STD_Callback(hObject, eventdata, handles)
% hObject    handle to ThreshTypeRB_STD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ThreshTypeRB_STD
if get(hObject,'Value'),
    set(handles.ThreshTypeRB_V, 'Value', 0);
    handles.SpikeExtract_ThreshType = 0;
    handles = PlotContinuousData(handles);
    guidata(hObject, handles);
end


% --- Executes on button press in ThreshTypeRB_V.
function ThreshTypeRB_V_Callback(hObject, eventdata, handles)
% hObject    handle to ThreshTypeRB_V (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ThreshTypeRB_V
if get(hObject,'Value'),
    set(handles.ThreshTypeRB_STD, 'Value', 0);
    handles.SpikeExtract_ThreshType = 1;
    handles = PlotContinuousData(handles);
    guidata(hObject, handles);
end

% --- Executes on button press in ReExtractWFButton.
function ReExtractWFButton_Callback(hObject, eventdata, handles)
% hObject    handle to ReExtractWFButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Currently selected variable
contents = cellstr(get(handles.VarPopup,'String'));
cur_var = contents{get(handles.VarPopup,'Value')};

%Size of the current variable
t_len = size(handles.DataMatObj, cur_var);

if isempty(handles.WFFilename),
    %Set waveform filename
    [pathstr, name, ext] = fileparts(handles.Filename);
    handles.WFFilename = sprintf('%s\\%s.wf', pathstr, name);
end

if exist(handles.WFFilename, 'file'),
    button = questdlg('Waveform file already exists. Overwrite?', 'Overwrite Waveform file?', 'Yes', 'No', 'Yes');
    if strcmpi(button, 'No'),
        return;
    end
end        
set(hObject, 'String', 'Extracting...'); guidata(hObject, handles); drawnow;
ExtractTrodeSpikeWF(handles.DataMatObj.(cur_var)(1:t_len(1), 1:t_len(2)), handles.WFFilename, 'Overwrite', 1, ...
    'SpikeExtract_ThreshType', handles.SpikeExtract_ThreshType, 'SpikeExtract_SpikeThresh', handles.SpikeExtract_SpikeThresh, ...
    'SpikeExtract_ArtifactThresh', handles.SpikeExtract_ArtifactThresh);
set(hObject, 'String', 'Re-Extract Waveforms'); guidata(hObject, handles); drawnow;



function ArtifactThreshEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ArtifactThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ArtifactThreshEdit as text
%        str2double(get(hObject,'String')) returns contents of ArtifactThreshEdit as a double

handles.SpikeExtract_ArtifactThresh = str2double(get(hObject, 'String'));
handles = PlotContinuousData(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ArtifactThreshEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ArtifactThreshEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on mouse press over axes background.
function AvgSpikeAxes_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to AvgSpikeAxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = PlotSummaryWFData(handles);
guidata(hObject, handles);


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over AvgSpikeTitle.
function AvgSpikeTitle_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to AvgSpikeTitle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = PlotSummaryWFData(handles);
guidata(hObject, handles);