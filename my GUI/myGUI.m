function varargout = myGUI(varargin)
% MYGUI MATLAB code for myGUI.fig
%      MYGUI, by itself, creates a new MYGUI or raises the existing
%      singleton*.
%
%      H = MYGUI returns the handle to a new MYGUI or the handle to
%      the existing singleton*.
%
%      MYGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MYGUI.M with the given input arguments.
%
%      MYGUI('Property','Value',...) creates a new MYGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before myGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to myGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help myGUI

% Last Modified by GUIDE v2.5 27-Feb-2017 11:02:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @myGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @myGUI_OutputFcn, ...
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


% --- Executes just before myGUI is made visible.
function myGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to myGUI (see VARARGIN)

% Choose default command line output for myGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% display text in terminal:
handles.teminalText.String = 'please choose CPG case and load data';

% % disable all irrelevant button groups:
handles.singleNNCB.Enable = 'off';
handles.ourMoECB.Enable = 'off';
handles.paperMoECB.Enable = 'off';
handles.hidNueronNum_slider.Enable = 'off';
handles.Experts_Num_slider.Enable = 'off';
handles.epochsNumSlider.Enable = 'off';
handles.hiddenNeuSliderText.Enable = 'off';
handles.NumOfExpertsSliderText.Enable = 'off';
handles.numOfEpochsSliderText.Enable = 'off';
handles.trainPB.Enable = 'off';
handles.shufflePB.Enable = 'off';
handles.resultsTable.Enable = 'off';
handles.perf_over_epochPB.Enable = 'off';
handles.gatePerfOverEpochPB.Enable = 'off';
handles.MoE_regressionGraphPB.Enable = 'off';
handles.NNregressionPB.Enable = 'off';
handles.toGraph0.Enable = 'off';
handles.toGraph1.Enable = 'off';
handles.toGraph2.Enable = 'off';
handles.toGraph3.Enable = 'off';
handles.refreshGrpahsPB.Enable = 'off';
handles.viewFitting2DCB.Enable = 'off';
handles.viewFitting3DCB.Enable = 'off';
handles.toggleToViewNNRB.Enable = 'off';
handles.toggleToViewMoERB.Enable = 'off';
handles.expertNumType.Enable = 'off';
handles.viewNNWeightsPB.Enable = 'off';
% UIWAIT makes myGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = myGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
disp('starting GUI...');

function teminalText_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of teminalText as text
%        str2double(get(hObject,'String')) returns contents of teminalText as a double


% --- Executes during object creation, after setting all properties.
function teminalText_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in twoNsymmRB.
function twoNsymmRB_Callback(hObject, eventdata, handles)
% hObject    handle to twoNsymmRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of twoNsymmRB


% --- Executes on button press in twoNgenRB.
function twoNgenRB_Callback(hObject, eventdata, handles)
% hObject    handle to twoNgenRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of twoNgenRB


% --- Executes on button press in fourNsymmRB.
function fourNsymmRB_Callback(hObject, eventdata, handles)
% hObject    handle to fourNsymmRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fourNsymmRB


% --- Executes on button press in fourNgenRB.
function fourNgenRB_Callback(hObject, eventdata, handles)
% hObject    handle to fourNgenRB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fourNgenRB


% --- Executes on button press in loadData_PB.
function loadData_PB_Callback(hObject, eventdata, handles)
% hObject    handle to loadData_PB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.teminalText.String = 'loading data now...';
targetCells = {'freq'};
fileName = handles.enterFileNameLine.String;
if handles.twoNsymmRB.Value
    parametersCells = {'tau','b','a','s'};
    seqOrder = {'tau','b','s','_','a'};
    myCode1=myCode(fileName,parametersCells,targetCells,seqOrder,2);
    myCode1.sizeOfCPG = 2;
elseif handles.twoNgenRB.Value        
    parametersCells = {'tau','b','w_{12}','w_{21}'};
    seqOrder = {'tau','b','c_1','c_2','w_{12}','w_{21}'};
    myCode1=myCode(fileName,parametersCells,targetCells,seqOrder,2,true);
    myCode1.sizeOfCPG = 2;
elseif handles.fourNsymmRB.Value
    parametersCells = {'tau','b','w_{12}','w_{13}','w_{14}','w_{23}','w_{24}','w_{34}'};
    seqOrder = {'tau','b','c','w_{12}','w_{13}','w_{14}','w_{23}','w_{24}','w_{34}'};
    myCode1=myCode(fileName,parametersCells,targetCells,seqOrder,4);
    myCode1.sizeOfCPG = 4;
elseif handles.fourNsymmRB.Value   
end

%NOTE: this doesn't change with shuffling!
handles.sampleSizeTable.Data{1,1} = size(myCode1.sampl_train,2);
handles.sampleSizeTable.Data{2,1} = size(myCode1.sampl_valid,2);
handles.sampleSizeTable.Data{3,1} = size(myCode1.sampl_test,2);
hObject.UserData = myCode1;

% enable the Training bottun group:
handles.singleNNCB.Enable = 'on';
handles.ourMoECB.Enable = 'on';
handles.paperMoECB.Enable = 'on';
handles.shufflePB.Enable = 'on';

handles.teminalText.String = 'data was loaded';


function enterFileNameLine_Callback(hObject, eventdata, handles)
% hObject    handle to enterFileNameLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of enterFileNameLine as text
%        str2double(get(hObject,'String')) returns contents of enterFileNameLine as a double


% --- Executes during object creation, after setting all properties.
function enterFileNameLine_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enterFileNameLine (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in useDefaultFileCB.
function useDefaultFileCB_Callback(hObject, eventdata, handles)
% hObject    handle to useDefaultFileCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    if handles.twoNsymmRB.Value
        handles.enterFileNameLine.String = 'MatsRandomRes_2Neurons_symm_trainData_narrow_range.mat';
    elseif handles.twoNgenRB.Value        
        handles.enterFileNameLine.String = 'MatsRandomRes_2Neurons_general_1to4_combained.mat';
    elseif handles.fourNsymmRB.Value
%         handles.enterFileNameLine.String = 'MatsRandomRes_4Neurons_symm.mat';
        handles.enterFileNameLine.String = 'MatsRandomRes_4Neurons_symm_all_samples.mat';
    elseif handles.fourNsymmRB.Value
        handles.enterFileNameLine.String = ' ';
    end
end


% --- Executes on button press in singleNNCB.
function singleNNCB_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.teminalText.String = 'chosen single NN for training';
    pause(1);
    handles.hidNueronNum_slider.Enable = 'on';
    handles.epochsNumSlider.Enable = 'on';
    handles.hiddenNeuSliderText.Enable = 'on';
    handles.numOfEpochsSliderText.Enable = 'on';
    handles.trainPB.Enable = 'on';
    handles.ShowNNGUI.Enable = 'on';
    handles.teminalText.String = 'choose addiional methods and/or choose training parameters';
else
    handles.ShowNNGUI.Enable = 'off';
end
if ~get(hObject,'Value') && ~handles.ourMoECB.Value && ~handles.paperMoECB.Value
    % if no botton is enable then deactivate all
    handles.hidNueronNum_slider.Enable = 'off';
    handles.epochsNumSlider.Enable = 'off';
    handles.hiddenNeuSliderText.Enable = 'off';
    handles.numOfEpochsSliderText.Enable = 'off';
    handles.trainPB.Enable = 'off';
end

% --- Executes on button press in ourMoECB.
function ourMoECB_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.teminalText.String = 'chosen our MoE method';
    pause(1);
    handles.hidNueronNum_slider.Enable = 'on';
    handles.Experts_Num_slider.Enable = 'on';
    handles.epochsNumSlider.Enable = 'on';
    handles.hiddenNeuSliderText.Enable = 'on';
    handles.NumOfExpertsSliderText.Enable = 'on';
    handles.numOfEpochsSliderText.Enable = 'on';
    handles.trainPB.Enable = 'on';
    handles.teminalText.String = 'choose addiional methods and/or choose training parameters';
end
if ~get(hObject,'Value') && ~handles.singleNNCB.Value && ~handles.paperMoECB.Value
    % if no botton is enable then deactivate all
    handles.hidNueronNum_slider.Enable = 'off';
    handles.epochsNumSlider.Enable = 'off';
    handles.hiddenNeuSliderText.Enable = 'off';
    handles.numOfEpochsSliderText.Enable = 'off';
    handles.trainPB.Enable = 'off';
end
    
% --- Executes on button press in paperMoECB.
function paperMoECB_Callback(hObject, eventdata, handles)
if get(hObject,'Value')
    handles.teminalText.String = 'chosen papers MoE method';
    pause(1);
    handles.hidNueronNum_slider.Enable = 'on';
    handles.Experts_Num_slider.Enable = 'on';
    handles.epochsNumSlider.Enable = 'on';
    handles.hiddenNeuSliderText.Enable = 'on';
    handles.NumOfExpertsSliderText.Enable = 'on';
    handles.numOfEpochsSliderText.Enable = 'on';
    handles.trainPB.Enable = 'on';
    handles.teminalText.String = 'choose addiional methods and/or choose training parameters';
end
if ~get(hObject,'Value') && ~handles.singleNNCB.Value && ~handles.ourMoECB.Value
    % if no botton is enable then deactivate all
    handles.hidNueronNum_slider.Enable = 'off';
    handles.epochsNumSlider.Enable = 'off';
    handles.hiddenNeuSliderText.Enable = 'off';
    handles.numOfEpochsSliderText.Enable = 'off';
    handles.trainPB.Enable = 'off';
end

% --- Executes on button press in shufflePB.
function shufflePB_Callback(hObject, eventdata, handles)
% hObject    handle to shufflePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.teminalText.String = 'shuffling the samples...';
if ~isempty(handles.loadData_PB.UserData)
    myCode1 = handles.loadData_PB.UserData;
    myCode1 = myCode1.shuffle_samples();
    handles.loadData_PB.UserData = myCode1;
end
handles.teminalText.String = 'samples shuffled';

% --- Executes on slider movement.
function hidNueronNum_slider_Callback(hObject, eventdata, handles)
% hObject    handle to hidNueronNum_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = get(hObject,'Value');
sliderValue = floor(sliderValue);
set(hObject,'Value',sliderValue);
set(handles.hiddenNeuSliderText,'String',num2str(sliderValue));

% --- Executes during object creation, after setting all properties.
function hidNueronNum_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hidNueronNum_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function hiddenNeuSliderText_Callback(hObject, eventdata, handles)
% hObject    handle to hiddenNeuSliderText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hiddenNeuSliderText as text
%        str2double(get(hObject,'String')) returns contents of hiddenNeuSliderText as a double

% --- Executes during object creation, after setting all properties.
function hiddenNeuSliderText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hiddenNeuSliderText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function Experts_Num_slider_Callback(hObject, eventdata, handles)
% hObject    handle to Experts_Num_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = get(hObject,'Value');
sliderValue = floor(sliderValue);
set(hObject,'Value',sliderValue);
set(handles.NumOfExpertsSliderText,'String',num2str(sliderValue));

% --- Executes during object creation, after setting all properties.
function Experts_Num_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Experts_Num_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function NumOfExpertsSliderText_Callback(hObject, eventdata, handles)
% hObject    handle to NumOfExpertsSliderText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NumOfExpertsSliderText as text
%        str2double(get(hObject,'String')) returns contents of NumOfExpertsSliderText as a double


% --- Executes during object creation, after setting all properties.
function NumOfExpertsSliderText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NumOfExpertsSliderText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function epochsNumSlider_Callback(hObject, eventdata, handles)
% hObject    handle to epochsNumSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
sliderValue = get(hObject,'Value');
sliderValue = floor(sliderValue);
set(hObject,'Value',sliderValue);
set(handles.numOfEpochsSliderText,'String',num2str(sliderValue));

% --- Executes during object creation, after setting all properties.
function epochsNumSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to epochsNumSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function numOfEpochsSliderText_Callback(hObject, eventdata, handles)
% hObject    handle to numOfEpochsSliderText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numOfEpochsSliderText as text
%        str2double(get(hObject,'String')) returns contents of numOfEpochsSliderText as a double


% --- Executes during object creation, after setting all properties.
function numOfEpochsSliderText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numOfEpochsSliderText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in ShowNNGUI.
function ShowNNGUI_Callback(hObject, eventdata, handles)

% Hint: get(hObject,'Value') returns toggle state of ShowNNGUI


% --- Executes on button press in trainPB.
function trainPB_Callback(hObject, eventdata, handles)
% hObject    handle to trainPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.loadData_PB.UserData)
    return;
end

myCode1=handles.loadData_PB.UserData;
numOfEpochs = handles.epochsNumSlider.Value;
numOfExperts = handles.Experts_Num_slider.Value;
numOfHiddenNeurons= handles.hidNueronNum_slider.Value;
myCode1.disp_information = false;

if handles.singleNNCB.Value
    disp('training NN...');
    NNtrainingGUI_flag = handles.ShowNNGUI.Value;
    myCode1 = myCode1.Set('NN',numOfHiddenNeurons,numOfEpochs);
    myCode1 = myCode1.trainNN(NNtrainingGUI_flag);
    handles.resultsTable.Data{1,1} = num2str(myCode1.NN.MSE_test_perf);
    handles.resultsTable.Data{2,1} = num2str(myCode1.NN.RsquarTest);
end

if handles.ourMoECB.Value
    disp('training MoE...');
    competetiveflag = 3;
    myCode1 = myCode1.Set('our_MoE',floor(numOfEpochs/5),...
        numOfExperts,numOfHiddenNeurons,[2],5,competetiveflag);
    myCode1 = myCode1.my_MoE_train();
    handles.resultsTable.Data{1,2} = num2str(myCode1.my_MoE_out.Moe_MSE_on_test);
    handles.resultsTable.Data{2,2} = num2str(myCode1.my_MoE_out.RsquarTest);
end

if handles.paperMoECB.Value
    disp('training papers MoE...');
    myCode1 = myCode1.Set('paper_MoE',numOfEpochs,numOfExperts,0.005,0.995);
    myCode1 = myCode1.paper_MoE_train();
    handles.resultsTable.Data{1,3} = myCode1.paper_MoE_out.Moe_perf_over_iter(1,end);
end

hObject.UserData = myCode1;

handles.resultsTable.Enable = 'on';
handles.toGraph0.Enable = 'on';
handles.toGraph1.Enable = 'on';
handles.toGraph2.Enable = 'on';
handles.toGraph3.Enable = 'on';
handles.refreshGrpahsPB.Enable = 'on';
handles.viewFitting2DCB.Enable = 'on';
handles.viewFitting3DCB.Enable = 'on';
if handles.singleNNCB.Value
    handles.NNregressionPB.Enable = 'on';
    handles.toggleToViewNNRB.Enable = 'on';
end
if handles.ourMoECB.Value
    handles.MoE_regressionGraphPB.Enable = 'on';
    handles.perf_over_epochPB.Enable = 'on';
    handles.gatePerfOverEpochPB.Enable = 'on';
    handles.toggleToViewMoERB.Enable = 'on';
    handles.expertNumType.Enable = 'on';
end
handles.viewNNWeightsPB.Enable = 'on';


function ax = whichAxis(handles)
% check on which axis we want to plot
if handles.toGraph0.Value
    ax = handles.graph0;
elseif handles.toGraph1.Value
    ax = handles.graph1;
elseif handles.toGraph2.Value
    ax = handles.graph2;
elseif handles.toGraph3.Value
    ax = handles.graph3;
end

% --- Executes on button press in perf_over_epochPB.
function perf_over_epochPB_Callback(hObject, eventdata, handles)
if handles.refreshGrpahsPB.Value
    % if we want to refresh the graphs after retraining
    ax = hObject.UserData;
else
    ax = whichAxis(handles);
    hObject.UserData = ax; %which graph this plot was on last
end
obj = handles.trainPB.UserData;
iterNum = 1:obj.numOfIteretions;
% total MSE error over #iteration
plot(ax,iterNum,obj.my_MoE_out.Moe_perf_over_iter,'b-o'); hold on;
ax.XLabel.String='#iteration'; ax.YLabel.String='MoE MSE error';
ax.Title.String='total MSE error over #iteration';
ax.UserData = hObject;
hold off



% --- Executes on button press in gatePerfOverEpochPB.
function gatePerfOverEpochPB_Callback(hObject, eventdata, handles)
if handles.refreshGrpahsPB.Value
    % if we want to refresh the graphs after retraining
    ax = hObject.UserData;
else
    ax = whichAxis(handles);
    hObject.UserData = ax; %which graph this plot was on last
end 
obj = handles.trainPB.UserData;
iterNum = 1:obj.numOfIteretions;
% gateNet perf over #interation
plot(ax,iterNum,obj.my_MoE_out.gateTraniData.gateNN_perf_vec,'-o'); hold on
ax.XLabel.String='#iteration'; ax.YLabel.String='performance [MSE]';
ax.Title.String='gateNet perf (MSE) over #interation';
ax.UserData = hObject;
hold off


% --- Executes on button press in MoE_regressionGraphPB.
function MoE_regressionGraphPB_Callback(hObject, eventdata, handles)
% hObject    handle to MoE_regressionGraphPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.refreshGrpahsPB.Value
    % if we want to refresh the graphs after retraining
    ax = hObject.UserData;
else
    ax = whichAxis(handles);
    hObject.UserData = ax; %which graph this plot was on last
end
cla(ax);
ax.NextPlot= 'add';
obj = handles.trainPB.UserData;
switch obj.expertCount
    case {2,3} % in case of small number of expert, make colors clear:
        colors = [1,0,0;0,1,0;0,0,1];
    otherwise
        colors = rand(obj.expertCount,3);
end
gateOut = obj.my_MoE_out.gateNet(obj.sampl_test);
ouputs = obj.my_MoE_out.out_from_test;
targets = obj.targ_test;
[g_max,g_max_ind] = max(gateOut,[],1);
for j=1:obj.expertCount
    for i=1:size(ouputs,2)
        if g_max_ind(1,i) == j
            if g_max(1,i) > 0.5
                plot(ax,targets(1,i),ouputs(1,i),'k-o','MarkerFaceColor',colors(j,:));
            else
                plot(ax,targets(1,i),ouputs(1,i),'k-o');
            end
        end
    end
end
ax.XLabel.String='target'; ax.YLabel.String='output';
ax.Title.String={'regression graph with different color for every dominant expert';'empty circle mean g<0.5'};
ax.UserData = hObject;
ax.NextPlot= 'replace';

% --- Executes on button press in NNregressionPB.
function NNregressionPB_Callback(hObject, eventdata, handles)
% hObject    handle to NNregressionPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function expertNumType_Callback(hObject, eventdata, handles)
% hObject    handle to expertNumType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of expertNumType as text
%        str2double(get(hObject,'String')) returns contents of expertNumType as a double


% --- Executes during object creation, after setting all properties.
function expertNumType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to expertNumType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in viewNNWeightsPB.
function viewNNWeightsPB_Callback(hObject, eventdata, handles)
% hObject    handle to viewNNWeightsPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.refreshGrpahsPB.Value
    % if we want to refresh the graphs after retraining
    ax = hObject.UserData;
else
    ax = whichAxis(handles);
    hObject.UserData = ax; %which graph this plot was on last
end
axes(ax); %make 'ax' the current axis
obj = handles.trainPB.UserData;
if handles.toggleToViewNNRB.Value
    tempNet = obj.NN.net;
    title_temp = 'NN weights';
elseif handles.toggleToViewMoERB.Value
    expertNum = str2double(get(handles.expertNumType,'String'));
    tempNet = obj.my_MoE_out.expertsNN{1,expertNum};
    title_temp = ['expert #',num2str(expertNum)];
end

parametersCells = obj.inputsNames;
weightsInput = tempNet.IW;  
weightsInput = cell2mat(weightsInput);
weightsOutput = tempNet.LW;
weightsOutput = cell2mat(weightsOutput);
bias = tempNet.b;
bias = cell2mat(bias);

weights = horzcat(diag(weightsOutput)*weightsInput);  % multiply the output weights with the neurons outputs
[x,y] = meshgrid(1:size(weights,2),1:size(weights,1));   %# Create x and y coordinates for the strings
HiddenNumCells = num2cell((1:size(weights,1)));
textStrings = num2str(weights(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding

imagesc(abs(weights));
colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                         %#   black and lower values are white)
ax.Title.String = title_temp;
hStrings = text(x(:),y(:),textStrings(:),'HorizontalAlignment','center');
midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(abs(weights(:)) > midValue,1,3);%# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca,'XTick',1:size(weights,2),...     %# Change the axes tick marks
        'XTickLabel',parametersCells,...  %#   and tick labels
        'YTick',1:size(weights,1),...
        'YTickLabel',HiddenNumCells,...
        'TickLength',[0 0]);
ax.UserData = hObject;
 

% --- Executes on button press in refreshGrpahsPB.
function refreshGrpahsPB_Callback(hObject, eventdata, handles)
% hObject    handle to refreshGrpahsPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% gatePerfOverEpochPB_Callback(handles.gatePerfOverEpochPB,eventdata,handles)
if ~isempty(handles.graph0.UserData)
    handles.graph0.UserData.Callback(handles.graph0.UserData,eventdata);
end
if ~isempty(handles.graph1.UserData)
    handles.graph1.UserData.Callback(handles.graph1.UserData,eventdata);
end
if ~isempty(handles.graph2.UserData)
    handles.graph2.UserData.Callback(handles.graph2.UserData,eventdata);
end
if ~isempty(handles.graph3.UserData)
    handles.graph3.UserData.Callback(handles.graph3.UserData,eventdata);
end


% --- Executes on button press in viewFitting2DCB.
function viewFitting2DCB_Callback(hObject, eventdata, handles)
% unselect the other one
handles.viewFitting3DCB.Value = false;
handles.chooseOneMore.Enable = 'off';
handles.text9.Enable = 'off';
handles.ViewFittingPB.Enable = 'on';

% --- Executes on button press in viewFitting3DCB.
function viewFitting3DCB_Callback(hObject, eventdata, handles)
% unselect the other one
handles.viewFitting2DCB.Value = false;
handles.chooseOneMore.Enable = 'on';
handles.text9.Enable = 'on';
handles.ViewFittingPB.Enable = 'on';

% --- Executes on selection change in fittingViewOptions.
function fittingViewOptions_Callback(hObject, eventdata, handles)
% Hints: contents = cellstr(get(hObject,'String')) returns fittingViewOptions contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fittingViewOptions


% --- Executes during object creation, after setting all properties.
function fittingViewOptions_CreateFcn(hObject, eventdata, handles)
% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [param,chosenParam] = chooseParam(objectHandle,handles)
% get the right CPG parameter from the list for fitting plot
contents = cellstr(get(objectHandle,'String'));
chosenParam = contents{get(objectHandle,'Value')};
obj = handles.trainPB.UserData;
switch chosenParam
    case 'tau'
        param = obj.sampl_test(strcmp('tau',obj.inputsNames),:);
    case 'b'
        param = obj.sampl_test(strcmp('b',obj.inputsNames),:);
    case 'a'
        param = obj.sampl_test(strcmp('a',obj.inputsNames),:);
    case 's'
        param = obj.sampl_test(strcmp('s',obj.inputsNames),:);
    case 'c1'
        param = obj.sampl_test(strcmp('c_1',obj.inputsNames),:);
    case 'c2'
        param = obj.sampl_test(strcmp('c_2',obj.inputsNames),:);
    case 'c3'
        param = obj.sampl_test(strcmp('c_3',obj.inputsNames),:);
    case 'c4'
        param = obj.sampl_test(strcmp('c_4',obj.inputsNames),:);
    case 'w12'
        param = obj.sampl_test(strcmp('w_{12}',obj.inputsNames),:);
    case 'w13'
        param = obj.sampl_test(strcmp('w_{13}',obj.inputsNames),:);
    case 'w14'
        param = obj.sampl_test(strcmp('w_{14}',obj.inputsNames),:);
    case 'w21'
        param = obj.sampl_test(strcmp('w_{21}',obj.inputsNames),:);
    case 'w23'
        param = obj.sampl_test(strcmp('w_{23}',obj.inputsNames),:);
    case 'w24'
        param = obj.sampl_test(strcmp('w_{24}',obj.inputsNames),:);
    case 'w31'
        param = obj.sampl_test(strcmp('w_{31}',obj.inputsNames),:);
    case 'w32'
        param = obj.sampl_test(strcmp('w_{32}',obj.inputsNames),:);
    case 'w34'
        param = obj.sampl_test(strcmp('w_{34}',obj.inputsNames),:);
    case 'w41'
        param = obj.sampl_test(strcmp('w_{41}',obj.inputsNames),:);
    case 'w42'
        param = obj.sampl_test(strcmp('w_{42}',obj.inputsNames),:);
    case 'w43'
        param = obj.sampl_test(strcmp('w_{43}',obj.inputsNames),:);
end

% --- Executes on button press in ViewFittingPB.
function ViewFittingPB_Callback(hObject, eventdata, handles)
obj = handles.trainPB.UserData;
if handles.refreshGrpahsPB.Value
    % if we want to refresh the graphs after retraining
    ax = hObject.UserData;
else
    ax = whichAxis(handles);
    hObject.UserData = ax; %which graph this plot was on last
end
cla(ax);
ax.NextPlot= 'add';

[param1,param1Name] = chooseParam(handles.fittingViewOptions,handles);
if handles.viewFitting3DCB.Value
    % if we want a 3D plot
    [param2,param2Name] = chooseParam(handles.chooseOneMore,handles);
end

if handles.viewFitting2DCB.Value
    scatter(ax,param1,obj.targ_test,'o');
else
    scatter3(ax,param1,param2,obj.targ_test,'o');
end

Legend={'targets'};
if handles.singleNNCB.Value
    output_from_NN = obj.NN.out_from_test;
    if handles.viewFitting2DCB.Value
        scatter(ax,param1,output_from_NN,'d');
    else
        scatter3(ax,param1,param2,output_from_NN,'d');
    end
    Legend{end+1} = 'NN estimation';
end
if handles.ourMoECB.Value
    output_from_ourMoE = obj.my_MoE_out.out_from_test;
    if handles.viewFitting2DCB.Value
        scatter(ax,param1,output_from_ourMoE,'d');
    else
        scatter3(ax,param1,param2,output_from_ourMoE,'d');
    end
    Legend{end+1} = 'our MoE estimation';
end

if handles.twoNsymmRB.Value
    % if 2N CPG than show Matsuoka's estimation 
    freq_Matsuoka_est = zeros(1,size(obj.targ_test,2));
    tau = obj.sampl_test(strcmp('tau',obj.inputsNames),:);
    b = obj.sampl_test(strcmp('b',obj.inputsNames),:);
    a = obj.sampl_test(strcmp('a',obj.inputsNames),:);
    T = 5.*tau;
    for j=1:size(obj.targ_test,2);
        freq_Matsuoka_est(1,j) = obj.MatsuokaEstimation(tau(1,j),T(1,j),b(1,j),a(1,j));
    end
    
    if handles.viewFitting2DCB.Value
        scatter(ax,param1,freq_Matsuoka_est,'x');
    else
        scatter3(ax,param1,param2,freq_Matsuoka_est,'x');
    end
    Legend{end+1}='Matsuoka est';
end

ax.XLabel.String=param1Name;   grid on;
ax.Title.String=['frequency over ',param1Name];
if handles.viewFitting3DCB.Value
    ax.YLabel.String=param2Name;
    ax.ZLabel.String='freq [Hz]';
    ax.Title.String=['frequency over ',param1Name, 'and ',param2Name];
    rotate3d(ax,'on');
else
    view(ax,0,90);
    ax.YLabel.String='freq [Hz]';
end
legend(ax,Legend);
ax.UserData = hObject;

% --- Executes on selection change in chooseOneMore.
function chooseOneMore_Callback(hObject, eventdata, handles)
% hObject    handle to chooseOneMore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns chooseOneMore contents as cell array
%        contents{get(hObject,'Value')} returns selected item from chooseOneMore


% --- Executes during object creation, after setting all properties.
function chooseOneMore_CreateFcn(hObject, eventdata, handles)
% hObject    handle to chooseOneMore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

