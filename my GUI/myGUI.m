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

% Last Modified by GUIDE v2.5 22-Feb-2017 19:35:54

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
        handles.enterFileNameLine.String = 'MatsRandomRes_4Neurons_symm.mat';
    elseif handles.fourNsymmRB.Value
        handles.enterFileNameLine.String = ' ';
    end
end


% --- Executes on button press in singleNNCB.
function singleNNCB_Callback(hObject, eventdata, handles)
% hObject    handle to singleNNCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of singleNNCB


% --- Executes on button press in ourMoECB.
function ourMoECB_Callback(hObject, eventdata, handles)
% hObject    handle to ourMoECB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ourMoECB


% --- Executes on button press in paperMoECB.
function paperMoECB_Callback(hObject, eventdata, handles)
% hObject    handle to paperMoECB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of paperMoECB


% --- Executes on button press in spareCB.
function spareCB_Callback(hObject, eventdata, handles)
% hObject    handle to spareCB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of spareCB


% --- Executes on button press in shufflePB.
function shufflePB_Callback(hObject, eventdata, handles)
% hObject    handle to shufflePB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.loadData_PB.UserData)
    myCode1 = handles.loadData_PB.UserData;
    myCode1 = myCode1.shuffle_samples();
    handles.loadData_PB.UserData = myCode1;
end

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
    myCode1 = myCode1.Set('NN',numOfHiddenNeurons,numOfEpochs);
    myCode1 = myCode1.trainNN(0);
    handles.resultsTable.Data{1,1} = num2str(myCode1.NN.MSE_test_perf);
end

if handles.ourMoECB.Value
    disp('training MoE...');
    competetiveflag = 3;
    myCode1 = myCode1.Set('our_MoE',floor(numOfEpochs/5),...
        numOfExperts,numOfHiddenNeurons,[2],5,competetiveflag);
    myCode1 = myCode1.my_MoE_train();
    handles.resultsTable.Data{1,2} = num2str(myCode1.my_MoE_out.Moe_MSE_on_test);
    
%     plot_myMoE(myCode1,handles)
end

if handles.paperMoECB.Value
    disp('training papers MoE...');
    myCode1 = myCode1.Set('paper_MoE',numOfEpochs,numOfExperts,0.005,0.995);
    myCode1 = myCode1.paper_MoE_train();
    handles.resultsTable.Data{1,3} = myCode1.paper_MoE_out.Moe_perf_over_iter(1,end);
end

if handles.spareCB.Value
    %%%%%%%%%
end

hObject.UserData = myCode1;

function plot_myMoE(obj,handles)
% iterNum = 1:obj.numOfIteretions;
% % total MSE error over #iteration
% plot(handles.graph0,iterNum,obj.my_MoE_out.Moe_perf_over_iter,'b-o');
% xlabel('#iteration'); ylabel('MoE MSE error');
% title('total MSE error over #iteration');
% 
% % gateNet perf over #interation
% plot(handles.graph1,iterNum,obj.my_MoE_out.gateTraniData.gateNN_perf_vec,'-o');
% title('gateNet perf (MSE) over #interation');
% xlabel('#iteretion');   ylabel('performance [crossentropy]');
% 
% plotregression(obj.targ_test,obj.my_MoE_out.out_from_test,'test');

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
elseif handles.toGraph4.Value
    ax = handles.graph4;   
end

% --- Executes on button press in perf_over_epochPB.
function perf_over_epochPB_Callback(hObject, eventdata, handles)
ax = whichAxis(handles);   
obj = handles.trainPB.UserData;
iterNum = 1:obj.numOfIteretions;
% total MSE error over #iteration
plot(ax,iterNum,obj.my_MoE_out.Moe_perf_over_iter,'b-o'); hold on;
ax.XLabel.String='#iteration'; ax.YLabel.String='MoE MSE error';
ax.Title.String='total MSE error over #iteration';
hold off


% --- Executes on button press in gatePerfOverEpochPB.
function gatePerfOverEpochPB_Callback(hObject, eventdata, handles)
ax = whichAxis(handles);   
obj = handles.trainPB.UserData;
iterNum = 1:obj.numOfIteretions;
% gateNet perf over #interation
plot(ax,iterNum,obj.my_MoE_out.gateTraniData.gateNN_perf_vec,'-o'); hold on
ax.XLabel.String='#iteration'; ax.YLabel.String='performance [MSE]';
ax.Title.String='gateNet perf (MSE) over #interation';
hold off


% --- Executes on button press in MoE_regressionGraphPB.
function MoE_regressionGraphPB_Callback(hObject, eventdata, handles)
% hObject    handle to MoE_regressionGraphPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ax = whichAxis(handles);
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
ax.NextPlot= 'replace';

% --- Executes on button press in NNregressionPB.
function NNregressionPB_Callback(hObject, eventdata, handles)
% hObject    handle to NNregressionPB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton19.
function pushbutton19_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
