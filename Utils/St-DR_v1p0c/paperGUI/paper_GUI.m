% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function varargout = paper_GUI(varargin)
%PAPER_GUI M-file for paper_GUI.fig
%      PAPER_GUI, by itself, creates a new PAPER_GUI or raises the existing
%      singleton*.
%
%      H = PAPER_GUI returns the handle to a new PAPER_GUI or the handle to
%      the existing singleton*.
%
%      PAPER_GUI('Property','Value',...) creates a new PAPER_GUI using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to paper_GUI_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      PAPER_GUI('CALLBACK') and PAPER_GUI('CALLBACK',hObject,...) call the
%      local function named CALLBACK in PAPER_GUI.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help paper_GUI

% Last Modified by GUIDE v2.5 19-Nov-2007 17:15:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @paper_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @paper_GUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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


% --- Executes just before paper_GUI is made visible.
function paper_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for paper_GUI
handles.output = hObject;

% GUI initialisation
clear global

global pp sh figure_model figure_3D rulings rulingHandles rulingPositions isDown handlesG

handlesG = handles;
rulings = [];
rulingHandles = [];
rulingPositions = [];
isDown = 0;

set(handles.output,'position',[0.5000 30.1667 73.1667 31.8333])
pp.r = str2double(get(handles.model_editRatio,'string'));
sh = make_shape(pp.r);
figure_model = figure('position',[4 4 443 276],'name','Model Parameters','numbertitle','off');
figure_3D = figure('position',[449 2 573 687],'name','3D Model','numbertitle','off');
hold on

set(handles.output,'deleteFcn','closeAll_guiPaper');

figure(figure_model)
plot(sh.poly(1,:),sh.poly(2,:));
axis off
set(figure_model,'color','w');
axis equal;
set(figure_model,'pointer','fleur')
hold on

set(figure_model,'WindowButtonMotionFcn','moveParameter')
set(figure_model,'WindowButtonDownFcn','downParameter')
set(figure_model,'WindowButtonUpFcn','upParameter')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes paper_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = paper_GUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function model_editRatio_Callback(hObject, eventdata, handles)
% hObject    handle to model_editRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of model_editRatio as text
%        str2double(get(hObject,'String')) returns contents of model_editRatio as a double
global pp sh rulingHandles rulingPositions rulings figure_model
pp.r = str2double(get(handles.model_editRatio,'string'));
sh = make_shape(pp.r);
figure(figure_model)
hold off
plot(sh.poly(1,:),sh.poly(2,:));
axis off
axis equal;
hold on

% plot the parameters
for i=1:size(rulings,2)
    m = al2m(rulings(1:2,i)',sh);
    be = m(:,1) + diff(m,1,2)*(rulings(3,i)/(2*pi)+.5);
    rulingPositions(:,i) = [m(:) ; be];
    rulingHandles(1,i) = plot(m(1,1),m(2,1),'k.','markersize',25);
    rulingHandles(2,i) = plot(m(1,2),m(2,2),'k.','markersize',25);
    rulingHandles(3,i) = plot(be(1),be(2),'k.','markersize',25);
    rulingHandles(4,i) = plot(m(1,:),m(2,:),'k','markersize',25);
end
plotPaper

% --- Executes during object creation, after setting all properties.
function model_editRatio_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_editRatio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function model_editExtra_Callback(hObject, eventdata, handles)
% hObject    handle to model_editExtra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of model_editExtra as text
%        str2double(get(hObject,'String')) returns contents of model_editExtra as a double
plotPaper

% --- Executes during object creation, after setting all properties.
function model_editExtra_CreateFcn(hObject, eventdata, handles)
% hObject    handle to model_editExtra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in model_pushbuttonAdd.
function model_pushbuttonAdd_Callback(hObject, eventdata, handles)
% hObject    handle to model_pushbuttonAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global figure_model fsPoint
figure(figure_model)
fsPoint = 'f'; % first point
set(figure_model,'pointer','cross')
set(figure_model,'WindowButtonMotionFcn','movePointerAdd')
set(figure_model,'WindowButtonDownFcn','selectPointAdd')


% --- Executes on button press in model_pushbuttonDelete.
function model_pushbuttonDelete_Callback(hObject, eventdata, handles)
% hObject    handle to model_pushbuttonDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global figure_model fsPoint
figure(figure_model)
set(figure_model,'WindowButtonMotionFcn','moveDelete')
set(figure_model,'WindowButtonDownFcn','downDelete')



% --- Executes on button press in model_pushbuttonLoad.
function model_pushbuttonLoad_Callback(hObject, eventdata, handles)
% hObject    handle to model_pushbuttonLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global rulings rulingPositions rulingHandles pp figure_model
dir('*_gui.mat');
s = input('name of the file: ','s');
load(s);

set(handles.model_editRatio,'string',num2str(pp.r));
sh = make_shape(pp.r);
figure(figure_model)
hold off
plot(sh.poly(1,:),sh.poly(2,:));
axis off
axis equal;
hold on

% plot the parameters
for i=1:size(rulings,2)
    m = al2m(rulings(1:2,i)',sh);
    be = m(:,1) + diff(m,1,2)*(rulings(3,i)/(2*pi)+.5);
    rulingPositions(:,i) = [m(:) ; be];
    rulingHandles(1,i) = plot(m(1,1),m(2,1),'k.','markersize',25);
    rulingHandles(2,i) = plot(m(1,2),m(2,2),'k.','markersize',25);
    rulingHandles(3,i) = plot(be(1),be(2),'k.','markersize',25);
    rulingHandles(4,i) = plot(m(1,:),m(2,:),'k','markersize',25);
end
plotPaper

warning('the number of extra rulings is not the good one!')

function export_editName_Callback(hObject, eventdata, handles)
% hObject    handle to export_editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of export_editName as text
%        str2double(get(hObject,'String')) returns contents of export_editName as a double


% --- Executes during object creation, after setting all properties.
function export_editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to export_editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in export_pushbuttonExport.
function export_pushbuttonExport_Callback(hObject, eventdata, handles)
% hObject    handle to export_pushbuttonExport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global paper randomPoints_2D randomPoints_3D rulings rulingHandles rulingPositions pp mesh3D mesh2D T
suffix = get(handles.export_editName,'string');
if isempty(dir([suffix '_*.mat']))
    save([suffix '_model.mat'],'paper');
    save([suffix '_randomPoints.mat'],'randomPoints_2D','randomPoints_3D');
    save([suffix '_mesh.mat'],'mesh3D','mesh2D','T');
    save([suffix '_gui.mat'],'rulings','rulingPositions', 'pp');
    fprintf('Data saved.\n')
else
    fprintf('Files exists - change the name.\n')
end

function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in generation_pushbuttonGenerate.
function generation_pushbuttonGenerate_Callback(hObject, eventdata, handles)
% hObject    handle to generation_pushbuttonGenerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global randomPoints_2D
randomPoints_2D = rand(2,str2double(get(handles.generation_editNumber,'string')));
plotPaper

% --- Executes on button press in display_radiobuttonRulings.
function display_radiobuttonRulings_Callback(hObject, eventdata, handles)
% hObject    handle to display_radiobuttonRulings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of display_radiobuttonRulings
plotPaper

% --- Executes on button press in display_radiobuttonPoints.
function display_radiobuttonPoints_Callback(hObject, eventdata, handles)
% hObject    handle to display_radiobuttonPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of display_radiobuttonPoints
plotPaper

% --- Executes on button press in display_radiobuttonMesh.
function display_radiobuttonMesh_Callback(hObject, eventdata, handles)
% hObject    handle to display_radiobuttonMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of display_radiobuttonMesh
plotPaper


function display_editMesh_Callback(hObject, eventdata, handles)
% hObject    handle to display_editMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of display_editMesh as text
%        str2double(get(hObject,'String')) returns contents of display_editMesh as a double
plotPaper

% --- Executes during object creation, after setting all properties.
function display_editMesh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to display_editMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




