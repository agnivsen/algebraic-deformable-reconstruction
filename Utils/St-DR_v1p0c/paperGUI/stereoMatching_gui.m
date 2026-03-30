% This function goes with the St-DR package.
% Please cite our paper on this topic that you shall find on my web page if
% you use this package. Adrien Bartoli.

function varargout = stereoMatching_gui(varargin)
% STEREOMATCHING_GUI M-file for stereoMatching_gui.fig
%      STEREOMATCHING_GUI, by itself, creates a new STEREOMATCHING_GUI or raises the existing
%      singleton*.
%
%      H = STEREOMATCHING_GUI returns the handle to a new STEREOMATCHING_GUI or the handle to
%      the existing singleton*.
%
%      STEREOMATCHING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STEREOMATCHING_GUI.M with the given input arguments.
%
%      STEREOMATCHING_GUI('Property','Value',...) creates a new STEREOMATCHING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before stereoMatching_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to stereoMatching_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help stereoMatching_gui

% Last Modified by GUIDE v2.5 22-Nov-2007 15:26:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @stereoMatching_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @stereoMatching_gui_OutputFcn, ...
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


% --- Executes just before stereoMatching_gui is made visible.
function stereoMatching_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to stereoMatching_gui (see VARARGIN)

% Choose default command line output for stereoMatching_gui
handles.output = hObject;

clear global
setappdata(0,'UseNativeSystemDialogs',false); % to get a working selection gui
global sift_th1 sift_th2 epipolar_th matching_th pth
sift_th1 = str2double(get(handles.image1_editSift,'string'));
sift_th2 = str2double(get(handles.image2_editSift,'string'));
epipolar_th = str2double(get(handles.editMatchingTh,'string'));
matching_th = str2double(get(handles.editMatchingTh,'string'));

pth = '';

set(handles.output,'deleteFcn','closeAll_guiStereoMatching');

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes stereoMatching_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = stereoMatching_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in image1_pushbuttonCamera.
function image1_pushbuttonCamera_Callback(hObject, eventdata, handles)
% hObject    handle to image1_pushbuttonCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam1 pth
[fl,pth] = uigetfile('*.mat','choose the calibration file for camera 1',[pth '*']);
load([pth fl])
cam1.R = Rc_1;
cam1.K = KK;
cam1.C = -Rc_1'*Tc_1;
cam1.mat = cam1.K * cam1.R * [eye(3) -cam1.C];
set(handles.output,'deleteFcn','closeAll_guiStereoMatching');
fprintf('Camera 1=\n\t%s\n',[pth fl]);

% --- Executes on button press in image1_pushbuttonImage.
function image1_pushbuttonImage_Callback(hObject, eventdata, handles)
% hObject    handle to image1_pushbuttonImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im1 f1 pth
[fl,pth] = uigetfile('*.*','choose image 1',[pth '*']);
im1 = imread([pth fl]);
fprintf('Image 1=\n\t%s\n',[pth fl]);
f1 = figure('name','Image 1','numbertitle','off');
imshow(im1)
hold on

function image1_editSift_Callback(hObject, eventdata, handles)
% hObject    handle to image1_editSift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image1_editSift as text
%        str2double(get(hObject,'String')) returns contents of image1_editSift as a double
global sift_th1
sift_th1 = str2double(get(handles.image1_editSift,'string'));

% --- Executes during object creation, after setting all properties.
function image1_editSift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image1_editSift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in image1_pushbuttonCompute.
function image1_pushbuttonCompute_Callback(hObject, eventdata, handles)
% hObject    handle to image1_pushbuttonCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DSC1 M1 im1 sift_th1 f1 p1
[M1,DSC1] = sift(rgb2gray(im1),'Threshold',sift_th1);
fprintf('Image 1:\n\t%d sift points\n',size(M1,2));
figure(f1)
if ishandle(p1)
    delete(p1)
end
p1 = plot(M1(1,:),M1(2,:),'y.','markersize',10);
drawnow

%%
% --- Executes on button press in image2_pushbuttonCamera.
function image2_pushbuttonCamera_Callback(hObject, eventdata, handles)
% hObject    handle to image2_pushbuttonCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global cam2 pth
[fl,pth] = uigetfile('*.mat','choose the calibration file for camera 2',[pth '*']);
load([pth fl])
cam2.R = Rc_1;
cam2.K = KK;
cam2.C = -Rc_1'*Tc_1;
cam2.mat = cam2.K * cam2.R * [eye(3) -cam2.C];
fprintf('Camera 2=\n\t%s\n',[pth fl]);

% --- Executes on button press in image2_pushbuttonImage.
function image2_pushbuttonImage_Callback(hObject, eventdata, handles)
% hObject    handle to image2_pushbuttonImage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global im2 f2 pth
[fl,pth] = uigetfile('*.*','choose image 1',[pth '*']);
im2 = imread([pth fl]);
fprintf('Image 2=\n\t%s\n',[pth fl]);
f2 = figure('name','Image 2','numbertitle','off');
imshow(im2)
hold on

function image2_editSift_Callback(hObject, eventdata, handles)
% hObject    handle to image2_editSift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of image2_editSift as text
%        str2double(get(hObject,'String')) returns contents of image2_editSift as a double
global sift_th2
sift_th2 = str2double(get(handles.image2_editSift,'string'));

% --- Executes during object creation, after setting all properties.
function image2_editSift_CreateFcn(hObject, eventdata, handles)
% hObject    handle to image2_editSift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in image2_pushbuttonCompute.
function image2_pushbuttonCompute_Callback(hObject, eventdata, handles)
% hObject    handle to image2_pushbuttonCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global DSC2 M2 im2 sift_th2 f2 p2
[M2,DSC2] = sift(rgb2gray(im2),'Threshold',sift_th2);
fprintf('Image 2:\n\t%d sift points\n',size(M2,2));
figure(f2)
if ishandle(p2)
    delete(p2)
end
p2 = plot(M2(1,:),M2(2,:),'y.','markersize',10);
drawnow

%%
function editMatchingTh_Callback(hObject, eventdata, handles)
% hObject    handle to editMatchingTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMatchingTh as text
%        str2double(get(hObject,'String')) returns contents of editMatchingTh as a double
global matching_th
matching_th = str2double(get(handles.editMatchingTh,'string'));

% --- Executes during object creation, after setting all properties.
function editMatchingTh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMatchingTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editEpipolarTh_Callback(hObject, eventdata, handles)
% hObject    handle to editEpipolarTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editEpipolarTh as text
%        str2double(get(hObject,'String')) returns contents of editEpipolarTh as a double
global epipolar_th
epipolar_th = str2double(get(handles.editMatchingTh,'string'))

% --- Executes during object creation, after setting all properties.
function editEpipolarTh_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editEpipolarTh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbuttonMatch.
function pushbuttonMatch_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global matching_th epipolar_th M1 DSC1 M2 DSC2 m1 dsc1 m2 dsc2 q1 q2 cam1 cam2 h_matching im1 im2
m1 = M1;
dsc1 = DSC1;
m2 = M2;
dsc2 = DSC2;

indm = siftmatch(dsc1,dsc2,matching_th);
m1 = m1(:,indm(1,:));
dsc1 = dsc1(:,indm(1,:));
m2 = m2(:,indm(2,:));
dsc2 = dsc2(:,indm(2,:));

% stereo consistency
c = checkStereo(m1(1:2,:),m2(1:2,:),cam1.mat,cam2.mat,epipolar_th);
m1 = m1(:,c);
dsc1 = dsc1(:,c);
m2 = m2(:,c);
dsc2 = dsc2(:,c);

% points
q1 = m1(1:2,:);
q2 = m2(1:2,:);

fprintf('Matching:\n\t%d matches\n',size(q1,2));

% plot the matches
if ishandle(h_matching)
    figure(h_matching)
else
    h_matching = figure('name','Matching','numbertitle','off');
end
plotMatchesMP(im1,im2,q1,q2,h_matching);
drawnow


% --- Executes on button press in pushbuttonTriangulate.
function pushbuttonTriangulate_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonTriangulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global q1 q2 cam1 cam2 h_3D Q
Q = zeros(3,size(q1,2));
hwb = waitbar(0,'Triangulation Progress');
hAxes = findobj(hwb,'type','axes');
hTitle = get(hAxes,'title');


for i=1:size(q1,2)
    Qtmp = SfM_TriPoint({cam1.mat,cam2.mat}, [q1(:,i) q2(:,i)]);
    Q(:,i) = Qtmp(1:3)/Qtmp(4);
    waitbar(i/size(q1,2),hwb);
    set(hTitle,'string',['Triangulation Progress: ' num2str(i/size(q1,2)*100,'%.0f') ' %']);    
end
close (hwb);
% plot the 3D points cloud
if ishandle(h_3D)
    figure(h_3D)
else
    h_3D = figure('name','3D points','numbertitle','off');
end
plot3(Q(1,:),Q(2,:),Q(3,:),'b.','markersize',10);
axis equal
drawnow

% --- Executes on button press in pushbuttonSave.
function pushbuttonSave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global q1 q2 m1 m2 dsc1 dsc2 Q cam1 cam2 pth

[fl,pth] = uiputfile('*.mat','Save as',[pth 'dataStereo.mat']);

save([pth fl],'q1','q2','m1','m2','dsc1','dsc2','Q','cam1','cam2');
fprintf('Data saved as:\n\t%s\n\n',[pth fl]);


