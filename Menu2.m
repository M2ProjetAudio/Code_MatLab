function varargout = Menu2(varargin)
% MENU2 MATLAB code for Menu2.fig
%      MENU2, by itself, creates a new MENU2 or raises the existing
%      singleton*.
%
%      H = MENU2 returns the handle to a new MENU2 or the handle to
%      the existing singleton*.
%
%      MENU2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MENU2.M with the given input arguments.
%
%      MENU2('Property','Value',...) creates a new MENU2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Menu2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Menu2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu2.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Menu2

% Last Modified by GUIDE v2.5 10-Jan-2017 18:52:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Menu2_OpeningFcn, ...
                   'gui_OutputFcn',  @Menu2_OutputFcn, ...
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


% --- Executes just before Menu2 is made visible.
function Menu2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Menu2 (see VARARGIN)

% Choose default command line output for Menu2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Menu2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Menu2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in gauss.
function gauss_Callback(hObject, eventdata, handles)
% hObject    handle to gauss (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gauss


% --- Executes on button press in parole.
function parole_Callback(hObject, eventdata, handles)
% hObject    handle to parole (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of parole


% --- Executes on button press in hyp1.
function hyp1_Callback(hObject, eventdata, handles)
% hObject    handle to hyp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hyp1


% --- Executes on button press in hyp2.
function hyp2_Callback(hObject, eventdata, handles)
% hObject    handle to hyp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hyp2


% --- Executes on button press in hyp3.
function hyp3_Callback(hObject, eventdata, handles)
% hObject    handle to hyp3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hyp3


% --- Executes on button press in run_button.
function run_button_Callback(hObject, eventdata, handles)
% hObject    handle to run_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
runProgramme(hObject, eventdata, handles);


% --- Executes on button press in Afficher_button.
function Afficher_button_Callback(hObject, eventdata, handles)
% hObject    handle to Afficher_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1);
runAffichage(hObject, eventdata, handles);


function nom_Callback(hObject, eventdata, handles)
% hObject    handle to nom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nom as text
%        str2double(get(hObject,'String')) returns contents of nom as a double


% --- Executes during object creation, after setting all properties.
function nom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function runProgramme(hObject, eventdata, handles)
parole=get(handles.parole,'Value')
gauss=get(handles.gauss,'Value')
wav_name=get(handles.nom,'String')
hyp1=get(handles.hyp1,'Value')
hyp2=get(handles.hyp2,'Value')
hyp3=get(handles.hyp3,'Value')
cas_de_figure=1*(gauss*hyp1)+2*gauss*hyp2+3*gauss*hyp3+...
    4*parole*hyp1+5*parole*hyp2+6*parole*hyp3
fps=evalin('base','fps');


run algorithme_final_bloc1.m


save('result_data.mat','outputSignal','Jpos','D','azimuth','maxcrit_exp'...
    ,'thetaArg','theta_mle','duree_son','Nb_Loca','wav_name','fs',...
    'parole','gauss');

function runAffichage(hObject, eventdata, handles)
load result_data
run afficher_bloc1.m
clear resul_data



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
fps=str2double(get(hObject,'String'));
assignin('base','fps',fps);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
fps=10;
assignin('base','fps',fps);


% --- Executes on button press in pushbutton_enregistrer.
function pushbutton_enregistrer_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_enregistrer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
load result_data
run enregistrer_bloc1.m
clear resul_data

