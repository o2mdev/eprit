function varargout = IMRTGUI_Launch_RPT(varargin)
% IMRTGUI_LAUNCH_RPT MATLAB code for IMRTGUI_Launch_RPT.fig
%      IMRTGUI_LAUNCH_RPT, by itself, creates a new IMRTGUI_LAUNCH_RPT or raises the existing
%      singleton*.
%
%      H = IMRTGUI_LAUNCH_RPT returns the handle to a new IMRTGUI_LAUNCH_RPT or the handle to
%      the existing singleton*.
%
%      IMRTGUI_LAUNCH_RPT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMRTGUI_LAUNCH_RPT.M with the given input arguments.
%
%      IMRTGUI_LAUNCH_RPT('Property','Value',...) creates a new IMRTGUI_LAUNCH_RPT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before IMRTGUI_Launch_RPT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to IMRTGUI_Launch_RPT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help IMRTGUI_Launch_RPT

% Last Modified by GUIDE v2.5 24-Oct-2016 15:26:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @IMRTGUI_Launch_RPT_OpeningFcn, ...
                   'gui_OutputFcn',  @IMRTGUI_Launch_RPT_OutputFcn, ...
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


% --- Executes just before IMRTGUI_Launch_RPT is made visible.
function IMRTGUI_Launch_RPT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to IMRTGUI_Launch_RPT (see VARARGIN)

% Choose default command line output for IMRTGUI_Launch_RPT
handles.output = hObject;
% Add handle of calling object
handles.hh = varargin{1};
handles.Project_name = arbuz_get(handles.hh, 'FILENAME'); 
IMRTGUI([],handles.Project_name )
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes IMRTGUI_Launch_RPT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = IMRTGUI_Launch_RPT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
