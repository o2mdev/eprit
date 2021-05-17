function varargout = ImageTransformPLG(varargin)
% IMAGETRANSFORMPLG M-file for ImageTransformPLG.fig
%      IMAGETRANSFORMPLG, by itself, creates a new IMAGETRANSFORMPLG or raises the existing
%      singleton*.
%
%      H = IMAGETRANSFORMPLG returns the handle to a new IMAGETRANSFORMPLG or the handle to
%      the existing singleton*.
%
%      IMAGETRANSFORMPLG('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGETRANSFORMPLG.M with the given input arguments.
%
%      IMAGETRANSFORMPLG('Property','Value',...) creates a new IMAGETRANSFORMPLG or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageTransformPLG_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageTransformPLG_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageTransformPLG

% Last Modified by GUIDE v2.5 09-Jun-2015 11:36:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageTransformPLG_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageTransformPLG_OutputFcn, ...
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

% --------------------------------------------------------------------
function ImageTransformPLG_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<INUSL>
handles.output = hObject;

% Add handle of calling object
handles.hh = varargin{1};

plist = {'FullName'};
find_list3D    = arbuz_FindImage(handles.hh, 'master', 'ImageType', '3DEPRI', plist);
find_list3DPO2 = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'PO2_pEPRI', plist);
find_list3DAMP = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'AMP_pEPRI', plist);
find_list3DMRI = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'MRI', plist);
find_list3DAMIRA = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'AMIRA3D', plist);
find_list3DDICOM = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'DICOM3D', plist);
find_listFITRESULT = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'FITRESULT', plist);
find_listBIN = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'BIN', plist);

for ii=1:length(find_list3DPO2), find_list3D{end+1} = find_list3DPO2{ii}; end
for ii=1:length(find_list3DAMP), find_list3D{end+1} = find_list3DAMP{ii}; end
for ii=1:length(find_list3DMRI), find_list3D{end+1} = find_list3DMRI{ii}; end
for ii=1:length(find_list3DAMIRA), find_list3D{end+1} = find_list3DAMIRA{ii}; end
for ii=1:length(find_list3DDICOM), find_list3D{end+1} = find_list3DDICOM{ii}; end
for ii=1:length(find_listFITRESULT), find_list3D{end+1} = find_listFITRESULT{ii}; end
for ii=1:length(find_listBIN), find_list3D{end+1} = find_listBIN{ii}; end

handles.find_list = find_list3D;
str = cell(1, length(find_list3D));
for ii=1:length(find_list3D)
  str{ii}=find_list3D{ii}.FullName;
end
set(handles.pmSource, 'String', str);
pos1 = get(handles.pmSource, 'Value');
fill_PopupMenu(handles.pmSourceMask, handles, handles.find_list{pos1}, true);
set(handles.pmDestination, 'String', str);
% Update handles structure
guidata(hObject, handles);

% --------------------------------------------------------------------
function varargout = ImageTransformPLG_OutputFcn(hObject, eventdata, handles)  %#ok<INUSL>
% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function pbDone_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
the_type = get(handles.pmSourceObject, 'value');

pos1 = get(handles.pmSource, 'Value');
pos2 = get(handles.pmDestination, 'Value');

str = get(handles.pmSource, 'String');
image1.Image = str{pos1};
str = get(handles.pmDestination, 'String');
image2.Image = str{pos2}; 
image2.Slave = get(handles.eDestinationMask, 'String');
if the_type == 1
else
  pos3 = get(handles.pmSourceMask, 'Value');
  str = get(handles.pmSourceMask, 'String');
  image1.Slave = str{pos3};
end

pars.cut_off = str2double(get(handles.eCutOff, 'String'));
pars.dilate_factor = str2double(get(handles.eMaskDilate, 'String'));
pars.algorithm = iff(get(handles.pmTransAlgorithm, 'value'), ...
  'reslice', 'voxel_project');

set(handles.eLog, 'string', 'Processing...', 'backgroundcolor', [1 0 0]);
drawnow;

new_image = arbuz_util_transform(handles.hh, image1, image2, pars);

set(handles.eLog, 'string', 'Finished!', 'backgroundcolor', [0 1 0]);
arbuz_AddImage(handles.hh, new_image, new_image.Image);
arbuz_UpdateInterface(handles.hh);
% --------------------------------------------------------------------
function pmSource_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
pos1 = get(handles.pmSource, 'Value');
fill_PopupMenu(handles.pmSourceMask, handles, handles.find_list{pos1}, true);

% --------------------------------------------------------------------
function pmSourceObject_Callback(hObject, eventdata, handles)
pos1 = get(handles.pmSource, 'Value');
fill_PopupMenu(handles.pmSourceMask, handles, handles.find_list{pos1}, true);

% --------------------------------------------------------------------
function fill_PopupMenu(pm, handles, im_name, use_itself)
ProxyList = arbuz_FindImage(handles.hh, {im_name}, '', '', {'SlaveList'});
find_mask3D = arbuz_FindImage(handles.hh, ProxyList{1}.SlaveList, ...
  'ImageType', '3DMASK', {'Name'});

the_type = get(handles.pmSourceObject, 'value');

if the_type == 2, use_itself = false; end

itself_pos = iff(use_itself, 1, 0);
str = cell(length(find_mask3D)+itself_pos, 1);
if use_itself, str{1} = 'itself'; end
for ii=1:length(find_mask3D), str{ii+itself_pos}=find_mask3D{ii}.Name; end
if isempty(str), str = {'None'}; end
  
set(pm, 'String', str, 'Value', 1)
set(handles.eDestinationMask, 'String', str{1});

set(handles.eLog, 'string', 'Ready.', 'backgroundcolor', [1 1 1]);

% --------------------------------------------------------------------
function pmSourceMask_Callback(hObject, eventdata, handles) %#ok<INUSL,DEFNU>
str = get(handles.pmSourceMask, 'String');
val = get(handles.pmSourceMask, 'Value');
set(handles.eDestinationMask, 'String', str{val});
