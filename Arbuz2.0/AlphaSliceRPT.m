function varargout = AlphaSliceRPT(varargin)
% ALPHASLICERPT M-file for AlphaSliceRPT.fig
% Co-Registration GUI and plug-ins

% University of Chicago Medical Center
% Department of Radiation and Cellular Oncology
% B. Epel    2007

% Last Modified by GUIDE v2.5 22-Aug-2019 14:28:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
  'gui_Singleton',  gui_Singleton, ...
  'gui_OpeningFcn', @AlphaSliceRPT_OpeningFcn, ...
  'gui_OutputFcn',  @AlphaSliceRPT_OutputFcn, ...
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
function AlphaSliceRPT_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.hh = varargin{1};
handles.options.SliceDir = 'Z';
handles.options.images = {};

% Update handles structure
guidata(hObject, handles);
pbUpdateImageList_Callback(hObject, eventdata, handles)

% --------------------------------------------------------------------
function varargout = AlphaSliceRPT_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --------------------------------------------------------------------
function ForceRedraw(h)

handles = guidata(h);

if get(handles.cbRedraw, 'Value')
  pbShow_Callback([], [], handles);
end

% --------------------------------------------------------------------
function pbShow_Callback(~, ~, handles)
tic

mode = 'together';
if handles.cbShowAll.Value, mode = 'separately'; end

% find slices to show
Slices = eval(get(handles.eSlices, 'String'));
nSlices = length(Slices);
nImages = length(handles.options.images);

plist = {'Ashow','data','FullName','Color','Bbox','Mask'};
find_list = cell(nImages, 1);
for ii=1:nImages
  res = arbuz_FindImage(handles.hh, 'master', 'Name', handles.options.images{ii}.name, plist);
  find_list{ii} = res{1};
end

% create figure layout
FigN = str2double(get(handles.eFigureN, 'String'));
figure(FigN); clf;
set(FigN,'Name',[arbuz_get(handles.hh, 'FILENAME'),' -- ',get(handles.eBoundsX, 'String'), ';', get(handles.eBoundsY, 'String'), ';', get(handles.eSlices, 'String')])
switch mode
    case 'together', ax_coordinates = epr_CalcAxesPos(1, nSlices, [0.0025 0.005], [0.01 0.04]);
    case 'separately', ax_coordinates = epr_CalcAxesPos(nImages, nSlices, [0.0025 0.005], [0.01 0.04]);
end
h = zeros(1, nSlices);

% This is the bounding box
boundsX = eval(get(handles.eBoundsX, 'String'));
boundsY = eval(get(handles.eBoundsY, 'String'));
% Slices - Z of the slice space
switch handles.options.SliceDir
  case 'X', [x,y,z] = meshgrid(Slices, boundsX, boundsY);
  case 'Y', [x,y,z] = meshgrid(boundsY, Slices, boundsX);
  case 'Z', [x,y,z] = meshgrid(boundsX, boundsY, Slices);
end
xyzvol = [x(:) y(:) z(:)];

% This is the destination frame
str = get(handles.pmReferenceFrame, 'string');
idx = get(handles.pmReferenceFrame, 'value');
find_listREF = arbuz_FindImage(handles.hh, 'master', 'Name', str{idx}, {'Anative','Ashow'});
if isempty(find_listREF), disp('No reference image.'); return; end
Aref = find_listREF{1}.Ashow \ find_listREF{1}.Anative; % it is inv(A)*B but faster

for jj=1:length(find_list)
  AA = inv(find_list{jj}.Ashow*Aref);  % A2frame
  xyzvt = htransform_vectors(AA, xyzvol);
  
  xt = reshape(xyzvt(:,1), size(x));
  yt = reshape(xyzvt(:,2), size(y));
  zt = reshape(xyzvt(:,3), size(z));
  find_list{jj}.im  = interp3(cast(find_list{jj}.data,'double'), xt,yt,zt, 'linear',0);
  switch handles.options.SliceDir
    case 'X', find_list{jj}.im = permute(find_list{jj}.im, [1,3,2]);
    case 'Y', find_list{jj}.im = permute(find_list{jj}.im, [2,3,1]);
  end
  %
  PostProcessing = safeget(find_list{jj}.Color, 'PostProcessing','');
  if ~isempty(PostProcessing)
    img = find_list{jj}.im; %#ok<NASGU>
    find_list{jj}.im = eval(PostProcessing);
  end
  %
  the_mask = safeget(handles.options.images{jj}, 'mask', 'native');
  the_scale = safeget(handles.options.images{jj}, 'scale', [0,1]);
  switch the_mask
    case {''}, find_list{jj}.im_mask = [];
    case {'native'}
      if ~isempty(find_list{jj}.Mask)
        find_list{jj}.im_mask  = interp3(cast(find_list{jj}.Mask,'double'), xt,yt,zt, 'linear',0) > 0.5;
        switch handles.options.SliceDir
          case 'X', find_list{jj}.im_mask = permute(find_list{jj}.im_mask, [1,3,2]);
          case 'Y', find_list{jj}.im_mask = permute(find_list{jj}.im_mask, [2,3,1]);
        end
        SliceErode=safeget(find_list{jj}.Color, 'SliceErode', 0);
        if SliceErode > 0
          se = strel('disk', SliceErode);
          find_list{jj}.im_mask = imerode(find_list{jj}.im_mask,se);
        end
        if safeget(find_list{jj}.Color, 'SliceLargest', 0)
          for pp=1:nSlices
            [labeled_mask, nMask] = bwlabel(find_list{jj}.im_mask(:,:,pp));
            mask_square = zeros(nMask, 1);
            for ii=1:nMask
              mask_square(ii) = numel(find(labeled_mask == ii));
            end
            [~, max_idx] = max(mask_square);
            find_list{jj}.im_mask(:,:,pp) = labeled_mask == max_idx;
          end
        end
      else
        find_list{jj}.im_mask = [];
      end
    case 'minmax'
      find_list{jj}.im_mask = find_list{jj}.im >= min(the_scale) & find_list{jj}.im <= max(the_scale);
    case 'minmax_alpha'
      the_mask = find_list{jj}.im >= min(the_scale) & find_list{jj}.im <= max(the_scale);
      
      %       find_list{jj}.im_mask = (find_list{jj}.im - min(the_scale))/diff(the_scale);
      
      find_list{jj}.im_mask = sin((find_list{jj}.im - min(the_scale))/diff(the_scale)*pi/2);
      
      find_list{jj}.im_mask(~the_mask) = 0;
    otherwise
      % get the mask
      find_mask = arbuz_FindImage(handles.hh, 'all', 'Name', the_mask, {'data'});
      pidx = -1;
      for ii=1:length(find_mask)
        if strcmp(find_list{jj}.Image, find_mask{ii}.Image), pidx = ii; break; end
      end
      if pidx > 0 && isequal(size(find_mask{pidx}.data), size(find_list{jj}.data))
        if ~isempty(find_list{jj}.Mask)
          find_mask{pidx}.data = find_mask{pidx}.data & find_list{jj}.Mask;
        end
        find_list{jj}.im_mask  = interp3(cast(find_mask{pidx}.data,'double'), xt,yt,zt, 'linear',0) > 0.5;
        switch handles.options.SliceDir
          case 'X', find_list{jj}.im_mask = permute(find_list{jj}.im_mask, [1,3,2]);
          case 'Y', find_list{jj}.im_mask = permute(find_list{jj}.im_mask, [2,3,1]);
        end
      else
        find_list{jj}.im_mask = [];
      end
  end
  
  HighCutOff = safeget(find_list{jj}.Color, 'HighCutOff', 0.95);
  LowCutOff = safeget(find_list{jj}.Color, 'LowCutOff', 0.05);
  CutOffAbsoluteScale = safeget(find_list{jj}.Color, 'CutOffAbsoluteScale', 0);
  
  mx =  cast(max(find_list{jj}.im(:)), 'double');
  
  find_list{jj}.scale = [LowCutOff,HighCutOff];
  find_list{jj}.clevel = safeget(find_list{jj}.Color, 'ContourThreshold', 0.5);
  
  if ~CutOffAbsoluteScale
    find_list{jj}.scale = find_list{jj}.scale * mx;
    find_list{jj}.clevel = find_list{jj}.clevel * mx;
  end
  
  % generate contours
  the_contour = safeget(handles.options.images{jj}, 'contour', '');
  how_to_mask = safeget(handles.options.images{jj}, 'how_to_mask', 1);
  contour_level = safeget(handles.options.images{jj}, 'contour_level', [0.5,0,0]);
  contour_index = safeget(handles.options.images{jj}, 'contour_index', [true,false,false]);
  contour_level = contour_level(contour_index);
  
  image_for_cc = [];
  find_list{jj}.cbar_contours = {};
  if ~isempty(the_contour)
    idxContour = find(contour_index);
    nLevels = numel(idxContour);
    if nLevels == 0, break; end
    find_list{jj}.cc = cell(nSlices, nLevels);
    switch the_contour
      case 'self'
        image_for_cc = find_list{jj}.im;
      otherwise
        % get the mask
        find_mask = arbuz_FindImage(handles.hh, 'all', 'Name', the_contour, {'data'});
        pidx = -1;
        for ii=1:length(find_mask)
          if strcmp(find_list{jj}.Image, find_mask{ii}.Image), pidx = ii; break; end
        end
        if pidx > 0 && isequal(size(find_mask{pidx}.data), size(find_list{jj}.data))
          if ~isempty(find_list{jj}.Mask)
            find_mask{pidx}.data = find_mask{pidx}.data & find_list{jj}.Mask;
          end
          image_for_cc  = interp3(cast(find_mask{pidx}.data,'double'), xt,yt,zt, 'linear',0) > 0.5;
          switch handles.options.SliceDir
            case 'X', image_for_cc = permute(image_for_cc, [1,3,2]);
            case 'Y', image_for_cc = permute(image_for_cc, [2,3,1]);
          end
        end
    end
  end
  
  if ~isempty(image_for_cc)
    find_list{jj}.clevel = contour_level;
    find_list{jj}.c_width = ...
      safeget(handles.options.images{jj}, 'contour_linewidth', 0.5);
    find_list{jj}.c_color = ...
      safeget(handles.options.images{jj},'contour_colors', [1,1,1;1,1,1;1,1,1]);
    find_list{jj}.c_color = find_list{jj}.c_color(contour_index,:);
    
    % colorbar settings
    find_list{jj}.cbar_contours = cell(1, length(find_list{jj}.clevel));
    for kk=1:length(find_list{jj}.clevel)
      find_list{jj}.cbar_contours{kk}.level = find_list{jj}.clevel(kk);
      find_list{jj}.cbar_contours{kk}.color = find_list{jj}.c_color(kk,:);
      find_list{jj}.cbar_contours{kk}.linewidth = find_list{jj}.c_width;
    end
    
    image_for_cc = double(image_for_cc);
    if ~isempty(find_list{jj}.im_mask)
      switch how_to_mask
        case 1, image_for_cc(~(find_list{jj}.im_mask > 0)) = 1E6;
        case 2, image_for_cc(~(find_list{jj}.im_mask > 0)) = -1E6;
      end
    end
    
    for pp=1:nSlices
      if any(any(image_for_cc(:,:,pp)))
        if numel(contour_level) > 1
%           cmask  = image_for_cc(:,:,pp);
%           target = cmask < contour_level(clev);
          
          find_list{jj}.cc{pp} = contourc(boundsX, boundsY, image_for_cc(:,:,pp),...
            contour_level);
        else
          find_list{jj}.cc{pp} = contourc(boundsX, boundsY, image_for_cc(:,:,pp),...
            contour_level*[1,1]);
        end
      end
    end
  end
end

% draw all

for jj=1:nImages
  for ii = 1:nSlices
    hh = ii;
    gen_axis = jj==1;
    if strcmp(mode,'separately') == 1, hh = ii+(jj-1)*nSlices; gen_axis=true; end
    % create axis
    if gen_axis
      current_ax = axes('Position', ax_coordinates(hh,:));
      h(hh) = current_ax;
      title(h(ii), sprintf('%s: %4.1f', handles.options.SliceDir, Slices(ii)));
    else
      current_ax = h(ii);
    end
    
    hold(current_ax, 'on')
    the_colormap = safeget(handles.options.images{jj}, 'colormap', 'jet');
    the_scale = safeget(handles.options.images{jj}, 'scale', [0,1]);
    rgb = arbuz_ind2rgb(find_list{jj}.im(:,:,ii), the_colormap, the_scale);
    hhh = image(boundsX, boundsY, rgb, 'Parent', current_ax, 'CDataMapping','scaled');
    set(current_ax, 'Box', 'on', 'DataAspectRatio', [1,1,1])
    set(current_ax, 'yTick', []);
    alpha = safeget(handles.options.images{jj}, 'alpha', 0.5);
    if strcmp(mode,'separately') == 1, alpha=1; end
    axis(current_ax, 'image');
    
    label_color = [1,1,1];
    if ~isempty(find_list{jj}.im_mask)
      set(hhh,'AlphaData',alpha*find_list{jj}.im_mask(:,:,ii));
      
    label_color = [0,0,0];
    else
      alpha_mask = alpha*ones(size(find_list{jj}.im(:,:,ii)));
      set(hhh,'AlphaData',alpha_mask);
    end
    
    % label
    if ii==1 && handles.cbShowLabels.Value == 1
      text(0.05, 0.1, sprintf('%s',handles.options.images{jj}.name), 'Parent', current_ax, 'units', 'normalized', 'Color', label_color)    
    end
    
    % colorbar
    the_colorbar = safeget(handles.options.images{jj}, 'colorbar', 1);
    cbar_ax = get(current_ax,'Position');
    switch the_colorbar
      case 2
        cbar_ax = cbar_ax + [4.8*cbar_ax(3)/6, 1*cbar_ax(4)/6, -5*cbar_ax(3)/6, -2*cbar_ax(4)/6];
        colorbar_ax = axes('Position', cbar_ax);
        ShowColorbarFIG(colorbar_ax, 'fixed',  handles.options.images{jj}.scale, 'text', 'center',...
          'colormap', the_colormap, 'contours', find_list{jj}.cbar_contours)
      case 3
        cbar_ax = cbar_ax + [0.2*cbar_ax(3)/6, 1*cbar_ax(4)/6, -5*cbar_ax(3)/6, -2*cbar_ax(4)/6];
        colorbar_ax = axes('Position', cbar_ax);
        ShowColorbarFIG(colorbar_ax, 'fixed', handles.options.images{jj}.scale, 'text', 'center', ...
          'colormap', the_colormap, 'contours', find_list{jj}.cbar_contours)
    end
  end
end

% draw contours
for jj=1:nImages
    for ii = 1:nSlices
        for kk=1:nImages
            hh = h(ii);
            if strcmp(mode,'separately') == 1, hh = h(ii+(kk-1)*nSlices); end
            if isfield(find_list{jj}, 'cc')
                contour_data=find_list{jj}.cc{ii};
                ll=1;
                while ll < size(contour_data,2)
                    p = contour_data(2,ll); lev = contour_data(1,ll);
                    r1 = contour_data(:,ll+1:ll+p); ll = ll+1+p;
                    lev_idx = find_list{jj}.clevel == lev;
                    lev_color = find_list{jj}.c_color(lev_idx,:);
                    plot(r1(1,:),r1(2,:),'Parent', hh, ...
                        'Color', lev_color, 'Linewidth', find_list{jj}.c_width);
                end
            end
        end
    end
end

% set(h(:), 'XTick', [], 'YTick', [])
guidata(handles.figure1, handles)
toc

% --- Executes on button press in cbShowAll.
function cbShowAll_Callback(hObject, eventdata, handles)
% hObject    handle to cbShowAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbShowAll

% --------------------------------------------------------------------
function cache = move_to_cache(im, x, y, data, AA)

cache.x = x;
cache.y = y;
cache.im = data;
cache.AA = AA;
cache.ImageIdx = im.ImageIdx;
cache.ProxyIdx = im.ProxyIdx;

% --------------------------------------------------------------------
function ret = is_nochange(im1, im2)

ret = 1;
if isempty(im1) || isempty(im2) || ...
    im1.ImageIdx ~= im2.ImageIdx || im1.ProxyIdx ~= im2.ProxyIdx || ...
    ~isequal(im1.AA, im2.AA)
  ret = 0;
end

% --------------------------------------------------------------------
function pbBringUp_Callback(hObject, eventdata, handles)
FigN = str2double(get(handles.eFigureN, 'String'));
figure(FigN);

% --------------------------------------------------------------------
function pbUpdateImageList_Callback(hObject, eventdata, handles)
[~, str] = GetImageList(handles);

set(handles.pmReferenceFrame, 'string', str);

% --------------------------------------------------------------------
function listbox1_Callback(hObject, eventdata, handles)
sel = get(handles.listbox1, 'value');
if sel < 0, return; end;

alpha = safeget(handles.options.images{sel}, 'alpha', 0.5);
set(handles.eAlpha, 'string', num2str(alpha));
the_mask = safeget(handles.options.images{sel}, 'mask', '');
set(handles.eMask, 'string', the_mask);
scale = safeget(handles.options.images{sel}, 'scale', [0,1]);
set(handles.eImageMin, 'string', num2str(scale(1)));
set(handles.eImageMax, 'string', num2str(scale(2)));

the_colormap = safeget(handles.options.images{sel}, 'colormap', 'jet');
str = get(handles.pmColormap, 'string');
for ii=1:length(str)
  if strcmp(str{ii}, the_colormap)
    set(handles.pmColormap, 'value', ii);
    break;
  end
end
the_colorbar = safeget(handles.options.images{sel}, 'colorbar', 1);
set(handles.pmColorbar, 'Value', the_colorbar);

the_contour = safeget(handles.options.images{sel}, 'contour', '');
set(handles.eContour, 'string', the_contour);
contour_level = safeget(handles.options.images{sel}, 'contour_level', [0.5,0,0]);
contour_index = safeget(handles.options.images{sel}, 'contour_index', [true,false,false]);
contour_controls = [handles.eContourLevel1, handles.eContourLevel2, handles.eContourLevel3];
contour_color_controls = [handles.pbContourColor1, handles.pbContourColor2, handles.pbContourColor3];
contour_color = safeget(handles.options.images{sel},'contour_colors', [1,1,1;1,1,1;1,1,1]);
for ii=1:3
  if contour_index(ii)
    set(contour_controls(ii), 'string', num2str(contour_level(ii)));
  else
    set(contour_controls(ii), 'string', '');
  end
  set(contour_color_controls(ii), 'BackgroundColor', contour_color(ii,:));
end
the_linewidth = safeget(handles.options.images{sel}, 'contour_linewidth', 0.5);
set(handles.eLinewidth, 'string', num2str(the_linewidth));
how_to_mask = safeget(handles.options.images{sel}, 'how_to_mask', 1);
set(handles.pmHowToMask, 'Value', how_to_mask);


% --------------------------------------------------------------------
function pbAdd_Callback(hObject, eventdata, handles)
[~, str] = GetImageList(handles);
[sel,isOk] = listdlg('PromptString','Select an image:',...
  'SelectionMode','single',...
  'ListString',str);
if isOk
  handles.options.images{end+1}.name = str{sel};
  UpdateImageList(handles);
  guidata(handles.figure1, handles);
end

% --------------------------------------------------------------------
function UpdateImageList(handles)
str = {};
for ii=1:length(handles.options.images)
  str{end+1} = handles.options.images{ii}.name;
end
val = get(handles.listbox1, 'value');
if isempty(str)
  set(handles.listbox1, 'string', {'Empty list'}, 'value', 1);
else
  set(handles.listbox1, 'string', str, 'value', min(length(str), val));
end
% --------------------------------------------------------------------
function pbRemove_Callback(hObject, eventdata, handles)
sel = get(handles.listbox1, 'value');
if sel < 0, return; end;
handles.options.images = handles.options.images([1:sel-1,sel+1:end]);
UpdateImageList(handles);
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function ImageData_Callback(hObject, eventdata, handles)
sel = get(handles.listbox1, 'value');
if sel < 0, return; end;

alpha = eval(get(handles.eAlpha, 'string'));
handles.options.images{sel}.alpha = alpha;
handles.options.images{sel}.mask = get(handles.eMask, 'string');
handles.options.images{sel}.scale(1) = str2double(get(handles.eImageMin, 'string'));
handles.options.images{sel}.scale(2) = str2double(get(handles.eImageMax, 'string'));
str = get(handles.pmColormap, 'String');
val = get(handles.pmColormap, 'Value');
handles.options.images{sel}.colormap = str{val};
handles.options.images{sel}.colorbar = get(handles.pmColorbar, 'Value');

handles.options.images{sel}.contour = strtrim(get(handles.eContour, 'string'));
contour_controls = [handles.eContourLevel1, handles.eContourLevel2, handles.eContourLevel3];
contour_color_controls = [handles.pbContourColor1, handles.pbContourColor2, handles.pbContourColor3];
for ii=1:3
  val = strtrim(get(contour_controls(ii), 'string'));
  if isempty(val)
    handles.options.images{sel}.contour_index(ii) = false;
  else
    handles.options.images{sel}.contour_index(ii) = true;
    handles.options.images{sel}.contour_level(ii) = str2double(val);
  end
  handles.options.images{sel}.contour_colors(ii,:) = ...
    get(contour_color_controls(ii), 'BackgroundColor');
end
handles.options.images{sel}.contour_linewidth = ...
  str2double(get(handles.eLinewidth, 'string'));
handles.options.images{sel}.how_to_mask = get(handles.pmHowToMask, 'Value');

guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function pbSelectMask_Callback(hObject, eventdata, handles)
sel = get(handles.listbox1, 'value');
if sel < 0, return; end;

find_master    = arbuz_FindImage(handles.hh, 'master', 'Name', handles.options.images{sel}.name, {'SlaveList'});
find_mask    = arbuz_FindImage(handles.hh, find_master{1}.SlaveList, 'ImageType', '3DMASK', {});

str = {'native', 'minmax', 'minmax_alpha'};
for ii=1:length(find_mask)
  str{end+1}= find_mask{ii}.Slave;
end

[res,isOk] = listdlg('PromptString','Select an image:',...
  'SelectionMode','single', 'ListSize',[160 150],...
  'ListString',str);
if isOk
  set(handles.eMask, 'string', str{res});
  ImageData_Callback(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function [find_list3D, str] = GetImageList(handles)
find_list3D    = arbuz_FindImage(handles.hh, 'master', 'ImageType', '3DEPRI', {'SlaveList'});
find_listFIT = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'FITRESULT', {'SlaveList'});
find_list3DPO2 = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'PO2_pEPRI', {'SlaveList'});
find_list3DAMP = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'AMP_pEPRI', {'SlaveList'});
find_list3DMRI = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'MRI', {'SlaveList'});
find_list3DRAW = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'RAW', {'SlaveList'});
find_list3DAMIRA = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'AMIRA3D', {'SlaveList'});
find_list3DDICOM = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'DICOM3D', {'SlaveList'});
find_list3DIDL = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'IDL', {'SlaveList'});
find_list3DBIN = arbuz_FindImage(handles.hh, 'master', 'ImageType', 'BIN', {'SlaveList'});

for ii=1:length(find_list3DPO2), find_list3D{end+1} = find_list3DPO2{ii}; end
for ii=1:length(find_listFIT), find_list3D{end+1} = find_listFIT{ii}; end
for ii=1:length(find_list3DAMP), find_list3D{end+1} = find_list3DAMP{ii}; end
for ii=1:length(find_list3DMRI), find_list3D{end+1} = find_list3DMRI{ii}; end
for ii=1:length(find_list3DRAW), find_list3D{end+1} = find_list3DRAW{ii}; end
for ii=1:length(find_list3DAMIRA), find_list3D{end+1} = find_list3DAMIRA{ii}; end
for ii=1:length(find_list3DDICOM), find_list3D{end+1} = find_list3DDICOM{ii}; end
for ii=1:length(find_list3DIDL), find_list3D{end+1} = find_list3DIDL{ii}; end
for ii=1:length(find_list3DBIN), find_list3D{end+1} = find_list3DBIN{ii}; end

str = {};
for ii=1:length(find_list3D), str{end+1} = find_list3D{ii}.Image; end

[str, idx]=sort(str);
find_list3D = find_list3D(idx);

% for ii=1:length(find_list3D)
%   find_mask    = arbuz_FindImage(handles.hh, find_list3D{ii}.SlaveList, 'ImageType', '3DMASK', {});
%
% end

% --------------------------------------------------------------------
function pbLoad_Callback(hObject, eventdata, handles)
str = arbuz_GetSaveParList(handles.hh, 'AlphaSliceRPT');
if isempty(str), return; end
[sel,isOk] = listdlg('PromptString','Select an image:',...
  'SelectionMode','single',...
  'ListString',str);
if isOk
  set(handles.figure1, 'name', ['AlphaSliceRPT: ', str{sel}]);
  handles.options = arbuz_GetSavePar(handles.hh, str{sel}, 'AlphaSliceRPT');
  UpdateImageList(handles);
  listbox1_Callback(hObject, eventdata, handles);
  
  % update reference frame
  set(handles.eFigureN, 'String', handles.options.FigN);
  
  % This is the bounding box
  set(handles.eSlices, 'String',safeget(handles.options, 'Slices', '-4:4'));
  set(handles.eBoundsX, 'String',handles.options.boundsX);
  set(handles.eBoundsY, 'String',handles.options.boundsY);
  
  SliceDir = ['X', 'Y', 'Z'];
  handles.options.SliceDir = safeget(handles.options, 'SliceDir', 'Z');
  for ii=1:length(SliceDir)
    if strcmp(SliceDir(ii), handles.options.SliceDir)
      set(handles.pmXYZ, 'value', ii);
      break;
    end
  end
  
  % This is the destination frame
  str = get(handles.pmReferenceFrame, 'string');
  for ii=1:length(str)
    if strcmp(str{ii}, handles.options.REFimg)
      set(handles.pmReferenceFrame, 'value', ii);
      break;
    end
  end
  
  guidata(handles.figure1, handles);
end

% --------------------------------------------------------------------
function pbSave_Callback(hObject, eventdata, handles)
save_pars = arbuz_GetSaveParList(handles.hh, 'AlphaSliceRPT');
handles.options.renderer = 'AlphaSlice';

if isempty(save_pars)
  answer=inputdlg({'Parameter''s name'}, 'Input the name', 1, {'AlphaSlice'});
  arbuz_SetSavePar(handles.hh, answer, 'AlphaSliceRPT', handles.options);
else
  save_pars = [{'New name'}, save_pars{:}];
  [sel,isOk] = listdlg('PromptString','Select an image:',...
    'SelectionMode','single',...
    'ListString',save_pars);
  if isOk
    if sel == 1
      new_name=inputdlg({'Parameter''s name'}, 'Input the name', 1, {'AlphaSlice'});
      new_name = new_name{1};
    else
      new_name = save_pars{sel};
    end
    arbuz_SetSavePar(handles.hh, new_name, 'AlphaSliceRPT', handles.options);
  end
end

% --------------------------------------------------------------------
function ReferenceFrame_Callback(hObject, eventdata, handles)
% create figure layout
handles.options.FigN = get(handles.eFigureN, 'String');

% This is the bounding box
handles.options.Slices = get(handles.eSlices, 'String');
handles.options.boundsX = get(handles.eBoundsX, 'String');
handles.options.boundsY = get(handles.eBoundsY, 'String');

idx = get(handles.pmXYZ, 'Value');
SliceDir = ['X', 'Y', 'Z'];
handles.options.SliceDir = SliceDir(idx);

% This is the destination frame
str = get(handles.pmReferenceFrame, 'string');
idx = get(handles.pmReferenceFrame, 'value');
handles.options.REFimg = str{idx};
guidata(handles.figure1, handles);

% --------------------------------------------------------------------
function pbSelectContour_Callback(hObject, eventdata, handles)
sel = get(handles.listbox1, 'value');
if sel < 0, return; end;

find_master    = arbuz_FindImage(handles.hh, 'master', 'Name', handles.options.images{sel}.name, {'SlaveList'});
find_mask    = arbuz_FindImage(handles.hh, find_master{1}.SlaveList, 'ImageType', '3DMASK', {});

str = {'self'};
for ii=1:length(find_mask)
  str{end+1}= find_mask{ii}.Slave;
end

[res,isOk] = listdlg('PromptString','Select an image:',...
  'SelectionMode','single','ListSize',[160 150],...
  'ListString',str);
if isOk
  set(handles.eContour, 'string', str{res});
  ImageData_Callback(hObject, eventdata, handles);
end

% --------------------------------------------------------------------
function pbContourColors_Callback(hObject, eventdata, handles)
contour_color_controls = [handles.pbContourColor1, handles.pbContourColor2, handles.pbContourColor3];
idx = find(contour_color_controls == hObject, 1);
if isempty(idx), return; end
the_color = get(hObject, 'BackgroundColor');
new_color = uisetcolor(the_color);
set(hObject, 'BackgroundColor', new_color);
ImageData_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function pbMinMaxShow_Callback(hObject, eventdata, handles)
sel = get(handles.listbox1, 'value');
if sel < 0, return; end;

image_name = handles.options.images{sel}.name;
res = arbuz_FindImage(handles.hh, 'master', 'Name', image_name, {'data'});
dmin = min(res{1}.data(:));
dmax = max(res{1}.data(:));
set(handles.eImageMin, 'string', num2str(dmin));
set(handles.eImageMax, 'string', num2str(dmax));
handles.options.images{sel}.scale = [dmin, dmax];
guidata(handles.figure1, handles);


% --- Executes on button press in pbUp.
function pbUp_Callback(hObject, eventdata, handles)
sel = get(handles.listbox1, 'value');
if sel < 2, return; end
handles.options.images([sel-1,sel])  = handles.options.images([sel,sel-1]);
guidata(handles.figure1, handles);
UpdateImageList(handles);
listbox1_Callback(hObject, eventdata, handles);

% --- Executes on button press in pbDn.
function pbDn_Callback(hObject, eventdata, handles)
nImages = length(handles.options.images);
sel = get(handles.listbox1, 'value');
if sel+1 > nImages, return; end
handles.options.images([sel,sel+1])  = handles.options.images([sel+1,sel]);
guidata(handles.figure1, handles);
UpdateImageList(handles);
listbox1_Callback(hObject, eventdata, handles);
