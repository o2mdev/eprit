function [s,X,pars_out]=read_dicom_image_file(name)

pars_out = dicominfo(name, 'UseDictionaryVR', true);

[fpath, ~, fext] = fileparts(name);
if isempty(fext), fext = '.'; end
all_files = dir(fullfile(fpath, ['*', fext]));
dirFlags = [all_files.isdir];
files_only = all_files(~dirFlags);
all_names = {files_only.name};

% exclude DICOMDIR
idx = cellfun(@(x) not(contains(x, 'DICOMDIR')),all_names);
all_names = all_names(idx);

% [all_names_sorted] = sort(all_names);
if isfield(pars_out, 'Manufacturer') || ...
    contains(pars_out.Manufacturer, 'Agilent')
      % image resolution       
      SliceThickness = pars_out.SliceThickness;
      PixelSpacing = pars_out.PixelSpacing(:);
      
      VS = [PixelSpacing; SliceThickness];
      
      N = length(all_names);
      
      % slice size and datatype
      S = dicomread(name);
      sz = size(S);
      tp = class(S);
    % pre-allocate data
      VT = zeros([sz N], tp);
      X  = zeros([sz N]);
      POS = zeros(N,2);
      % load each slice and its properties
      for i=1:N
        ifname = fullfile(fpath, all_names{i});
%         fprintf('Reading %s\n', ifname);
        VT(:,:,i) = squeeze(dicomread(ifname));
        info = dicominfo(ifname, 'UseDictionaryVR', true);
        if isfield(info, 'ImagePositionPatient')
          POS(i,:) = [info.ImagePositionPatient(3) i];
        else
          POS(i,:) = [i i];
        end
      end
      % resort the slices according to the image position
      POS = sortrows(POS,1);
      for i=1:N
        X(:,:,i) = VT(:,:,POS(i,2));
      end
      
      pars_out.ImagePositionPatient = []; % block use of position for DICOM with wrong position
      s.bbox = [VS(1),sz(1)*VS(1);VS(2),sz(2)*VS(2);SliceThickness,N*SliceThickness];
      s.dims = [sz(:)',N];
      s.origin = sum(s.bbox,2)/2;
      s.pixsize = VS;
      s.pixel_to_world = hmatrix_translate([-1 -1 -1])*hmatrix_scale(s.pixsize'.*[1 1 1])*...
        hmatrix_translate(-transpose(diff(s.bbox, 1, 2))/2);
      return
end

StationName = safeget(pars_out, 'StationName', '');
if contains(StationName, 'PHILIPS-IJI1EMU')
  DICOMfrom = 'MRI3T';
elseif contains(StationName, 'XRAD-UofC')
      DICOMfrom = 'IMRT_Chuck';
elseif isfield(pars_out, 'Manufacturer')
  switch pars_out.Manufacturer
    case 'Molecubes NV'
      DICOMfrom = 'generic-processed';
    case 'Precision X-Ray'
      DICOMfrom = 'IMRT_Chuck';
    otherwise
      DICOMfrom = 'USI';
  end
else
  DICOMfrom = 'IMRT_Chuck';
end

switch DICOMfrom
  case 'IMRT_Chuck'
    [fpath, fname, fext] = fileparts(name);
    
    pos = strfind(fname, '_');
    if isempty(pos), return; end
    
    fname = [fname(1:pos(end)),'%04i', fext];
    
    n1 = double(safeget(pars_out, 'Rows', 1));
    n2 = double(safeget(pars_out, 'Columns', 1));
    nSlice = safeget(pars_out, 'ImagesInAcquisition', 1);
    res12 = safeget(pars_out, 'PixelSpacing', 1);
    resS = safeget(pars_out, 'SliceThickness', 1);
    
    X = zeros(n1, n2, nSlice);
    slice_pos = zeros(nSlice, 3);
    
    for ii=0:(nSlice-1)
      fname_dicom = fullfile(fpath, sprintf(fname, ii));
      X(:,:,ii+1) = dicomread(fname_dicom);
      pars_out = dicominfo(fname_dicom);
      
      if isfield(pars_out, 'ImagePositionPatient')
        slice_pos(ii+1, :) = pars_out.ImagePositionPatient(:)';
      end
    end
    if mean(diff(slice_pos(:,3))) < 0
      X = flip(X, 3);
    end
    
    s.bbox = [1*res12(1),n1*res12(1);1*res12(2),n2*res12(2);1*resS,nSlice*resS];
    s.dims = [n1,n2,nSlice];
    s.origin = sum(s.bbox,2)/2;
    s.pixsize = [res12(1); res12(2); resS];
    s.pixel_to_world = hmatrix_translate([-1 -1 -1])*hmatrix_scale(s.pixsize'.*[1 1 1])*...
      hmatrix_translate(-diff(s.bbox', [], 1)/2);
  case 'MRI3T'
    X = dicomread(name);
    X = flip(squeeze(X(:,:,1,:)),1);
    
    h = figure(1);
    imagesc(sum(X, 3)); axis image
    title('Select the area of interest.');
    a = gca;
    k = waitforbuttonpress;
    point1 = a.CurrentPoint;    % button down detected
    finalRect = rbbox;          % return figure units
    point2 = a.CurrentPoint;    % button down detected
    close(h);
    
    point1 = fix(point1(1,1:2));            % extract x and y
    point2 = fix(point2(1,1:2));
    yy = sort([point1(1),point2(1)]);
    xx = sort([point1(2),point2(2)]);
    X = X(xx(1):xx(2), yy(1):yy(2), :);
    
    vx = 0.234375;
    a = inputdlg({'Resolution X';'Resolution Y'},'Resolution',1,{num2str(vx);num2str(vx)});
    if ~isempty(a)
      vx = str2double(a{1});
    end
    
    sl = pars_out.SpacingBetweenSlices;
    
    s.bbox = [1*vx,size(X,1)*vx; 1*vx,size(X,2)*vx; 1*sl,size(X,3)*sl];
    s.dims = size(X);
    s.origin = sum(s.bbox,2)/2;
    s.pixsize = [vx; vx; sl];
    s.pixel_to_world = hmatrix_translate([-1 -1 -1])*hmatrix_scale(s.pixsize'.*[1 1 1])*...
      hmatrix_translate(-diff(s.bbox', [], 1)/2);
  case 'generic-processed' % Molecubes
    info = dicominfo(name);
    SliceThickness = info.SliceThickness;
    NumberOfFrames = info.NumberOfFrames;
    X = dicomread(name);
    X = double(squeeze(X));
    n = size(X);
    
    s.bbox = [SliceThickness,n(1)*SliceThickness;SliceThickness,n(2)*SliceThickness;SliceThickness,n(3)*SliceThickness];
    s.dims = n;
    s.origin = sum(s.bbox,2)/2;
    s.pixsize = [SliceThickness; SliceThickness; SliceThickness];
    s.pixel_to_world = hmatrix_translate([-1 -1 -1])*hmatrix_scale(s.pixsize'.*[1 1 1])*...
      hmatrix_translate(-diff(s.bbox', [], 1)/2);

  case 'USI'
    X = dicomread(name);
    X = flip(squeeze(X(:,:,1,:)),1);
    
    Regions = pars_out.SequenceOfUltrasoundRegions.Item_1;
    yy = [Regions.RegionLocationMinX0, Regions.RegionLocationMaxX1];
    xx = [Regions.RegionLocationMinY0, Regions.RegionLocationMaxY1];
    
%     h = figure;
%     imagesc(sum(X, 3)); axis image
%     title('Select the area of interest.');
%     a = gca;
%     k = waitforbuttonpress;
%     point1 = a.CurrentPoint;    % button down detected
%     finalRect = rbbox;          % return figure units
%     point2 = a.CurrentPoint;    % button down detected
%     close(h);
%     point1 = fix(point1(1,1:2));            % extract x and y
%     point2 = fix(point2(1,1:2));
%     yy = sort([point1(1),point2(1)]);
%     xx = sort([point1(2),point2(2)]);

    
    X = X(xx(1):xx(2), yy(1):yy(2), :);
    
    vx = pars_out.PixelSpacing(1);
%     vx = pars_out.PixelSpacing(1);
    dx = 25/246;
    
    a = inputdlg('Input step size', 'USI parameters', 1, {num2str(dx)});
    if ~isempty(a)
      dx = str2double(a{1});
    end
    
    s.bbox = [1*vx,size(X,1)*vx; 1*vx,size(X,2)*vx; 1*dx,size(X,3)*dx];
    s.dims = size(X);
    s.origin = sum(s.bbox,2)/2;
    s.pixsize = [vx; vx; dx];
    s.pixel_to_world = hmatrix_translate([-1 -1 -1])*hmatrix_scale(s.pixsize'.*[1 1 1])*...
      hmatrix_translate(-diff(s.bbox', [], 1)/2);

    % im.Size = size(im.a).*[vx,vx,sl];
    % ibGUI(im)
end
