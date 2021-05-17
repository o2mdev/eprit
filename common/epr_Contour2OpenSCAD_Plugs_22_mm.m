function epr_Contour2OpenSCAD_Plugs_22_mm(mod_name, fname, vectors, scale, UVShift)

endl = sprintf('\n');
str = ['//',mod_name,'(10);', endl, endl];
str = [str,'module ',mod_name,'(h)', endl];
str = [str,'{',endl];
str = [str,sprintf('%s','translate([',num2str(UVShift(1)),',',num2str(UVShift(2)),',','0])'),endl] ;
str = [str, sprintf(' scale([%g, %g, %g])', scale(1), scale(2), scale(3)),endl];
str = [str, sprintf('rotate([0,0,270])',endl)] ;
str = [str,' union() {', endl];

for k = 1:length(vectors)

  str = [str,'linear_extrude(height=h) polygon([', endl];
  vector = vectors{k};
  
  separator = '';
  for i=1:size(vector,1)
    str = [str,separator,sprintf('[%g,%g]',vector(i,1),vector(i,2))];
    separator = ',';
  end
  str = [str,']);', endl];
end

str = [str,'}',endl];
str = [str,'}',endl];

str = [str,' ',endl]; %add lines for readablity
str = [str,' ',endl];
str = [str,' ',endl];

str = [str,' // special variables for the properties of circular objects ',endl];
str = [str,' $fa = 0.01; ',endl];
str = [str,' $fs = 0.5;',endl];

%  //constants 
%  Plug_height = 10; 
%  Plug_OD = 18.0; 
%  Notch_parameter = 2.55; 
%  
%  //plug with arb hole for hypoxia targeting 
%  module Plug_model(){ 
%  cylinder(d = Plug_OD,h = Plug_height); 
% 
%  translate([-Plug_OD/2-1 ,-Notch_parameter/2,2]) 
%   cube([Notch_parameter, Notch_parameter, 4],true);
% 
% }

str = [str,' //constants ',endl];
str = [str,' Plug_height = 9.8; ',endl];
%Changed plug height 5/2/16 -MM
%str = [str,' Plug_OD = 15; ',endl];
str = [str,' Plug_OD = 22.4; ',endl];
str = [str,' Notch_parameter = 4.3; ',endl];
str = [str,' ',endl];
str = [str,' //plug with arb hole for hypoxia targeting ',endl];
str = [str,' module Plug_model(){ ',endl];
str = [str,' cylinder(d = Plug_OD,h = Plug_height); ',endl];
str = [str,' translate([-Plug_OD/2-1 ,0,2]) ',endl];
str = [str,'  cube([Notch_parameter, Notch_parameter, 4],true); } ',endl];
str = [str,' ',endl];

% str = [str,' mirror([0,0,0]) ',endl];
% str = [str,' translate([(Plug_OD)-4,0,-Plug_height]) ',endl];
% str = [str,' linear_extrude(height = Plug_height, scale = 1.0) ',endl];
% str = [str,' circle(Plug_OD/4,$fn = 6, center = true);',endl];


str = [str,' ',endl];
str = [str,' ',endl];
str = [str,' difference(){ ',endl];
str = [str,' Plug_model();',endl];
str = [str,' translate([0,0,-1]) ',endl];
str = [str,' Plug_cut(20); ',endl];
str = [str,'}',endl];



fid = fopen(fname, 'w+');

if fid == -1
    slashes = strfind(fname,'\');
    path = fname(1:slashes(end))
    status = mkdir(path)
    fid = fopen(fname, 'w+');
end
    
fwrite(fid, str);
fclose(fid);
