function plugin_data = arbuz_GetSavePar(hGUI, the_save_name, plugin_name)

% load the project
prj = getappdata(hGUI, 'project');

plugin_data = {};

if isfield(prj, 'saves')
  for ii=1:length(prj.saves)
    if strcmp(prj.saves{ii}.name, the_save_name)
      plugin_data = prj.saves{ii}; break;
    end
  end
end