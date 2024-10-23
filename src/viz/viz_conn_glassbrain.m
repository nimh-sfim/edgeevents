function conn_display = viz_conn_glassbrain(indat,edg_color,rois)

global CONN_gui; CONN_gui.usehighres=true; %#ok<GVMIS>

if ~isfield(rois,'sph_r')
    disp('making roi sizes according to degree')
    rois.sph_r = normalize(sum(indat), 'range',[2,5]); % Add Radius based on degree
else
    disp('aready provides node sz')
end

% Prepare input for plotting function
% -----------------------------------
disp('++ Opening the 3D Window')
tic
conn_display = conn_mesh_display('','','',rois,indat) ; 
toc
disp('++ Set connections color')
tic
conn_display('con_color',edg_color)
toc
disp('Set ROI Transparency')
tic
conn_display('roi_transparency',0.7)
toc
disp('Set Brain Transparency')
tic
conn_display('brain_transparency',0.1)
toc
disp('Unset Subcorical surface')
tic
conn_display('sub_transparency',0)
toc
conn_display('con_bundling',1)
conn_display('con_width',0.5)
conn_display('con_transparency',.3)
conn_display('view',[0,0,1])
conn_display('background',[1,1,1])
conn_display('material',[]) % Equivalent to selecting Flat
