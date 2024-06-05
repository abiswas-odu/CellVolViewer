
function show_volume_tif(file_name)

    % Scaling factor for z axis
    z_axis_scale = 10;

    [file_dir,name,ext] = fileparts(file_name);
    
    V = tiffreadVolume(file_name);
    labels = unique(V);
    cmap = jet(size(labels,1));
    
    % Loop over the cells and draw each mesh
    for i=2:size(labels,1)
        A = V;
        A = flip(A,3);
        A(A~=labels(i))=0;
        h = patch(isosurface(A,0));
        h.EdgeColor = 'none';
        h.FaceColor = cmap(i, :);
        h.FaceAlpha = 0.7;
    end
    
    view(3);
    camlight 
    lighting gouraud
    zlim([0 z_axis_scale]);
    savefig(fullfile(file_dir,strcat(name,'_3D_volume.fig')));
end
