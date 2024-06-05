function gen_3d_vol(config_dir)
    % Read input configuration 
    config_file = fullfile(config_dir,'config.txt');
    [opt, ~, ~, symbolicDefs] = readConfig(config_file);
    output_dir = fullfile(config_dir,'output');
    if ~exist(output_dir, 'dir')
       mkdir(output_dir)
    end
    [cell_summary, mask_slice_labels, Mask_3D] = gen_cell_summary(config_dir);

    header = {'z_id','cell_id','centroid_x','centroid_y','cell_area'};
    cell_summary = [header; num2cell(cell_summary)];
    writecell(cell_summary,fullfile(output_dir,'cell_summary.csv'),'Delimiter','comma');
    save(fullfile(output_dir,'mask_view.mat'),'mask_slice_labels','cell_summary','Mask_3D');

    exePath = fullfile(pwd,'gentracemap.exe');
    exePath = ['"' exePath '"'];
    output_dir = ['"' output_dir '"'];
    commandComp = {exePath,output_dir};
    commandStr = strjoin(commandComp," ");
    system(commandStr);
end