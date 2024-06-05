function [cell_summary, mask_slice_labels, Mask_3D] = gen_cell_summary(config_dir)
    % Read input configuration 
    config_file = fullfile(config_dir,'config.txt');
    [opt, ~, ~, symbolicDefs] = readConfig(config_file);
    hand_corrected_mask_dir = config_dir;
    file_names = dir(hand_corrected_mask_dir);
    dirFlags = [file_names.isdir];
    sub_folders = file_names(dirFlags);
    sub_folders(ismember( {sub_folders.name}, {'.', '..', 'output'})) = [];
    sub_folder_names = {sub_folders.name};
    r = opt.x_px;
    c = opt.y_px;
    x_min = opt.x_placode_min;
    x_max = opt.x_placode_max;
    y_min = opt.y_placode_min;
    y_max = opt.y_placode_max;
    try
       max_cell_area = opt.max_cell_size;
    catch exception
       max_cell_area = 5000;
    end
    try
       min_cell_area = opt.min_cell_size;
    catch exception
       min_cell_area = 50;
    end

    p = length(sub_folder_names);
    Mask_3D = zeros(r,c,p);

    for i = 1:length(sub_folder_names)
        fullpath = fullfile(hand_corrected_mask_dir,num2str(i),'handCorrection.tif');
        Mask_3D(:,:,i) = im2bw(imread(char(fullpath)),0.9);
    end

    % Thicken mask to 3px
    for inp = 1:p;
        inp_mask = uint8(Mask_3D(:,:,inp));
        for inr = 3:r-2;
            for inc = 3:c-2;
                if (inp_mask(inr,inc)==1)
                    for i = inr-1:inr+1;
                        for j = inc-1:inc+1;
                            inp_mask(i,j) = 3;
                        end
                    end
                end
            end
        end
        Mask_3D(:,:,inp) = logical(inp_mask);
    end
    Mask_3D = logical(Mask_3D);
    cell_summary = [];
    mask_slice_labels = zeros(r,c,p);
    for pid = 1:p
        L_cent = bwlabeln(imcomplement(Mask_3D(:,:,pid)),4);
        L_cent = L_cent + 100;
        s = regionprops(L_cent, 'centroid');
        disp_map = [];
        for cell_id = 1:length(s)
            cell_id = cell_id + 100;
            cellsToColor = (L_cent==cell_id);
            cell_area_img = zeros(r,c);
            cell_area_img(cellsToColor) = 1;
            cell_area = sum(cell_area_img,'all');
            if(cell_area > 50 && cell_area < 5000)
                e = s(cell_id).Centroid;
                test_p1 = imfill(Mask_3D(:,:,pid), ceil([e(1,2) e(1,1)]));
                test_p2 = imsubtract(test_p1, Mask_3D(:,:,pid));
                cell_centroid_area = sum(test_p2,'all');
                if(e(1)>=x_min && e(1)<=x_max && e(2)>=y_min && e(2)<=y_max)
                    this_cell = [pid cell_id e(1) e(2) cell_centroid_area];
                    cell_summary = [cell_summary;this_cell];
                    disp_map = [disp_map 1];
                else
                    disp_map = [disp_map 0];
                end 
            else
                disp_map = [disp_map 0];
            end
        end
        %figure;
        vislabels(L_cent,disp_map)
        %savefig(fullfile(output_dir,strcat(num2str(pid),'_frame_cell_map.fig')));
        mask_slice_labels(:,:,pid) = L_cent;
        close;
    end
    header = {'z_id','cell_id','centroid_x','centroid_y','cell_area'};
    cell_summary = [header; num2cell(cell_summary)];
    writecell(cell_summary,fullfile(config_dir,'cell_summary.csv'),'Delimiter','comma');
end