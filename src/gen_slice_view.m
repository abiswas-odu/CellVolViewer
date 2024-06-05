function gen_slice_view(config_dir)
    config_file = fullfile(config_dir,'config.txt');
    % Read input configuration 
    [opt, ~, ~, symbolicDefs] = readConfig(config_file);
    output_dir = fullfile(config_dir,'output');
    if ~exist(output_dir, 'dir')
       fprintf('No previously generated 3D volume found. Please run gen_3d_vol first.\n');
       return;
    end
    % Generate new Masks and cell summaries
    %[new_cell_summary, new_mask_slice_labels, New_Mask_3D] = gen_cell_summary(config_dir);
    %header = {'z_id','cell_id','centroid_x','centroid_y','cell_area'};
    %new_cell_summary = [header; num2cell(new_cell_summary)];
    
    % Load previous environemnt
    load(fullfile(output_dir,'mask_view.mat'));
    [r, c, p] = size(mask_slice_labels);
    
    % Load the cell tracks
    cell_tracks = csvread(fullfile(output_dir,'cell_tracks.csv'),1);
    % Remove the cell_track_id in it exists; gets added by gen_cell_metrics
    while p < size(cell_tracks,2)
        cell_tracks(:,end) = [];
    end
    
    % Generate updated mask and environemnt 
    %[cell_summary, mask_slice_labels] = combine_cells(cell_summary, mask_slice_labels, ...
    %                                        new_cell_summary, new_mask_slice_labels);
    %cell_summary = [header; cell_summary];                                    

    inc_frame_no=1;
    inc_cell_id_list=[]; 
    %inc_cell_id_list=[];
    bg = [1 1 1; 0 0 0; 1 1 0; 0 1 0];
    cmap = distinguishable_colors(size(cell_tracks,1),bg);
    cmap = [0 0 0; 1 1 1;0 1 0; cmap];
    mask_slice_color = zeros(r,c,p);
    inc_cell_track_bitmap = zeros(size(cell_tracks,1));
    % Color each cell track across frames  
    for i=1:size(cell_tracks,1)
        this_track = cell_tracks(i,:);
        cell_track_area = zeros(1,p);
        % Assign a color index
        track_color_index = i+3;
        % Create a slice view of only that cell
        cell_slice_color = zeros(r,c,p);

        % Check cell is in the inclusion list 
        if (isempty(inc_cell_id_list) || ismember(this_track(inc_frame_no),inc_cell_id_list))
            inc_cell_track_bitmap(i)=1;
            % Loop over all the frames and see if the cell exists in the track
            for j=1:length(this_track)
                cell_id = this_track(j);
                if(cell_id ~= 0)
                    L_cent = mask_slice_labels(:,:,j);
                    L_cent(L_cent ~= cell_id) = 0;
                    cell_track_area(j) = sum(L_cent>0,'all');
                    L_cent(L_cent == cell_id) = track_color_index;
                    cell_slice_color(:,:,j)=L_cent;
                end
            end
            mask_slice_color(cell_slice_color == track_color_index) = track_color_index; 
        else
            % Loop over all the frames and see if the cell exists in the track
            for j=1:length(this_track)
                cell_id = this_track(j);
                if(cell_id ~= 0)
                    L_cent = mask_slice_labels(:,:,j);
                    L_cent(L_cent ~= cell_id) = 0;
                    cell_track_area(j) = sum(L_cent>0,'all');
                    L_cent(L_cent == cell_id) = track_color_index;
                    cell_slice_color(:,:,j)=L_cent;
                end
            end
            mask_slice_color(cell_slice_color == track_color_index) = 3; 
        end
    end
    % Color the excluded cells across frames
    % Compute excluded cells to color seperately
    cell_excluded_summary = cell_summary(1,:);
    for j=1:p
        this_frame = cell_tracks(:,j);
        for k=2:length(cell_summary)
            if(cell_summary{k,1}==j)
                cell_id = cell_summary{k,2};
                if ~any(this_frame(:) == cell_id)
                    cell_excluded_summary = [cell_excluded_summary; cell_summary(k,:)];
                    if(cell_id ~= 0)
                        L_cent = mask_slice_labels(:,:,j);
                        L_cent(L_cent ~= cell_id) = 0;
                        L_cent(L_cent == cell_id) = 2;
                        this_frame_slice = mask_slice_color(:,:,j);
                        this_frame_slice(L_cent == 2) = 2;
                        mask_slice_color(:,:,j) = this_frame_slice;
                    end
                end
            end
        end
    end
    
    % Label the cells with frame-cell ids 
    mask_slice_color_number = zeros(r,c,p,3);
    for j=1:p
        this_frame = cell_tracks(:,j);
        % Convert the index image to an RGB so that we can add text labels into it
        InsertedImage = ind2rgb(mask_slice_color(:,:,j),cmap);
        % Loop over all the cells
        for k=2:length(cell_summary)
            if(cell_summary{k,1}==j)
                cell_id = cell_summary{k,2};
                cent_x = cell_summary{k,3};
                cent_y = cell_summary{k,4};
                cell_str = [num2str(j) '-' num2str(cell_id)];
                
                % Insert text into the RGB image
                % White label for tracked images 
                % White 
                if ~any(this_frame(:) == cell_id)
                    InsertedImage = insertText(InsertedImage, [cent_x-10 cent_y-10], cell_str, ...
                            'BoxOpacity', 0.0, 'FontSize', 10, 'TextColor','black');                    
                else
                    InsertedImage = insertText(InsertedImage, [cent_x-10 cent_y-10], cell_str, ...
                            'BoxOpacity', 0.0, 'FontSize', 10, 'TextColor','white');
                end
            end
        end
        mask_slice_color_number(:,:,j,:) = InsertedImage;
    end
    figure;
    sliceViewer(mask_slice_color_number);
    writecell(cell_summary,fullfile(output_dir,'cell_summary.csv'),'Delimiter','comma');
    save(fullfile(output_dir,'mask_view.mat'),'mask_slice_labels','cell_summary','New_Mask_3D');
    writecell(cell_excluded_summary,fullfile(output_dir,'excluded_cells.csv'),'Delimiter','comma');
end

function [cell_summary,mask_slice_labels] = combine_cells(old_cell_summary, old_mask_slice_labels, ...
                                        new_cell_summary, new_mask_slice_labels)

    [r, c, p] = size(old_mask_slice_labels);
    mask_slice_labels = zeros(r,c,p);
    max_cell_ids = zeros(1,p);
    for i=2:length(old_cell_summary)
        z_id = old_cell_summary{i,1};
        if old_cell_summary{i,2} > max_cell_ids(z_id)
            max_cell_ids(z_id) = old_cell_summary{i,2}; 
        end
    end
    %Add a large offset
    max_cell_ids = max_cell_ids + 900;
    cell_summary = [];
    for i=2:length(new_cell_summary)
        new_cell  = new_cell_summary(i,:);
        match_found = 0; 
        for j=2:length(old_cell_summary)
            if (new_cell_summary{i,1}==old_cell_summary{j,1} && ...
                round(new_cell_summary{i,3},1)==round(old_cell_summary{j,3},1) && ...
                round(new_cell_summary{i,4},1)==round(old_cell_summary{j,4},1))
                
                old_cell  = old_cell_summary(j,:);
                z_id = old_cell{1};
                cell_id = old_cell{2};
                cell_summary = [cell_summary; old_cell];
                % Assign the cell_ids
                this_slice = mask_slice_labels(:,:,z_id);
                L_cent = old_mask_slice_labels(:,:,z_id);
                this_slice(L_cent == cell_id) = cell_id;
                mask_slice_labels(:,:,z_id) = this_slice;
                match_found = 1;
            end
        end
        % Adding as new cell
        if match_found == 0
            z_id = new_cell{1};
            cell_id = new_cell{2};
            % Assign the next cell id
            next_cell_id = max_cell_ids(z_id);
            next_cell_id = next_cell_id + 1;
            new_cell{2} = next_cell_id;
            max_cell_ids(z_id) = next_cell_id;
            cell_summary = [cell_summary; new_cell];
            
            % Assign the cell_ids
            this_slice = mask_slice_labels(:,:,z_id);
            L_cent = new_mask_slice_labels(:,:,z_id);
            this_slice(L_cent == cell_id) = next_cell_id;
            mask_slice_labels(:,:,z_id) = this_slice;
            
        end
    end
    % Check if new cells found 
end