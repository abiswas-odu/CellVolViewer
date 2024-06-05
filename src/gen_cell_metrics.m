function gen_cell_metrics(config_dir)

    % Specify tracks to display empty means ALL
    inc_cell_tracks=[];
    
    config_file = fullfile(config_dir,'config.txt');
    % Load configuration 
    [opt, ~, ~, symbolicDefs] = readConfig(config_file);
    intrep_frame_count = opt.inprep_frames;
    output_dir = fullfile(config_dir,'output');
    nhood_dialate = 9; % Default if the parameter is missing 
    if isfield(opt,'nhood_dialate')
        nhood_dialate = opt.nhood_dialate;
    end

    % Load environemnt
    load(fullfile(output_dir,'mask_view.mat'));
    [r c p] = size(mask_slice_labels);
    
    % Load the cell tracks
    cell_tracks = csvread(fullfile(output_dir,'cell_tracks.csv'),1);
    % Remove the cell_track_id if it exists; gets added by gen_cell_metrics
    while p < size(cell_tracks,2)
        cell_tracks(:,end) = [];
    end

    % Compute scaling factor
    px_area_scaling = (opt.x_real*opt.y_real)/(opt.x_px*opt.y_px);
    px_vol_scaling = (opt.x_real*opt.y_real*opt.z_real)/(opt.x_px*opt.y_px*(p-1)*intrep_frame_count);

    [cell_ctr frame_ctr] = size(cell_tracks);
    cell_areas = zeros(cell_ctr,frame_ctr+11);
    cell_plot = zeros(cell_ctr,frame_ctr);
    vol_interp_multi_cell = zeros(r,c,intrep_frame_count*p);
    neighbor_cell_count = [];
    neighbor_cell_track_ids = {};
    neighbor_cell_count_frame = [];
    for i=1:size(cell_tracks,1)
        fprintf('Processing cell: %d out of %d.\n', i, size(cell_tracks,1));
        this_track = cell_tracks(i,:);
        vol_interp = zeros(r,c,intrep_frame_count*p);
        Mask_3D_Test1 = zeros(r,c,p);
        track_color_index = i+1;
        valid_cell_plane = [];
        touching_tracks = [];
        touching_tracks_frame = [];
        % Fill cell body in the frame and make everything else 0
        for j = 1:p
            % Goto first plane
            cell_id = this_track(j);
            if(cell_id ~= 0 && cell_is_valid(cell_id, j, cell_summary))
                valid_cell_plane=[valid_cell_plane j];
                touching_tracks_j = find_touching_cells(mask_slice_labels, j, cell_id, cell_tracks, nhood_dialate);
                touching_tracks = [touching_tracks touching_tracks_j];
                touching_tracks_frame = [touching_tracks_frame length(touching_tracks_j)]; 
                L_cent = mask_slice_labels(:,:,j);
                L_cent(L_cent ~= cell_id) = 0;
                L_cent(L_cent == cell_id) = 1;
                curr_plane_area = sum(L_cent>0,'all');
                cell_areas(i,j + 10) = curr_plane_area * px_area_scaling;
                cell_plot(i,j) = curr_plane_area * px_area_scaling;
                Mask_3D_Test1(:,:,j) = L_cent;
                % Obtain the centroid of the cell from cell summary
%                 for k=2:length(cell_summary)
%                     if(cell_summary{k,1}==j && cell_summary{k,2}==cell_id)
%                         cent_x = cell_summary{k,3};
%                         cent_y = cell_summary{k,4};
%                         break;
%                     end
%                 end
%                 % Fill the cell body from the centroid 
%                 test_p1 = imfill(New_Mask_3D(:,:,j), ceil([cent_y cent_x]));
%                 test_p2 = imsubtract(test_p1,New_Mask_3D(:,:,j));
%                 Mask_3D_Test1(:,:,j) = test_p2;
            else
                touching_tracks_frame = [touching_tracks_frame 0]; 
            end
        end
        neighbor_cell_count = [neighbor_cell_count; length(unique(touching_tracks))];
        neighbor_cell_track_ids = [neighbor_cell_track_ids; {unique(touching_tracks)}];
        neighbor_cell_count_frame = [neighbor_cell_count_frame; touching_tracks_frame];
        % Check if any valid cells are found 
        if(isempty(valid_cell_plane))
            continue;
        end
        
        % Get volume centroid
        mid_frame = valid_cell_plane(ceil(end/2));
        mid_frame_cell_id = this_track(mid_frame);
        for k=2:length(cell_summary)
            if(cell_summary{k,1}==mid_frame && cell_summary{k,2}==mid_frame_cell_id)
                cell_areas(i,1) = cell_summary{k,3} * (opt.x_real/opt.x_px);
                cell_areas(i,2) = cell_summary{k,4} * (opt.y_real/opt.y_px);
                break;
            end
        end
        cell_areas(i,3) = mid_frame;
        cell_areas(i,4) = mid_frame_cell_id;
        
        % Get first plane centroid for height calculation 
        first_frame = valid_cell_plane(1);
        first_frame_cell_id = this_track(first_frame);
        for k=2:length(cell_summary)
            if(cell_summary{k,1}==first_frame && cell_summary{k,2}==first_frame_cell_id)
                cell_areas(i,5) = cell_summary{k,3} * (opt.x_real/opt.x_px);
                cell_areas(i,6) = cell_summary{k,4} * (opt.y_real/opt.y_px);
                break;
            end
        end
        cell_areas(i,7) = first_frame;
        % Get last plane centroid for height calculation 
        last_frame = valid_cell_plane(end);
        last_frame_cell_id = this_track(last_frame);
        for k=2:length(cell_summary)
            if(cell_summary{k,1}==last_frame && cell_summary{k,2}==last_frame_cell_id)
                cell_areas(i,8) = cell_summary{k,3} * (opt.x_real/opt.x_px);
                cell_areas(i,9) = cell_summary{k,4} * (opt.y_real/opt.y_px);
                break;
            end
        end
        cell_areas(i,10) = last_frame;
        % Generated interpolated view
        cell_terminated=0;
        for j = 1:p-1
            if(cell_terminated==1)
                break;
            end
            cell_id = this_track(j);
            if(cell_id ~= 0)
                curr_plane = Mask_3D_Test1(:,:,j);
%                 curr_plane_area = sum(curr_plane,'all');
%                 cell_areas(i,j+4) = curr_plane_area * px_area_scaling;
                vol_interp(:,:,(j-1)*intrep_frame_count+1) = curr_plane;

                if(this_track(j+1) ~= 0) % Cell exists in next plane
                    next_plane=Mask_3D_Test1(:,:,j+1);
                    for inr = 2:r;
                        found_curr_indx = find(curr_plane(inr,:)==1);
                        found_next_indx = find(next_plane(inr,:)==1);
                        if(numel(found_curr_indx)>0 && numel(found_next_indx)==0) % cell not yet found in next plane
                            max_pos = found_curr_indx(end);
                            min_pos = found_curr_indx(1);
                            mid_pos = min_pos + (max_pos-min_pos)/2;
                            start_pos = round(linspace(min_pos,mid_pos,9));
                            end_pos = round(linspace(max_pos,mid_pos,9));
                            for k = 2:intrep_frame_count
                                 for inc = start_pos(k-1):end_pos(k-1)
                                    vol_interp(inr,inc,(j-1)*intrep_frame_count+k) = 1;
                                 end
                            end
                        elseif(numel(found_curr_indx)>0 && numel(found_next_indx)>0)
                            start_pos = round(linspace(found_curr_indx(1),found_next_indx(1),9));
                            end_pos = round(linspace(found_curr_indx(end),found_next_indx(end),9));
                            for k = 2:intrep_frame_count
                                 for inc = start_pos(k-1):end_pos(k-1)
                                    vol_interp(inr,inc,(j-1)*intrep_frame_count+k) = 1;
                                 end
                            end
                        end
                    end
                else % Cell should end
                    max_row=1;
                    min_row=r;
                    % Find max min row range
                    for inr = 2:r;
                        found_indx = find(curr_plane(inr,:)==1);
                        if(numel(found_indx)>0)
                            if(inr>max_row)
                                max_row=inr;
                            end
                            if(inr<min_row)
                                min_row = inr;
                            end
                        end
                    end
                    mid_row = min_row + (max_row-min_row)/2;
                    start_row = round(linspace(min_row,mid_row,9));
                    start_row = start_row(1:end-uint8(length(start_row)/2));
                    end_row = round(linspace(max_row,mid_row,9));
                    end_row = end_row(1:end-uint8(length(end_row)/2));
                    for k = 2:uint8(intrep_frame_count/2);
                        for inr = start_row(k-1):end_row(k-1);
                            found_indx = find(curr_plane(inr,:)==1);
                            if(numel(found_indx)>0)
                                max_pos = found_indx(numel(found_indx));
                                min_pos = found_indx(1);
                                mid_pos = min_pos + (max_pos-min_pos)/2;
                                start_pos = round(linspace(min_pos,mid_pos,9));
                                end_pos = round(linspace(max_pos,mid_pos,9));
                                 for inc = start_pos(k-1):end_pos(k-1);
                                    vol_interp(inr,inc,(j-1)*intrep_frame_count+k) = 1;
                                 end
                            end
                        end
                    end
                    cell_terminated=1;
                end
            end
        end
        cell_areas(i,frame_ctr+11) = sum(vol_interp,'all') * px_vol_scaling;
        for j = 1:intrep_frame_count*p;
            for inr = 1:r;
                for inc = 1:c;
                    if vol_interp(inr,inc,j) == true;
                        vol_interp_multi_cell(inr,inc,j) = track_color_index;
                    end
                end
            end
        end
    end

    fprintf('Processing cell volumes complete.\n');

    % Compute distances of the cells from the center of the image
    
    x_cent = (opt.x_placode_max + opt.x_placode_min)/2 * (opt.x_real/opt.x_px);
    y_cent = (opt.y_placode_max + opt.y_placode_min)/2 * (opt.y_real/opt.y_px);
    x_um_centroid = cell_areas(:,1);
    y_um_centroid = cell_areas(:,2);
    center_distances = sqrt((x_cent - x_um_centroid) .^ 2 + (y_cent - y_um_centroid) .^ 2);
    
    % Compute hieght of the cells 
    x_first = cell_areas(:,5);
    y_first = cell_areas(:,6);
    z_first = cell_areas(:,7);
    x_last = cell_areas(:,8);
    y_last = cell_areas(:,9);
    z_last = cell_areas(:,10);
    hieght_distances = sqrt((x_first - x_last) .^ 2 + (y_first - y_last) .^ 2 + (z_first - z_last) .^ 2);
    
    % Compute ratio of apical:basal surface 
    apical_basal_ratio = [];
    for i=1:length(cell_areas)
        apical_surf_frame = cell_areas(i,7);
        apical_surf_area = cell_areas(i,apical_surf_frame+10);
        if(apical_surf_frame+1 < p)
            apical_surf_area_2 = cell_areas(i,apical_surf_frame+1+10);
            apical_surf_area = (apical_surf_area + apical_surf_area_2)/2;
        end
        
        basal_surf_frame = cell_areas(i,10);
        basal_surf_area = cell_areas(i,basal_surf_frame+10);
        if(basal_surf_frame-1 > 0)
            basal_surf_area_2 = cell_areas(i,basal_surf_frame-1+10);
            basal_surf_area = (basal_surf_area + basal_surf_area_2)/2;
        end
        apical_basal_ratio = [apical_basal_ratio (apical_surf_area/basal_surf_area)];
    end
    % Delete the first last frame coordinates. 
    cell_areas(:,10) = [];
    cell_areas(:,7) = [];
    
    % Add to output file
    cell_ids = 1:1:size(cell_areas,1);
    cell_areas = [cell_areas cell_ids' apical_basal_ratio' center_distances hieght_distances neighbor_cell_count];
    % Sort by distance from center 
    cell_areas = sortrows(cell_areas,size(cell_areas,2)-2);

    header = {'centroid_x(um)','centroid_y(um)','centroid_z(um)','centroid_z_cell_id','centroid_first_x(um)','centroid_first_y(um)','centroid_last_x(um)','centroid_last_y(um)',};
    for j = 1:p
        header_txt_area = strcat('area_frame_',num2str(j),'(um)');
        header = [header header_txt_area];
    end
    header = [header 'cell_volume(um)' 'cell_track_id' 'apical_basal_area_ratio(um)' 'center_distance(um)' 'hieght(um)' 'neighboring_cell_count'];

    cell_area_summary = [header; num2cell(cell_areas)];
    writecell(cell_area_summary,fullfile(output_dir,'cell_areas.csv'),'Delimiter','comma');

    cell_plot = [cell_areas(:,3:3) cell_plot];
    center_frame = 1 + floor((p-1)/2);
    centre_plot_lines = [];
    for i=1:length(cell_plot)
        plot_line = zeros(1,p);
        z_centroid = cell_plot(i,1) + 1;
        plot_line(center_frame) = cell_plot(i,z_centroid);
        plot_index = z_centroid;
        for j=center_frame-1:-1:1
            plot_index = plot_index - 1;
            if plot_index>1
                plot_line(j) = cell_plot(i,plot_index);
            end
        end
        plot_index = z_centroid;
        for j=center_frame+1:p
            plot_index = plot_index + 1;
            if plot_index <= size(cell_plot,2)
                plot_line(j) = cell_plot(i,plot_index);
            end
        end
        centre_plot_lines = [centre_plot_lines; plot_line];
    end

    centre_plot_lines(centre_plot_lines==0) = NaN;
    figure;hold on;
    for i=1:length(centre_plot_lines)
       h(i)= plot(centre_plot_lines(i,:)');
       % Fourth last column in cell_areas is the cell track id : cell_areas(:,size(cell_areas,2)-3)
       y2 = cell_areas(i,size(cell_areas,2)-4).* ones(length(centre_plot_lines(i,:)),1);
       row = dataTipTextRow('Cell Track',y2);
       h(i).DataTipTemplate.DataTipRows(end+1) = row;
    end
    cmap = autumn(length(centre_plot_lines));
    for i=1:length(centre_plot_lines)
       set(h(i),'Color',cmap(i,:));
    end
    xlabel('Frame Number')
    ylabel('Cross-sectional Area')

    savefig(fullfile(output_dir,'cross_section_plot.fig'));
    save(fullfile(output_dir,'volume_view.mat'), 'cell_tracks', 'cell_areas', 'vol_interp_multi_cell', 'New_Mask_3D');
    
    % Update cell tracks file with cell track id
    cell_tracks = [cell_tracks cell_ids'];
    header = {};
    for j = 1:p
        header_txt_area = strcat('Frame_',num2str(j));
        header = [header header_txt_area];
    end
    header = [header 'cell_track_id'];
    cell_tracks_cells = [header; num2cell(cell_tracks)];
    writecell(cell_tracks_cells,fullfile(output_dir,'cell_tracks.csv'),'Delimiter','comma');
    
    %Writing cell neighborhood counts in each frame
    neighbor_cell_count_frame = [cell_ids' neighbor_cell_count_frame];
    header = {};
    for j = 1:p
        header_txt_area = strcat('Frame_',num2str(j));
        header = [header header_txt_area];
    end
    header = ['cell_track_id' header];
    neighbor_cell_count_frame = [header; num2cell(neighbor_cell_count_frame)];
    writecell(neighbor_cell_count_frame,fullfile(output_dir,'neighboring_cell_counts_frame.csv'),'Delimiter','comma');
    
    %Writing cell neighborhood tracks
    neighbor_cell_track_ids = [num2cell(cell_ids') neighbor_cell_track_ids];
    header = {'cell_track_id' 'neighborhood_track_id_list'};
    neighbor_cell_track_ids = [header; neighbor_cell_track_ids];
    writecell(neighbor_cell_track_ids,fullfile(output_dir,'neighboring_cell_tracks.csv'),'Delimiter','comma');
end

function contact_tracks = find_touching_cells(mask_slice_labels, frame_id, cell_id, cell_tracks, nhood_dialate)
    contact_tracks=[];
    
    % Get the slice image for this frame 
    I = mask_slice_labels(:,:,frame_id);
    % Extract the cell body
    lo = I == cell_id;
    % Create a dialate mask
    o_dialate = ones(nhood_dialate);
    % Dialate and extract the touching cells 
    touching_cell_ids = unique(I(imdilate(lo,o_dialate) + lo == 1));
    
    % Extract the track_ids
    this_frame_cells = cell_tracks(:,frame_id);
    for i=1:length(this_frame_cells)
        if (ismember(this_frame_cells(i), touching_cell_ids) && this_frame_cells(i)~=0)
            contact_tracks = [contact_tracks i];
        end
    end
end

function isValidCellId = cell_is_valid(cell_id, frame_id, cell_summary)

    isValidCellId = false;
    for k=2:length(cell_summary)
        if(cell_summary{k,1}==frame_id && cell_summary{k,2}==cell_id)
            isValidCellId = true;
            break;
        end
    end
end