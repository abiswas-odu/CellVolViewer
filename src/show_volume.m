function show_volume(config_dir, reverse, rotate_view, gray_image)
    config_file = fullfile(config_dir,'config.txt');
    % Specify tracks to display empty means ALL
    inc_cell_tracks=[];
    
    %Check if reverse view
    if ~exist('reverse','var')
        reverse='false';
    end
    
    %Check if rotate view
    if ~exist('rotate_view','var')
        rotate_view = 'false';
    end
    
    %Check if gray image requested
    if ~exist('gray_image','var')
        gray_image = 'false';
    end
    
    % Load configuration 
    [opt, ~, ~, symbolicDefs] = readConfig(config_file);
    show_vol_axis = 1;
    if isfield(opt,'show_axis')
        show_vol_axis = opt.show_axis;
    end
    z_axis_scale = opt.x_px;
    if isfield(opt,'z_scale')
        z_axis_scale = opt.z_scale;
    end
    
    output_dir = fullfile(config_dir,'output');
    load(fullfile(output_dir,'volume_view.mat'));
    fprintf('Displaying cell volume...\n');
    bg = [1 1 1; 0 0 0];
    cmap = zeros(size(cell_tracks,1),3);
    % Load color file if available 
    if isfile(fullfile(output_dir,'cell_colors.csv'))
        cell_track_colors = readtable(fullfile(output_dir,'cell_colors.csv'));
        for i=1:size(cell_track_colors,1)
            track_id = cell_track_colors{i,1};
            track_color = rgb(char(cell_track_colors{i,2}));
            cmap(track_id,:) = track_color;
        end
    elseif strcmp(gray_image,'gray')
        % Setup gray random coloring
        % Fourth last column in cell_areas is the cell track id : cell_areas(:,size(cell_areas,2)-3)
        color_set = flipud(gray(size(cell_tracks,1)));
        rng(1023,'twister');
        dataset = 1:1:size(cell_tracks,1);
        r = randsample(dataset, size(cell_tracks,1), false);
        for i=1:size(cell_tracks,1)
            cmap(i,:) = color_set(r(i),:);
        end
    else
        % Setup rainbow coloring based of distance from center 
        % Fourth last column in cell_areas is the cell track id : cell_areas(:,size(cell_areas,2)-3)
        color_set = flipud(jet(size(cell_tracks,1)));
        for i=1:size(cell_tracks,1)
            track_dist_pos = find(cell_areas(:,size(cell_areas,2)-4)==i);
            track_color = color_set(track_dist_pos,:);
            cmap(i,:) = track_color;
        end
        % cmap = distinguishable_colors(size(cell_tracks,1),bg);
    end
    cmap = [0 0 0; cmap];
    
    intrep_frame_count = opt.inprep_frames;
    p = uint8(size(vol_interp_multi_cell,3)/intrep_frame_count);
    for j = 1:p
        tiff = vol_interp_multi_cell(:,:,(j-1)*intrep_frame_count+1);
        outputFileName = fullfile(output_dir,strcat(num2str(j,'%02.f'),'.tif'));
        imwrite(uint8(tiff),outputFileName,'WriteMode', 'overwrite');
    end
    
    % Display the 3D volume
    figure; hold on;
    for i=1:size(cell_tracks,1)
        is_in_view = cell_in_view(opt,cell_areas,i);
        if((isempty(inc_cell_tracks) || ismember(i, inc_cell_tracks)) && is_in_view)
            A = vol_interp_multi_cell;
            %Reverse the view to the required orientation 
            if strcmp(reverse,'true')
                A = flip(A,3);
            end
            A(A~=i+1)=0;
            h = patch(isosurface(A,1));
            h.EdgeColor = 'none';
            h.FaceColor = cmap(i+1,:);
            h.FaceAlpha = 0.7;
        end
    end
    %Rotate the view to the required orientation 
    if strcmp(rotate_view,'true')
        view([112.5 30]); 
    else
        view(3); 
    end
    if show_vol_axis == 1
        axis tight
    else
        axis off
    end
    set(gca,'YDir','reverse');
    camlight 
    lighting gouraud
    zlim([0 z_axis_scale]);
    savefig(fullfile(output_dir,strcat(datestr(now,'mm-dd-yyyy-HH-MM'),'_3D_volume_plot.fig')));
    fprintf('Cell volume display complete.\n');
end

function isIncluded = cell_in_view(opt,cell_areas,cell_track_id)

    isIncluded = true;
    if isfield(opt,'x_view_min_px') && isfield(opt,'x_view_max_px') && isfield(opt,'y_view_min_px') && isfield(opt,'y_view_max_px')
        x_min = opt.x_view_min_px * (opt.x_real/opt.x_px);
        x_max = opt.y_view_max_px * (opt.x_real/opt.x_px);
        y_min = opt.x_view_min_px * (opt.x_real/opt.x_px);
        y_max = opt.y_view_max_px * (opt.x_real/opt.x_px);
        % Fourth last column in cell_areas is the cell track id: cell_areas(:,size(cell_areas,2)-3)
        track_pos = find(cell_areas(:,size(cell_areas,2)-4)==cell_track_id);
        x_cent = cell_areas(track_pos,1);
        y_cent = cell_areas(track_pos,2);
        if(x_cent>=x_min && x_cent<=x_max && y_cent>=y_min && y_cent<=y_max)
            isIncluded = true;
        else
            isIncluded = false;
        end
    end
end