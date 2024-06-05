slice0 = zeros(1024);
% slice1 = imread("data/1/handCorrection.tif");
% slice2 = imread("data/2/handCorrection.tif");
% slice3 = imread("data/3/handCorrection.tif");
% slice4 = imread("data/4/handCorrection.tif");
% slice5 = imread("data/5/handCorrection.tif");
% slice6 = imread("data/6/handCorrection.tif");
% slice7 = imread("data/7/handCorrection.tif");
% slice8 = imread("data/8/handCorrection.tif");
% slice9 = imread("data/9/handCorrection.tif");
% slice10 = imread("data/10/handCorrection.tif");
% slice11 = imread("data/11/handCorrection.tif");



% for i = 1:p;
%     Mask_3D_Test1(:,:,i) = imcomplement(Mask_3D_Test1(:,:,i));
%     s = regionprops(Mask_3D_Test1(:,:,i),'centroid');
%     centroids = cat(1,s.Centroid);
%     [cr cc] = size(centroids);
%     for j = 1:cr;
%         row_coord = uint16(centroids(j,1));
%         col_coord = uint16(centroids(j,2));
%         Mask_3D_Test1(col_coord,row_coord,i)=0;
%     end
% end
Mask_3D_Test1(:,:,7) = imcomplement(Mask_3D_Test1(:,:,7));
s = regionprops(Mask_3D_Test1(:,:,7),'centroid');
centroids = cat(1,s.Centroid);
[cr cc] = size(centroids);
for i = 1:p;
    if(i~=7)
        Mask_3D_Test1(:,:,i) = imcomplement(Mask_3D_Test1(:,:,i));
    end
    for j = 1:cr;
        row_coord = uint16(centroids(j,1));
        col_coord = uint16(centroids(j,2));
        Mask_3D_Test1(col_coord+1,row_coord,i)=0;
        Mask_3D_Test1(col_coord+1,row_coord+1,i)=0;
        Mask_3D_Test1(col_coord+1,row_coord-1,i)=0;
        Mask_3D_Test1(col_coord,row_coord,i)=0;
        Mask_3D_Test1(col_coord-1,row_coord+1,i)=0;
        Mask_3D_Test1(col_coord-1,row_coord-1,i)=0;
        Mask_3D_Test1(col_coord-1,row_coord,i)=0;
    end
end

s = sliceViewer(Mask_3D_Test1);


for i = 1:p;
    Mask_3D_Test1(:,:,i) = imcomplement(Mask_3D_Test1(:,:,i));
    %imshow(imfill(Mask_3D(:,:,i),[480 500]));
    %bw2 = imcomplement(Mask_3D(:,:,i));
    test_p1 = imfill(Mask_3D(:,:,i),[499 478]);
    test_p2 = imsubtract(test_p1,Mask_3D(:,:,i));
    Mask_3D_Test1(:,:,i) = test_p2;

    test_p1 = imfill(Mask_3D(:,:,i),[530 500]);
    test_p2 = imsubtract(test_p1,Mask_3D(:,:,i));
    Mask_3D_Test2(:,:,i) = test_p2;
    Mask_3D_Test3(:,:,i) = imfill(Mask_3D(:,:,i),[530 500]);
end

wt1 = linspace(0.4,0,10);
wt2 = linspace(0.3,0.4,10);
wt3 = linspace(0.2,0.3,10);
wt4 = linspace(0.1,0.2,10);

for i = 1:p-3;
    im1 = Mask_3D_Test1(:,:,i);
    im2 = Mask_3D_Test1(:,:,i+1);
    im3 = Mask_3D_Test1(:,:,i+2); 
    im4 = Mask_3D_Test1(:,:,i+3);
    vol_interp(:,:,(i*10)+1) = wt1(1)*im1 + wt2(1)*im2 + wt3(1)*im3 + wt4(1)*im4;
    vol_interp(:,:,(i*10)+2) = wt1(2)*im1 + wt2(2)*im2 + wt3(2)*im3 + wt4(2)*im4;
    vol_interp(:,:,(i*10)+3) = wt1(3)*im1 + wt2(3)*im2 + wt3(3)*im3 + wt4(3)*im4;
    vol_interp(:,:,(i*10)+4) = wt1(4)*im1 + wt2(4)*im2 + wt3(4)*im3 + wt4(4)*im4;
    vol_interp(:,:,(i*10)+5) = wt1(5)*im1 + wt2(5)*im2 + wt3(5)*im3 + wt4(5)*im4;
    vol_interp(:,:,(i*10)+6) = wt1(6)*im1 + wt2(6)*im2 + wt3(6)*im3 + wt4(6)*im4;
    vol_interp(:,:,(i*10)+7) = wt1(7)*im1 + wt2(7)*im2 + wt3(7)*im3 + wt4(7)*im4;
    vol_interp(:,:,(i*10)+8) = wt1(8)*im1 + wt2(8)*im2 + wt3(8)*im3 + wt4(8)*im4;
    vol_interp(:,:,(i*10)+9) = wt1(9)*im1 + wt2(9)*im2 + wt3(9)*im3 + wt4(9)*im4;
    vol_interp(:,:,(i*10)+10) = wt1(10)*im1 + wt2(10)*im2 + wt3(10)*im3 + wt4(10)*im4;
    
    % Add the centroid
    for k = 1:10
        for j = 1:cr;
            row_coord = uint16(centroids(j,1));
            col_coord = uint16(centroids(j,2));
            vol_interp(col_coord+1,row_coord,(i*10)+k)  =1;
            vol_interp(col_coord+1,row_coord+1,(i*10)+k)=1;
            vol_interp(col_coord+1,row_coord-1,(i*10)+k)=1;
            vol_interp(col_coord,row_coord,(i*10)+k)    =1;
            vol_interp(col_coord-1,row_coord+1,(i*10)+k)=1;
            vol_interp(col_coord-1,row_coord-1,(i*10)+k)=1;
            vol_interp(col_coord-1,row_coord,(i*10)+k)  =1;
        end
    end
end

for i = 1:p-3;
    for j = 1:10;
        %vol_interp(:,:,(i*10)+j) = vol_interp(:,:,(i*10)+j) .* Mask_3D(:,:,i);
        %vol_interp(:,:,(i*10)+j) = im2bw(imoverlay(vol_interp(:,:,(i*10)+j),Mask_3D(:,:,i),'red'),0.4);
        for inr = 1:r;
            for inc = 1:c;
                if Mask_3D(inr,inc,i) == true
                    vol_interp(inr,inc,(i*10)+j)=1;
                end
            end
        end
    end
end

s = sliceViewer(vol_interp);


vol_interp(:,:,1) =  Mask_3D_Test1(:,:,1);
vol_interp(:,:,1) = vol_interp(:,:,1) | Mask_3D(:,:,1);
for i = 2:p-1;
    im1 = Mask_3D_Test1(:,:,i-1);
    im2 = Mask_3D_Test1(:,:,i);
    im3 = Mask_3D_Test1(:,:,i+1);
    end_indx = ((i-1)*10)+1;
    vol_interp(:,:,end_indx-9) = logical(round(0.6*im1 + 0.3*im2 + 0.1*im3));
    vol_interp(:,:,end_indx-9) = vol_interp(:,:,end_indx-9) | Mask_3D(:,:,i-1);
    
    vol_interp(:,:,end_indx-8) = 0.5375*im1 + 0.3375*im2 + 0.125*im3;
    vol_interp(:,:,end_indx-8) = vol_interp(:,:,end_indx-8) | Mask_3D(:,:,i-1);
    
    vol_interp(:,:,end_indx-7) = 0.475*im1 + 0.375*im2 + 0.15*im3;
    vol_interp(:,:,end_indx-7) = vol_interp(:,:,end_indx-7) | Mask_3D(:,:,i-1);
        
    vol_interp(:,:,end_indx-6) = 0.4125*im1 + 0.4125*im2 + 0.175*im3;
    vol_interp(:,:,end_indx-6) = vol_interp(:,:,end_indx-6) | Mask_3D(:,:,i-1);
        
    vol_interp(:,:,end_indx-5) = 0.35*im1 + 0.45*im2 + 0.2*im3;
    vol_interp(:,:,end_indx-5) = vol_interp(:,:,end_indx-5) | Mask_3D(:,:,i-1);
    
    vol_interp(:,:,end_indx-4) = 0.2875*im1 + 0.4875*im2 + 0.225*im3;
    vol_interp(:,:,end_indx-4) = vol_interp(:,:,end_indx-4) | Mask_3D(:,:,i-1);

    vol_interp(:,:,end_indx-3) = 0.225*im1 + 0.525*im2 + 0.25*im3;
    vol_interp(:,:,end_indx-3) = vol_interp(:,:,end_indx-3) | Mask_3D(:,:,i-1);
        
    vol_interp(:,:,end_indx-2) = 0.1625*im1 + 0.5625*im2 + 0.275*im3;
    vol_interp(:,:,end_indx-2) = vol_interp(:,:,end_indx-2) | Mask_3D(:,:,i-1);
        
    vol_interp(:,:,end_indx-1) = 0.1*im1 + 0.6*im2 + 0.3*im3;
    vol_interp(:,:,end_indx-1) = vol_interp(:,:,end_indx-1) | Mask_3D(:,:,i-1);
        
    vol_interp(:,:,end_indx) = im2;
    vol_interp(:,:,end_indx) = vol_interp(:,:,end_indx) | Mask_3D(:,:,i-1);
end
%vol_interp = round(vol_interp);
s = sliceViewer(vol_interp);
for inr = 2:r;
  for jnr = 2:c;
	  xi = Mask_3D(inr,jnr,:);
	  x = linspace(1,12,12);
	  xq = linspace(1,12,120);
	  xi = reshape(xi,[1,12]);
	  sumxi = sum(xi);
	  if (sumxi >= 1)
		vol_interp(inr,jnr,:) = ones(1,120)*10;
		vol_interp(inr-1,jnr,:) = ones(1,120)*10;
		vol_interp(inr,jnr-1,:) = ones(1,120)*10;
		vol_interp(inr-1,jnr-1,:) = ones(1,120)*10;
	  else
		vol_interp(inr,jnr,:) = zeros(1,120);
	  end
  end;
end;

vol_interp = logical(vol_interp);
for i = 1:120;
    vol_interp_1(:,:,i) = imfill(vol_interp(:,:,i),[480 500]);
    vol_interp_2(:,:,i) = imfill(vol_interp(:,:,i),[530 500]);
    vol_interp_3(:,:,i) = imfill(vol_interp(:,:,i),[530 500]);
end

%intensity = [-3024,-16.45,641.38,3071];
%alpha = [0, 0, 0.72, 0.72];
%color = ([0 0 0; 186 65 77; 231 208 141; 255 255 255]) ./ 255;

intensity = [0 20 40 120 220 1024];
alpha = [0 0 0.15 0.3 0.38 0.5];
color = ([0 0 0; 43 0 0; 103 37 20; 199 155 97; 216 213 201; 255 255 255]) ./ 255;
queryPoints = linspace(min(intensity),max(intensity),256);
alphamap = interp1(intensity,alpha,queryPoints)';
colormap = interp1(intensity,color,queryPoints);
figure
vol = volshow(vol_interp,'Colormap',colormap,'Alphamap',alphamap);
[L,NumLabels] = superpixels3(vol_interp,1000);
imSize = size(vol_interp);

% 
% intensity = [0,1];
% alpha = [0, 0.72];
% color = ([0 0 0; 255 255 255]) ./ 255;
% queryPoints = linspace(min(intensity),max(intensity),256);
% alphamap = interp1(intensity,alpha,queryPoints)';
% colormap = interp1(intensity,color,queryPoints);
% figure
% vol = volshow(vol_interp,'Colormap',colormap,'Alphamap',alphamap);
% 
% imPlusBoundaries = zeros(imSize(1),imSize(2),3,imSize(3),'uint8');
% for plane = 1:imSize(3)
%   BW = boundarymask(L(:, :, plane));
%   % Create an RGB representation of this plane with boundary shown
%   % in cyan.
%   imPlusBoundaries(:, :, :, plane) = imoverlay(vol_interp(:, :, plane), BW, 'cyan');
% end
% implay(imPlusBoundaries,5)
% 
% 
% BW = L(:, :, 10);
% BW(BW>200)=0;
% BW(BW<200)=0;
% figure
% imshow(logical(vol2dp));
% 
% 
% for cell = 1:256
% 	vol_interp_cell = zeros(r,c,10*p);
% 	BW = L(:, :, 10);
% 	BW(BW>cell)=0;
% 	BW(BW<cell)=0;
% 	for plane = 1:imSize(3)
% 		vol2dp = vol_interp(:, :, plane);
% 		vol2dp(BW==0) = 0;
% 		vol2dp(BW>0) = 1;
% 		vol_interp_cell(:, :, plane) = vol2dp;
% 	end
% 	vol_1 = volshow(vol_interp_1,'Colormap',colormap,'Alphamap',alphamap);
% 	f = fullfile('C:\\Users\\ab50\\Documents\\git\\unet\\dist\\test_batch\\',cell,'.')
% 	imwrite(vol_interp_cell,"")
% end
% 
% 
% 


cell=78
vol_interp_cell = zeros(r,c,10*p);
vol_interp_bck = zeros(r,c,10*p);
BW = L(:, :, 10);
BW(BW>cell)=0;
BW(BW<cell)=0;
for plane = 1:imSize(3)
    vol2dp = vol_interp(:, :, plane);
    vol2dp(BW==0) = 0;
    vol2dp(BW>0) = 750;
    vol_interp_cell(:, :, plane) = vol2dp;

    vol2dp = vol_interp(:, :, plane);
    vol2dp(BW>0) = 750;
    vol_interp_bck(:, :, plane) = vol2dp;
end
figure
%vol_cell = volshow(vol_interp_cell,'Colormap',colormap);
%vol_whole = volshow(vol_interp_bck,'Colormap',colormap);
vol_interp_bck = uint8(vol_interp_bck);
cmap = parula(256);
s = sliceViewer(vol_interp,'Colormap',cmap);













% 
% 
% 
% 
% 
% 
% 
% cell_id_list= [137 150 171 182 164 154 147 168];%[168 147 138 154 164];
% cell_count = length(cell_id_list);
% 
% vol_interp_multi_cell = zeros(r,c,intrep_frame_count*p);
% 
% for cell_id = 1:cell_count;
%     vol_interp = zeros(r,c,intrep_frame_count*p);
%     Mask_3D_Test1 = Mask_3D;
%     for i = 1:p;
%         e = s(cell_id_list(cell_id)).Centroid;
%         test_p1 = imfill(Mask_3D(:,:,i), ceil([e(1,2) e(1,1)]));
%         test_p2 = imsubtract(test_p1,Mask_3D(:,:,i));
%         Mask_3D_Test1(:,:,i) = test_p2;
%     end
% 
%     cell_px_area = sum(Mask_3D_Test1(:,:,7),'all');
%     cell_terminated=0;
%     prev_plane_area = sum(Mask_3D_Test1(:,:,1),'all');
%     for i = 1:p-1;
%         if(cell_terminated==1)
%             break;
%         end
%         curr_plane = Mask_3D_Test1(:,:,i);
%         curr_plane_area = sum(curr_plane,'all');
%         if(curr_plane_area > cell_px_area*3.0 || curr_plane_area<100)
%             for j = 1:intrep_frame_count;
%                 vol_interp(:,:,(i-1)*intrep_frame_count+j) = 0;
%             end
%             prev_plane_area = sum(Mask_3D_Test1(:,:,i),'all');
%         else
%             %assign current plane
%             vol_interp(:,:,(i-1)*intrep_frame_count+1) = curr_plane;
%             % zero added planes
%             for j = 2:intrep_frame_count;
%                 vol_interp(:,:,(i-1)*intrep_frame_count+j) = 0;
%             end
%             % interpolate to next plane
%             next_plane=Mask_3D_Test1(:,:,i+1);
%             next_plane_area = sum(next_plane,'all');
%             if(next_plane_area > cell_px_area*3 || next_plane_area<100 || next_plane_area > 3*prev_plane_area) % No cell identified in next plane
%                 max_row=1;
%                 min_row=r;
%                 % Find max min row range
%                 for inr = 2:r;
%                     found_indx = find(curr_plane(inr,:)==1);
%                     if(numel(found_indx)>0)
%                         if(inr>max_row)
%                             max_row=inr;
%                         end
%                         if(inr<min_row)
%                             min_row = inr;
%                         end
%                     end
%                 end
%                 mid_row = min_row + (max_row-min_row)/2;
%                 start_row = round(linspace(min_row,mid_row,9));
%                 start_row = start_row(1:end-uint8(length(start_row)/2));
%                 end_row = round(linspace(max_row,mid_row,9));
%                 end_row = end_row(1:end-uint8(length(end_row)/2));
%                 for j = 2:uint8(intrep_frame_count/2);
%                     for inr = start_row(j-1):end_row(j-1);
%                         found_indx = find(curr_plane(inr,:)==1);
%                         if(numel(found_indx)>0)
%                             max_pos = found_indx(numel(found_indx));
%                             min_pos = found_indx(1);
%                             mid_pos = min_pos + (max_pos-min_pos)/2;
%                             start_pos = round(linspace(min_pos,mid_pos,9));
%                             end_pos = round(linspace(max_pos,mid_pos,9));
%                              for inc = start_pos(j-1):end_pos(j-1);
%                                 vol_interp(inr,inc,(i-1)*intrep_frame_count+j) = 1;
%                              end
%                         end
%                     end
%                 end
%                 cell_terminated=1;
%             else
%                 for inr = 2:r;
%                     found_curr_indx = find(curr_plane(inr,:)==1);
%                     found_next_indx = find(next_plane(inr,:)==1);
%                     if(numel(found_curr_indx)>0 && numel(found_next_indx)==0) % cell not yet found in next plane
%                         max_pos = found_curr_indx(end);
%                         min_pos = found_curr_indx(1);
%                         mid_pos = min_pos + (max_pos-min_pos)/2;
%                         start_pos = round(linspace(min_pos,mid_pos,9));
%                         end_pos = round(linspace(max_pos,mid_pos,9));
%                         for j = 2:intrep_frame_count;
%                              for inc = start_pos(j-1):end_pos(j-1);
%                                 vol_interp(inr,inc,(i-1)*intrep_frame_count+j) = 1;
%                              end
%                         end
%                     elseif(numel(found_curr_indx)>0 && numel(found_next_indx)>0)
%                         start_pos = round(linspace(found_curr_indx(1),found_next_indx(1),9));
%                         end_pos = round(linspace(found_curr_indx(end),found_next_indx(end),9));
%                         for j = 2:intrep_frame_count;
%                              for inc = start_pos(j-1):end_pos(j-1);
%                                 vol_interp(inr,inc,(i-1)*intrep_frame_count+j) = 1;
%                              end
%                         end
%                     end
%                 end
%             end
%             prev_plane_area = next_plane_area;
%         end
%     end
%     for i = 1:intrep_frame_count*p;
%         for inr = 1:r;
%             for inc = 1:c;
%                 if vol_interp(inr,inc,i) == true
%                     vol_interp_multi_cell(inr,inc,i) = cell_id + 1;
%                 end
%             end
%         end
%     end
% end
% 
% % Display the 3D volume
% cmap = jet(cell_count+1);
% figure; hold on;
% for cell_id = 1:cell_count;
% 	A = vol_interp_multi_cell;
% 	A(A~=cell_id+1)=0;
% 	h = patch(isosurface(A,1));
% 	h.EdgeColor = 'none';
% 	h.FaceColor = cmap(cell_id+1,:);
% 	h.FaceAlpha = 0.7;
% end
% view(3); 
% axis tight
% camlight 
% lighting gouraud
% 
% %Display the slice view of mask
% vol_interp_mask = vol_interp_multi_cell;
% %Add the MASK
% for i = 1:p-1;
%     for j = 1:intrep_frame_count;
%         for inr = 1:r;
%             for inc = 1:c;
%                 if Mask_3D(inr,inc,i) == true
%                     vol_interp_mask(inr,inc,(i-1)*intrep_frame_count+j)=1;
%                 end
%             end
%         end
%     end
% end
% figure;
% sliceViewer(vol_interp_mask,'Colormap',cmap);
% 
% % Display the cell on the images
% slice_orig_1 = imread("orig_image/1.tif");
% slice_orig_2 = imread("orig_image/2.tif");
% slice_orig_3 = imread("orig_image/3.tif");
% slice_orig_4 = imread("orig_image/4.tif");
% slice_orig_5 = imread("orig_image/5.tif");
% slice_orig_6 = imread("orig_image/6.tif");
% slice_orig_7 = imread("orig_image/7.tif");
% slice_orig_8 = imread("orig_image/8.tif");
% slice_orig_9 = imread("orig_image/9.tif");
% slice_orig_10 = imread("orig_image/10.tif");
% slice_orig_11 = imread("orig_image/11.tif");
% 
% Orig_3D(:,:,1) = slice_orig_1;
% Orig_3D(:,:,2) = slice_orig_2;
% Orig_3D(:,:,3) = slice_orig_3;
% Orig_3D(:,:,4) = slice_orig_4;
% Orig_3D(:,:,5) = slice_orig_5;
% Orig_3D(:,:,6) = slice_orig_6;
% Orig_3D(:,:,7) = slice_orig_7;
% Orig_3D(:,:,8) = slice_orig_8;
% Orig_3D(:,:,9) = slice_orig_9;
% Orig_3D(:,:,10) = slice_orig_10;
% Orig_3D(:,:,11) = slice_orig_11;
% 
% [ro co po] = size(Orig_3D);
% vol_orig_intrep = zeros(r,c,intrep_frame_count*p);
% 
% for i = 1:po-1;
%     for j = 1:intrep_frame_count;
%         vol_orig_intrep(:,:,(i-1)*intrep_frame_count + j) = Orig_3D(:,:,i);
%     end
% end
% 
% for i = 1:po-1;
%     for j = 1:intrep_frame_count;
%         for inr = 1:ro;
%             for inc = 1:co;
%                 if vol_interp_mask(inr,inc,(i-1)*intrep_frame_count+j) > 0
%                     vol_orig_intrep(inr,inc,(i-1)*intrep_frame_count+j)=255;
%                 end
%             end
%         end
%     end
% end
% figure;
% cmap = jet(cell_count+1);
% sliceViewer(vol_orig_intrep,'Colormap',cmap);
% 
% 


