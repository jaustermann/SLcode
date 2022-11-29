function [lake_LR, lake] = calc_lake_HR(ice_HR,del_topo,lon_out, lat_out,LON_HR,LAT_HR,topo_HR,sampling_factor)

%ice_HR = interp2(lon_out,lat_out,ice_j_corr,LON_HR,LAT_HR);
del_topo_HR = interp2(lon_out,lat_out,del_topo,LON_HR,LAT_HR); 
topo_j_HR = topo_HR + del_topo_HR + ice_HR; %include ice in topography
below_0_HR = sign_01(topo_j_HR);

% do an adjusted ocean function that only removes the deep ocean
CC = bwconncomp(below_0_HR,8);
pixel_len = zeros(CC.NumObjects,1);
for j = 1:CC.NumObjects
    pixel_len(j) = length(CC.PixelIdxList{j});
end
% find the connected areas that are massive (Atlantic and Pacific)
ind = find(pixel_len > 1e6);
% group big patches
oc_HR = zeros(size(below_0_HR));
for i = 1:length(ind)
    mat = zeros(size(below_0_HR));
    mat(CC.PixelIdxList{ind(i)}) = 1;
    oc_HR = oc_HR+mat;
end


% Find outline of ice_j and generate a matrix that is zero when the pixel
% is not at the boundary and 1 when it is.
temp = sign(ice_HR-50);
temp(temp<0) = 0;
B = bwboundaries(temp);
boundary_mat = zeros(size(ice_HR));
for i = 1:length(B)
    for j = 1:length(B{i})
        boundary_mat(B{i}(j,1),B{i}(j,2)) = 1;
    end
end


% Calculate all depressions on continents
topo_cont = topo_j_HR .* (ones(size(oc_HR))-oc_HR); % limit to on land
% set values outside of the prescribed range to 0
% topo_cont(LAT_HR(:,1) < lat_limit(1) | LAT_HR(:,1) > lat_limit(2),:) = 0;
% topo_cont(:,LON_HR(1,:) < lon_limit(1) | LON_HR(1,:) > lon_limit(2)) = 0;
% 

topo_cont(isnan(topo_cont) ==1) = 0;
topo_filled = (imfill(topo_cont)-topo_cont) .* ... 
    (ones(size(oc_HR))-temp); % remove depressions on ice sheet

% Find depressions that are adjacent to ice sheets
% find patches that connected to each other
ice_lake_overlap = sign(topo_filled) + sign(boundary_mat);
CC = bwconncomp(ice_lake_overlap,8); % 4 means in 4 neighboring directions, other option is 8 (for diagonal connectivity)

% remove lakes less than 100/sampling_factor pixels because those aren't connected to the
% ice sheet anyways (because the pixels of the entire ice boundary add to the
% pixel count), unless the ice margin + lake have fewer than those pixels,
% if so assume that it's small enough to neglect. 
pixel_len = zeros(CC.NumObjects,1);
for j = 1:CC.NumObjects
    pixel_len(j) = length(CC.PixelIdxList{j});
end
ind = find(pixel_len>100/sampling_factor);

% group patches by connectivity
lake = zeros(size(oc_HR));
for i = 1:length(ind)
    mat = zeros(size(oc_HR));
    mat(CC.PixelIdxList{ind(i)}) = 1;
    
%     [r,c] = ind2sub(size(oc_HR),CC.PixelIdxList{j});
%     mat = zeros(size(oc_HR));
%     for i = 1:length(r)
%         mat(r(i),c(i)) = 1;
%     end
    % check that patch is at the ice margin but not just the ice
    % margin
    if max(max(mat + boundary_mat)) == 2 && ... % it contains ice boundary elements
            max(max(mat - boundary_mat))==1 % it still has elements if the boundary is removed
        
        % remove the ice margin part
        mat(mat + boundary_mat == 2) = 0;
        
        lake = lake + mat;
    end
end

% take ice boundary out of patch
lake = (lake - boundary_mat);
lake(lake<0) = 0; 
% multiply lake patch by its depth
lake = lake .* topo_filled;

lake_LR = interp2(LON_HR,LAT_HR,lake,lon_out,lat_out);
lake_LR(isnan(lake_LR) == 1) = 0;

%%

% [lat,lon] = gcwaypts(54,263-360,46,263-360,30);
% lon = lon+360;
% 
% clear ice_profile RSL_profile topo_profile lake_profile
% for i = 1:length(lat)
%     ice_profile(i) = interp2(LON_HR,LAT_HR,ice_HR,lon(i),lat(i));
%     
%     topo_profile(i) = interp2(LON_HR,LAT_HR,topo_cont-ice_HR,lon(i),lat(i));
%     
%     lake_profile(i) = interp2(LON_HR,LAT_HR,lake,lon(i),lat(i));
%     
% end
% 
% figure
% hold on
% plot(lat,lake_profile+topo_profile+ice_profile,'b')
% plot(lat,ice_profile+topo_profile,'c')
% plot(lat,topo_profile,'k')
% xlim([min(lat) max(lat)])
% xlabel('latitude')
% ylabel('elevation (m)')