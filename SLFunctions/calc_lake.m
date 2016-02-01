function lake = calc_lake(ice_j_corr,oc_j,topo_j,lat_out,lon_out)

% Find outline of ice_j and generate a matrix that is zero when the pixel
% is not at the boundary and 1 when it is.
B = bwboundaries(sign(ice_j_corr));
boundary_mat = zeros(size(oc_j));
for i = 1:length(B)
    for j = 1:length(B{i})
        boundary_mat(B{i}(j,1),B{i}(j,2)) = 1;
    end
end

% POTENTIALLY CALCULATE THIS PART ON HIGH RES!!
% Calculate all depressions on continents (only look at North
% America for now)
topo_cont = topo_j .* (ones(size(oc_j))-oc_j);
topo_cont(lat_out(:,1) < 30,:) = 0;
topo_cont(:,lon_out(1,:) < 200) = 0;
topo_filled = (imfill(topo_cont)-topo_cont) .* ...
    (ones(size(oc_j))-sign(ice_j_corr));

% Find depressions that are adjacent to ice sheets
ice_lake_overlap = sign(topo_filled) + sign(boundary_mat);
CC = bwconncomp(ice_lake_overlap);

% group patches by connectivity
lake = zeros(size(oc_j));
for j = 1:CC.NumObjects
    [r,c] = ind2sub(size(oc_j),CC.PixelIdxList{j});
    mat = zeros(size(oc_j));
    for i = 1:length(r)
        mat(r(i),c(i)) = 1;
    end
    % check that patch is at the ice margin but not just the ice
    % margin
    if max(max(mat + boundary_mat)) == 2 && ...
            max(max(mat-boundary_mat))==1;
        lake = lake + mat;
    end
end

% take ice boundary out of patch
lake = (lake - boundary_mat);
lake(lake<0) = 0;
% multiply lake patch by its depth
lake = lake .* topo_filled;