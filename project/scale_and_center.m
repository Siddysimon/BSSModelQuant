function [segmentationI, segmentationF, segmentationD, segmentationDDeletion, segmentationDAddition] = scale_and_center(headI, headF, headDiff, geometry, electCen, radNew, radOld)

%convert given mask array to collection of xyz coordinates within a certain
%radius of the electrode center
tissues = {'wm','gm','csf','bone','skin','eyes','air'};

%Calculating for the final version of skin mask, to use as scaling
%reference
skinMask = headF.skin;
[x,y,z] = ndgrid(1:size(skinMask, 1), 1:size(skinMask, 2), 1:size(skinMask, 3)); %get coordinates of all points, could also use meshgrid
skinReference(:,1) = x(skinMask == 1);
skinReference(:,2) = y(skinMask == 1); %these 3 lines also reshape the 3d array
skinReference(:,3) = z(skinMask == 1); %into column vectors

head_vox_len_x = max(skinReference(:,1)) - min(skinReference(:,1));
head_vox_len_y = max(skinReference(:,2)) - min(skinReference(:,2));
head_vox_len_z = max(skinReference(:,3)) - min(skinReference(:,3));

% 1.) Find the voxel range for the given head, using skin mask.
% 2.) Use this length to find maskChange factor. 
vox_x_weight = (max(geometry.node(:,1)) - min(geometry.node(:,1))) ./ head_vox_len_x;
vox_y_weight = (max(geometry.node(:,2)) - min(geometry.node(:,2))) ./ head_vox_len_y;
vox_z_weight = (max(geometry.node(:,3)) - min(geometry.node(:,3))) ./ head_vox_len_z;

skinReference(:,1) = skinReference(:,1) .* vox_x_weight;
skinReference(:,2) = skinReference(:,2) .* vox_y_weight;
skinReference(:,3) = skinReference(:,3) .* vox_z_weight;

skinReference(:, 1) = -skinReference(:, 1);

correctMax = max(skinReference) - max(geometry.node);
skinReference = skinReference - correctMax;
correctMin = min(skinReference) - min(geometry.node);
skinReference = skinReference - correctMin;

%Loop through all tissue
for t = 1:length(tissues)
    maskI = headI.(tissues{t});
    maskF = headF.(tissues{t});
    maskChange = headDiff.(tissues{t});
    maskDeletion = (maskChange == -1);
    maskAddition = (maskChange == 1);
    maskChange = abs(maskChange);

    
    [x,y,z] = ndgrid(1:size(maskI, 1), 1:size(maskI, 2), 1:size(maskI, 3)); %get coordinates of all points, could also use meshgrid
    segmentationI.(tissues{t})(:,1) = x(maskI == 1);
    segmentationI.(tissues{t})(:,2) = y(maskI == 1); %these 3 lines also reshape the 3d array
    segmentationI.(tissues{t})(:,3) = z(maskI == 1); %into column vectors


    [x,y,z] = ndgrid(1:size(maskF, 1), 1:size(maskF, 2), 1:size(maskF, 3)); %get coordinates of all points, could also use meshgrid
    segmentationF.(tissues{t})(:,1) = x(maskF == 1);
    segmentationF.(tissues{t})(:,2) = y(maskF == 1); %these 3 lines also reshape the 3d array
    segmentationF.(tissues{t})(:,3) = z(maskF == 1); %into column vectors


    [x,y,z] = ndgrid(1:size(maskChange, 1), 1:size(maskChange, 2), 1:size(maskChange, 3)); %get coordinates of all points, could also use meshgrid
    segmentationD.(tissues{t})(:,1) = x(maskChange == 1);
    segmentationD.(tissues{t})(:,2) = y(maskChange == 1); %these 3 lines also reshape the 3d array
    segmentationD.(tissues{t})(:,3) = z(maskChange == 1); %into column vectors
    
    %Segmentation Deletions and Additions separated are calculated too:
    [x,y,z] = ndgrid(1:size(maskDeletion, 1), 1:size(maskDeletion, 2), 1:size(maskDeletion, 3)); %get coordinates of all points, could also use meshgrid
    segmentationDDeletion.(tissues{t})(:,1) = x(maskDeletion == 1);
    segmentationDDeletion.(tissues{t})(:,2) = y(maskDeletion == 1); %these 3 lines also reshape the 3d array
    segmentationDDeletion.(tissues{t})(:,3) = z(maskDeletion == 1); %into column vectors
    
    [x,y,z] = ndgrid(1:size(maskAddition, 1), 1:size(maskAddition, 2), 1:size(maskAddition, 3)); %get coordinates of all points, could also use meshgrid
    segmentationDAddition.(tissues{t})(:,1) = x(maskAddition == 1);
    segmentationDAddition.(tissues{t})(:,2) = y(maskAddition == 1); %these 3 lines also reshape the 3d array
    segmentationDAddition.(tissues{t})(:,3) = z(maskAddition == 1); %into column vectors

    segmentationD.(tissues{t})(:,1) = segmentationD.(tissues{t})(:,1) .* vox_x_weight;
    segmentationD.(tissues{t})(:,2) = segmentationD.(tissues{t})(:,2) .* vox_y_weight;
    segmentationD.(tissues{t})(:,3) = segmentationD.(tissues{t})(:,3) .* vox_z_weight;
    
    
    segmentationDDeletion.(tissues{t})(:,1) = segmentationDDeletion.(tissues{t})(:,1) .* vox_x_weight;
    segmentationDDeletion.(tissues{t})(:,2) = segmentationDDeletion.(tissues{t})(:,2) .* vox_y_weight;
    segmentationDDeletion.(tissues{t})(:,3) = segmentationDDeletion.(tissues{t})(:,3) .* vox_z_weight;
    
    
    segmentationDAddition.(tissues{t})(:,1) = segmentationDAddition.(tissues{t})(:,1) .* vox_x_weight;
    segmentationDAddition.(tissues{t})(:,2) = segmentationDAddition.(tissues{t})(:,2) .* vox_y_weight;
    segmentationDAddition.(tissues{t})(:,3) = segmentationDAddition.(tissues{t})(:,3) .* vox_z_weight;

    segmentationI.(tissues{t})(:,1) = segmentationI.(tissues{t})(:,1) .* vox_x_weight;
    segmentationI.(tissues{t})(:,2) = segmentationI.(tissues{t})(:,2) .* vox_y_weight;
    segmentationI.(tissues{t})(:,3) = segmentationI.(tissues{t})(:,3) .* vox_z_weight;

    segmentationF.(tissues{t})(:,1) = segmentationF.(tissues{t})(:,1) .* vox_x_weight;
    segmentationF.(tissues{t})(:,2) = segmentationF.(tissues{t})(:,2) .* vox_y_weight;
    segmentationF.(tissues{t})(:,3) = segmentationF.(tissues{t})(:,3) .* vox_z_weight;

    segmentationI.(tissues{t})(:, 1) = -segmentationI.(tissues{t})(:, 1);
    segmentationF.(tissues{t})(:, 1) = -segmentationF.(tissues{t})(:, 1);
    segmentationD.(tissues{t})(:, 1) = -segmentationD.(tissues{t})(:, 1);
    segmentationDDeletion.(tissues{t})(:, 1) = -segmentationDDeletion.(tissues{t})(:, 1);
    segmentationDAddition.(tissues{t})(:, 1) = -segmentationDAddition.(tissues{t})(:, 1);
    
    segmentationI.(tissues{t}) = segmentationI.(tissues{t}) - correctMax;
    segmentationF.(tissues{t}) = segmentationF.(tissues{t}) - correctMax;
    segmentationD.(tissues{t}) = segmentationD.(tissues{t}) - correctMax;
    segmentationDDeletion.(tissues{t}) = segmentationDDeletion.(tissues{t}) - correctMax;
    segmentationDAddition.(tissues{t}) = segmentationDAddition.(tissues{t}) - correctMax;

    segmentationI.(tissues{t}) = segmentationI.(tissues{t}) - correctMin;
    segmentationF.(tissues{t}) = segmentationF.(tissues{t}) - correctMin;
    segmentationD.(tissues{t}) = segmentationD.(tissues{t}) - correctMin;
    segmentationDDeletion.(tissues{t}) = segmentationDDeletion.(tissues{t}) - correctMin;
    segmentationDAddition.(tissues{t}) = segmentationDAddition.(tissues{t}) - correctMin;

    %Calculating ROI of converted data
    ROI_centers_rep = repmat(electCen, length(segmentationD.(tissues{t})(:,1)), 1);
    flagROISphere = vecnorm(segmentationD.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew;
    flagROIShell = (radOld <= vecnorm(segmentationD.(tissues{t}) - ROI_centers_rep, 2, 2)) & (vecnorm(segmentationD.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew);
    segmentationD.shell.(tissues{t}) = segmentationD.(tissues{t})(flagROIShell, :);
    segmentationD.(tissues{t}) = segmentationD.(tissues{t})(flagROISphere, :);

    
    ROI_centers_rep = repmat(electCen, length(segmentationDAddition.(tissues{t})(:,1)), 1);
    flagROISphere = vecnorm(segmentationDAddition.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew;
    flagROIShell = (radOld <= vecnorm(segmentationDAddition.(tissues{t}) - ROI_centers_rep, 2, 2)) & (vecnorm(segmentationDAddition.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew);
    segmentationDAddition.shell.(tissues{t}) = segmentationDAddition.(tissues{t})(flagROIShell, :);
    segmentationDAddition.(tissues{t}) = segmentationDAddition.(tissues{t})(flagROISphere, :);

    ROI_centers_rep = repmat(electCen, length(segmentationDDeletion.(tissues{t})(:,1)), 1);
    flagROISphere = vecnorm(segmentationDDeletion.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew;
    flagROIShell = (radOld <= vecnorm(segmentationDDeletion.(tissues{t}) - ROI_centers_rep, 2, 2)) & (vecnorm(segmentationDDeletion.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew);
    segmentationDDeletion.shell.(tissues{t}) = segmentationDDeletion.(tissues{t})(flagROIShell, :);
    segmentationDDeletion.(tissues{t}) = segmentationDDeletion.(tissues{t})(flagROISphere, :);

    
    ROI_centers_rep = repmat(electCen, length(segmentationI.(tissues{t})(:,1)), 1);
    flagROISphere = vecnorm(segmentationI.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew;
    flagROIShell = (radOld <= vecnorm(segmentationI.(tissues{t}) - ROI_centers_rep, 2, 2)) & (vecnorm(segmentationI.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew);
    segmentationI.shell.(tissues{t}) = segmentationI.(tissues{t})(flagROIShell, :);
    segmentationI.(tissues{t}) = segmentationI.(tissues{t})(flagROISphere, :);

    
    ROI_centers_rep = repmat(electCen, length(segmentationF.(tissues{t})(:,1)), 1);
    flagROISphere = vecnorm(segmentationF.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew;
    flagROIShell = (radOld <= vecnorm(segmentationF.(tissues{t}) - ROI_centers_rep, 2, 2)) & (vecnorm(segmentationF.(tissues{t}) - ROI_centers_rep, 2, 2) < radNew);
    segmentationF.shell.(tissues{t}) = segmentationF.(tissues{t})(flagROIShell, :);
    segmentationF.(tissues{t}) = segmentationF.(tissues{t})(flagROISphere, :);

  

end

