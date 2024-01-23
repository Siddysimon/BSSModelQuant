 function [EAngle, EMagnitude, ROI, Segmentation] = Analyze(electrodes, ROIcenter, geometrySTR, pre, post, mask_prepI, mask_prepF, ROI_radii, name)


%Highest level field and voxel analysis program

%Simple run with no electric field consideration:
if(isempty(electrodes))
    warning('Only voxel differences will be calculated with no radius');
    EAngle = [];
    EMagnitude = [];
    ROI = [];
    [Segmentation(1).diff, Segmentation(1).final, Segmentation(1).initial] = get_voxel_data(mask_prepI, mask_prepF, geometrySTR);
    
else
    %Full run:
    
    Geometry = load(geometrySTR);
    if isfield(Geometry,'HeadModel')
        Geometry = Geometry.HeadModel;
    elseif isfield(Geometry,'Geometry')
        GeoTot = Geometry;
        Geometry = Geometry.Geometry;
    end
    
    %load voxel differences:
    electSize = size(electrodes);
    numElectrodes = electSize(1);
    element_centers = (Geometry.node(Geometry.cell(:,1),:) + Geometry.node(Geometry.cell(:,2),:) + Geometry.node(Geometry.cell(:,3),:) + Geometry.node(Geometry.cell(:,4),:))/4;
    
    numElements = length(element_centers);
    
    %Fill arrays with zeros for now:
    dist_to_center = zeros(numElements, numElectrodes);
    E_vector = load(pre, 'E_vector') %loading early to get size of E
    
    if isstruct(E_vector)
        E_vector = E_vector.E_vector;
    end
    
    EMagnitude = zeros(2, length(E_vector));
    EAngle = zeros(2, length(E_vector));
    
    %Function to select the voxel data
    [Segmentation.diff, Segmentation.final, Segmentation.initial] = get_voxel_data(mask_prepI, mask_prepF, Geometry);
    
    
    
    %Loop through each electrode configuration and collect field data
    for e = 1:numElectrodes
        
        %Load appropriate field data
        if (e == 1)
            Epre = E_vector;
            E_vector = load(post,'E_vector')
            if isstruct(E_vector)
                Epost = E_vector.E_vector;
            else
                Epost = E_vector;
            end
            E_vector_diff = Epost - Epre;
            EMagnitude(1,:) =  vecnorm(Epost, 2, 2) - vecnorm(Epre, 2, 2); %Get difference in field magnitudes
            EAngle(1,:) = acosd(dot(Epost,Epre, 2)./(vecnorm(Epost, 2, 2).*vecnorm(Epre, 2, 2))); %Get difference in field direction
            
        end
        
        
        for i = 1:length(ROI_radii)
            ROI(e,i).center = electrodes(e,:); %The electrode is the ROI center (for now)
            
            %Convert voxel data to 2D (nx3) arrays for specific ROIs
            if i == 1
               [ROI(e,i).TwoDSegmentationI, ROI(e,i).TwoDSegmentationF, ROI(e,i).TwoDSegmentationD, ROI(e,i).TwoDSegmentationDDel, ROI(e,i).TwoDSegmentationDAdd] = scale_and_center(Segmentation.initial, Segmentation.final, Segmentation.diff, Geometry, ROIcenter, ROI_radii(i), 0);
            else if i == 5
               [ROI(e,i).TwoDSegmentationI, ROI(e,i).TwoDSegmentationF, ROI(e,i).TwoDSegmentationD, ROI(e,i).TwoDSegmentationDDel, ROI(e,i).TwoDSegmentationDAdd] = scale_and_center(Segmentation.initial, Segmentation.final, Segmentation.diff, Geometry, ROIcenter, ROI_radii(i), ROI_radii(i - 1));
            else
                [ROI(e,i).TwoDSegmentationI, ROI(e,i).TwoDSegmentationF, ROI(e,i).TwoDSegmentationD, ROI(e,i).TwoDSegmentationDDel, ROI(e,i).TwoDSegmentationDAdd] = scale_and_center(Segmentation.initial, Segmentation.final, Segmentation.diff, Geometry, ROIcenter, ROI_radii(i), ROI_radii(i - 1));
            end
            %Generate electric field averages for different ROIs
            ROI_centers_repeated = repmat(ROI(e,i).center, numElements, 1);
            dist_elem_ROI_temp = element_centers - ROI_centers_repeated; % vector from each elem's center to the ROI center
            dist_to_center(:, i) = vecnorm(dist_elem_ROI_temp,2,2); % absolute distance from each elem's center to the ROI center
            if i == 1
               ROI(e,i).Eflag = dist_to_center(:, i) < ROI_radii(i); %Isolating points of field w/in ROI
            else
                ROI(e,i).Eflag = dist_to_center(:, i) >= ROI_radii(i-1) & dist_to_center(:, i) < ROI_radii(i); %Isolating points of field w/in ROI
            end
            ROI(e,i).noTissueE = mean(E_vector_diff(ROI(e,i).Eflag)); %Get average electric field vector w/in ROI
        end
        
    end
    
    %Visualize the results:
     Visualize(name, Epost, Epre, GeoTot, ROI, Segmentation, ROI_radii, electrodes)
    
    
    
    
end

for fn = fieldnames(Segmentation.diff)'
    Segmentation.voxChangeSum.(fn{1}) = sum(sum(sum(Segmentation.diff.(fn{1}) ~= 0)));
end



end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




