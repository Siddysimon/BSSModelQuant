function[] = brain_blocking(name, GeometryPost, GeometryPre, ROIAnalyze, configs, Epost, EpreFit, Epre)
datetimeRun = strrep(strrep(strrep(datestr(datetime('now')), ' ', '_'), '-', '_'), ':', '_');
titleRun = append('BB_', name, datetimeRun);
mkdir(titleRun);

brain_tissues = {'wm','gm'};
brain_tissue_idx = [1,2];
vox_corners = [];
vox_centers = [];
conductivities = [0.1400, 0.3300, 1.6540, 0.0100, 0.4630, 1.5000, 0];

for bt = 1:length(brain_tissues)
    vox_centers = [vox_centers; ROIAnalyze(20).TwoDSegmentationF.(brain_tissues{bt})];
end

dim_max = max(vox_centers);
dim_min = min(vox_centers);
% x_range = ceil(range(vox_centers(:,1)));
% y_range = ceil(range(vox_centers(:,2)));
% z_range = ceil(range(vox_centers(:,3)));
geoPost = GeometryPost.Geometry;
geoPre = GeometryPre.Geometry;
elem_centers_post = (geoPost.node(geoPost.cell(:,1),:) + geoPost.node(geoPost.cell(:,2),:) + geoPost.node(geoPost.cell(:,3),:) + geoPost.node(geoPost.cell(:,4),:))/4;
elem_centers_pre = (geoPre.node(geoPre.cell(:,1),:) + geoPre.node(geoPre.cell(:,2),:) + geoPre.node(geoPre.cell(:,3),:) + geoPre.node(geoPre.cell(:,4),:))/4;

if(~isfield(geoPost, 'field'))
    geoPost.field = GeometryPost.CI;
end

if(~isfield(geoPre, 'field'))
    geoPre.field = GeometryPre.CI;
end

%Trying to create grid with 5 cm spacing
x = dim_min(:,1):50:dim_max(:,1);
y = dim_min(:,2):50:dim_max(:,2);
z = dim_min(:,3):50:dim_max(:,3);

%Get ROI centers in grid
combined_placeholder = {x,y,z};
result_placeholder = combined_placeholder;
[result_placeholder{:}] = ndgrid(combined_placeholder{:});
ROI_centers = cell2mat(cellfun(@(m)m(:),result_placeholder,'uni',0));

%Find what ROI each elem_postent belongs in:

% %Get all elem_postents/ROI in brain
valid_tissues = [1:2];
valid_elem_post_tissue = ismember(geoPost.field, valid_tissues);
valid_elem_pre_tissue = ismember(geoPre.field, valid_tissues);
ROI_index_post = nan(length(elem_centers_post),1); % vector with index for each ROI
for i = 1:length(ROI_centers)
    ROI_index_post(vecnorm(elem_centers_post - ROI_centers(i,:),2,2) < 10 & valid_elem_post_tissue) = i; %Map elem_postent to ROI it's in
end

ROI_table = tabulate(ROI_index_post);
ROI_num_delete = find(ROI_table(:,2) < 1000); %Only use ROIs w more than 1k elem_post
ROI_index_post(ismember(ROI_index_post,ROI_num_delete)) = nan;
elem_post_ROI_idx = unique(ROI_index_post(~isnan(ROI_index_post))); %Collection of only valid ROIs
ROI_centers = ROI_centers(elem_post_ROI_idx,:);
ROI_index_pre = nan(length(elem_centers_pre),1); % vector with index for each ROI
elect_locs = cell2mat(configs(1)); %Finding distances from first (Front) electrode
for i = 1:length(ROI_centers)
    ROI_index_pre(vecnorm(elem_centers_pre - ROI_centers(i,:),2,2) < 10 & valid_elem_pre_tissue) = elem_post_ROI_idx(i); %Map elem_postent to ROI it's in
    ROI_dist_elect(i) = vecnorm(elect_locs(1,:) - ROI_centers(i,:),2,2);
end

% 
% for i = 1:length(ROI_index)
%     if isnan(ROI_index(i)) || ~valid_elem_post_tissue(i)
%         continue
%     else
%         brain_ROI_elem_post(i,:) = elem_centers_post(i,:);
%     end
% end
% brain_ROI_elem_post = brain_ROI_elem_post(~isnan(brain_ROI_elem_post),:);
% brain_ROI_elem_post = brain_ROI_elem_post(brain_ROI_elem_post nan(1,3));


%Finding volume of elem_postents
vol_all_elem_post = elemvolume(geoPost.node,geoPost.cell,'notsigned');
vol_all_elem_pre = elemvolume(geoPre.node,geoPre.cell,'notsigned');
vol_elem_post_by_ROI = zeros(length(elem_post_ROI_idx), 1);
vol_elem_pre_by_ROI = zeros(length(elem_post_ROI_idx)), 1;
for i= 1:length(elem_post_ROI_idx)
    vol_elem_post_by_ROI(i) = sum(vol_all_elem_post(elem_post_ROI_idx(i) == ROI_index_post));
    vol_elem_pre_by_ROI(i) = sum(vol_all_elem_pre(elem_post_ROI_idx(i) == ROI_index_pre));
end

%Now also filtering based on volume
ROI_num_delete = find(vol_elem_post_by_ROI < 1000); %Only use ROIs w more than 1k mm^3 of elem_post
ROI_index_post(ismember(ROI_index_post,ROI_num_delete)) = nan;
elem_post_ROI_idx = unique(ROI_index_post(~isnan(ROI_index_post))); %Collection of only valid ROIs

%Updating flags and elem_post array
flag_valid_elem_post = ~isnan(ROI_index_post) & valid_elem_post_tissue;
brain_ROI_elem_post = elem_centers_post .* flag_valid_elem_post;
elem_post_ROI_idx_valid = ROI_index_post .* flag_valid_elem_post;
flag_valid_elem_pre = ~isnan(ROI_index_pre) & valid_elem_pre_tissue;
elem_pre_ROI_idx_valid = ROI_index_pre .* flag_valid_elem_pre;
brain_ROI_elem_pre = elem_centers_pre .* flag_valid_elem_pre;

%Creating cylinder between electrode and ROIs to measure conductivities and Fc
%%Creating lines between each ROI center and electrode
lines = zeros(length(elect_locs(:,1)), length(ROI_centers), 3); %TODO: add this to the struct too
for i=1:length(elect_locs(:,1))
    for j=1:length(ROI_centers)
        lines(i,j,:) = (elect_locs(i,:) - ROI_centers(j,:))/vecnorm(elect_locs(i,:) - ROI_centers(j,:),2,2);
    end
end

% sphere_dist = 2.5;
% for i=1:length(elect_locs(:,1))
%     for j=1:length(ROI_centers)
%         center_dists = 0:sphere_dist:vecnorm(elect_locs(i,:) - ROI_centers(j,:),2,2);
%         incrementsRaw = lines(i,j,:) .* center_dists;
%         increments = reshape(incrementsRaw,[length(incrementsRaw(1,:,1)), 3]);
%         cylindersPost(i,j).spheres = ROI_centers(j,:) + increments;
%         %Todo: how to efficiently do subtraction here:
%         elements_cyl = [];
%         conds_cyl = [];
%         for k=1:length(cylindersPost(i,j).spheres(:,1))
%             elem_flag = vecnorm(cylindersPost(i,j).spheres(k,:) - elem_centers_post,2,2) > 5;
%             elements_cyl = [elements_cyl; elem_centers_post(elem_flag,:)];
%             conds_cyl = [conds_cyl; geoPost.field(elem_flag,:)];
%         end
%         cylindersPost(i,j).elements = elements_cyl;
%         cylindersPost(i,j).tissues = conds_cyl;
%     end
% end

sphere_dist = 2.5;
for i=1:length(elect_locs(:,1))
    for j=1:length(ROI_centers)
        center_dists = 0:sphere_dist:vecnorm(elect_locs(i,:) - ROI_centers(j,:),2,2);
        incrementsRaw = lines(i,j,:) .* center_dists;
        increments = reshape(incrementsRaw,[length(incrementsRaw(1,:,1)), 3]);
        cylindersPost(i,j).spheres = ROI_centers(j,:) + increments;
        cylindersPre(i,j).spheres = ROI_centers(j,:) + increments;
        elem_flag_post = zeros(size(elem_centers_post(:,1)));
        elem_flag_pre = zeros(size(elem_centers_pre(:,1)));
        for k=1:length(cylindersPost(i,j).spheres(:,1))
            elem_flag_post = elem_flag_post | vecnorm(cylindersPost(i,j).spheres(k,:) - elem_centers_post,2,2) < 5;
            elem_flag_pre = elem_flag_pre | vecnorm(cylindersPre(i,j).spheres(k,:) - elem_centers_pre,2,2) < 5;
        end
        cylindersPost(i,j).elements =  elem_centers_post(elem_flag_post,:);
        cylindersPost(i,j).conds = conductivities(geoPost.field(elem_flag_post,:))';
        cylindersPost(i,j).field = Epost(elem_flag_post, :);
        cylindersPost(i,j).avgDist = mean(vecnorm(elect_locs(i,:) - elem_centers_post(elem_flag_post,:),2,2));
        cylindersPre(i,j).field = Epre(elem_flag_pre, :);
        cylindersPre(i,j).elements =  elem_centers_pre(elem_flag_pre,:);
        cylindersPre(i,j).conds = conductivities(geoPre.field(elem_flag_pre,:))';
        cylindersPre(i,j).avgDist = mean(vecnorm(elect_locs(i,:) - elem_centers_pre(elem_flag_pre,:),2,2));
        cylindersDiff(i,j).conds = mean(cylindersPost(i,j).conds) - mean(cylindersPre(i,j).conds);
        cylindersDiff(i,j).condsAbs = abs(mean(cylindersPost(i,j).conds) - mean(cylindersPre(i,j).conds));
        cylindersDiff(i,j).fieldMag = mean(vecnorm(cylindersPost(i,j).field,2,2)) - mean(vecnorm(cylindersPre(i,j).field,2,2));
        cylindersDiff(i,j).field = mean(cylindersPost(i,j).field) - mean(cylindersPre(i,j).field);
        cylindersDiff(i,j).condRel = cylindersDiff(i,j).conds / mean(cylindersPre(i,j).conds);
        cylindersDiff(i,j).condRelAbs = cylindersDiff(i,j).condsAbs / mean(cylindersPre(i,j).conds);
        cylindersDiff(i,j).fieldRel = cylindersDiff(i,j).fieldMag / vecnorm(mean(cylindersPre(i,j).field),2,2);
        cylindersDiff(i,j).avgDist = mean([cylindersPost(i,j).avgDist; cylindersPre(i,j).avgDist]);

    end
end

% for i=1:length(elect_locs(:,1))
%     for j=1:length(ROI_centers)
%         center_dists = 0:sphere_dist:vecnorm(elect_locs(i,:) - ROI_centers(j,:),2,2);
%         incrementsRaw = lines(i,j,:) .* center_dists;
%         increments = reshape(incrementsRaw,[length(incrementsRaw(1,:,1)), 3]);
%         cylindersPost(i,j).spheres = ROI_centers(j,:) + increments;
%         elem_flag = zeros(size(elem_centers_post(:,1)));
%         for k=1:length(cylindersPost(i,j).spheres(:,1))
%             elem_flag = elem_flag | vecnorm(cylindersPost(i,j).spheres(k,:) - elem_centers_post,2,2) &lt; 5;
%         end
%         cylindersPost(i,j).elements =  elem_centers_post(elem_flag,:);
%         cylindersPost(i,j).tissues = geoPost.field(elem_flag,:);
%     end
% end

EMag = abs(vecnorm(Epost, 2,2) - vecnorm(EpreFit, 2,2));
ERel = 100 * (EMag ./ vecnorm(Epost, 2,2));
placeholderEMag = EMag .* flag_valid_elem_post;
placeholderERel = ERel .* flag_valid_elem_post;

%Collecting Field change and conductivity
collectEDiffMag = zeros(length(elem_post_ROI_idx_valid),1);
collectEDiffRel = zeros(length(elem_post_ROI_idx_valid),1);
collectCondPost = zeros(length(elem_post_ROI_idx_valid),1);
tissue_cond_post = GeometryPost.CL .* flag_valid_elem_post;
numPostElem_per_ROI = tabulate(elem_post_ROI_idx_valid);
for i = 1:length(elem_post_ROI_idx_valid)
    ROI_num = elem_post_ROI_idx_valid(i,:);
    if ~isnan(ROI_num)
        numElem = numPostElem_per_ROI(ROI_num,2);
        collectEDiffMag(ROI_num,:) = collectEDiffMag(ROI_num,:) + placeholderEMag(i);
        collectEDiffRel(ROI_num,:) = collectEDiffRel(ROI_num,:) + placeholderERel(i);
        collectCondPost(ROI_num,:) = (numElem * collectCondPost(elem_post_ROI_idx_valid(i),:) + tissue_cond_post(i))/numElem;
    end
end
%collectCondPost = collectCondPost(~(collectCondPost == 0)); %%Check if needed and if 159 - 174 is in correct place

closestElectIdx = zeros(length(ROI_centers));
closestDist2Elect = zeros(length(ROI_centers));
closestElectIdx = zeros(length(ROI_centers),1);
closestDist2Elect = zeros(length(ROI_centers),1);

for i=1:length(elect_locs(:,1))
    for j=1:length(ROI_centers)
        dist2Elect = vecnorm(elect_locs(i,:) - ROI_centers(j,:),2,2);
        lines(i,j,:) = (elect_locs(i,:) - ROI_centers(j,:)) / dist2Elect;

        if(closestDist2Elect(j) == 0 || dist2Elect < closestDist2Elect(j))
            closestDist2Elect(j) = dist2Elect;
            closestElectIdx(j) = i;
        end
    end
end


i=1;j=length(ROI_centers);
plot3(cylindersPost(i,j).elements(:,1), cylindersPost(i,j).elements(:, 2), cylindersPost(i,j).elements(:,3),'o','MarkerSize',10)
hold on
% test = ROI_centers(1:2, :);
% pROI = scatter3(test(:, 1), test(:, 2), test(:,3),'rs');

pROI = plot3(ROI_centers(j, 1), ROI_centers(j, 2), ROI_centers(j,3),'p','MarkerSize',50);
pElect = scatter3(elect_locs(:, 1), elect_locs(:, 2), elect_locs(:,3),'g*');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/Cylinder_Check', name))
saveas(gcf,append(titleRun, '/Cylinder_Check', name, '.png'))
close all
% p.MarkerSize = 10;
% p.Color = "blue";

p1 = scatter3(brain_ROI_elem_post(:,1), brain_ROI_elem_post(:,2), brain_ROI_elem_post(:,3),'g');
scatter3(elem_centers_post(valid_elem_post_tissue,1), elem_centers_post(valid_elem_post_tissue,2), elem_centers_post(valid_elem_post_tissue,3),'.');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/ROI_Check', name))
saveas(gcf,append(titleRun, '/ROI_Check', name, '.png'))
close all



conds = [];
field = [];
for j = 1:length(ROI_centers)
    conds = [conds; cylindersDiff(closestElectIdx(j),j).conds];
    field = [field; cylindersDiff(closestElectIdx(j),j).fieldMag];

    if closestElectIdx(j) == 1
        scatter(cylindersDiff(closestElectIdx(j),j).conds, cylindersDiff(closestElectIdx(j),j).fieldMag, 'DisplayName','Front Electrode', MarkerEdgeColor='r')
    else
        scatter(cylindersDiff(closestElectIdx(j),j).conds, cylindersDiff(closestElectIdx(j),j).fieldMag, 'DisplayName','Back Electrode', MarkerEdgeColor='b')
    end
    hold on
end
%scatter(conds, field)

xlabel('Average cylinder conductivity change (S/m)');
ylabel('Average magnitude of field change (V/m)');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/Cylinder_Cond_vs_field_change', name))
saveas(gcf,append(titleRun, '/Cylinder_Cond_vs_field_change', name, '.png'))
close all

conds = [];
field = [];
for j = 1:length(ROI_centers)
    conds = [conds; cylindersDiff(closestElectIdx(j),j).condRel];
    field = [field; cylindersDiff(closestElectIdx(j),j).fieldRel];

    if closestElectIdx(j) == 1
        scatter(cylindersDiff(closestElectIdx(j),j).condRel, cylindersDiff(closestElectIdx(j),j).fieldRel, 'DisplayName','Front Electrode', 'MarkerFaceColor', 'r')
    else
        scatter(cylindersDiff(closestElectIdx(j),j).condRel, cylindersDiff(closestElectIdx(j),j).fieldRel, 'DisplayName','Back Electrode', 'MarkerFaceColor', 'b')
    end
    hold on
end


% for j = 1:length(ROI_centers)

% conds = [conds; cylindersDiff(closestElectIdx(j),j).condRel];


% field = [field; cylindersDiff(closestElectIdx(j),j).fieldRel];

% end
% scatter(conds, field)
%legend

xlabel('Average relative cylinder conductivity change');
ylabel('Average magnitude of relative field change');

hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/Rel_Cylinder_Cond_vs_field_change', name))
saveas(gcf,append(titleRun, '/Rel_Cylinder_Cond_vs_field_change', name, '.png'))
close all



conds = [];
field = [];
dists = [];
for j = 1:length(ROI_centers)
    conds = [conds; cylindersDiff(closestElectIdx(j),j).conds];
    field = [field; cylindersDiff(closestElectIdx(j),j).fieldMag];
    dists = [dists; cylindersDiff(closestElectIdx(j),j).avgDist];
end
%scatter(conds, field)
scatter(conds, field, [],dists, "filled")
c = colorbar;
c.Label.String = "Distance to closest electrode (mm)";
c.Label.Rotation = 270;
xlabel('Average cylinder conductivity change (S/m)');
ylabel('Average magnitude of field change (V/m)');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/Closest_elect_cylinder_Cond_vs_field_change', name))
saveas(gcf,append(titleRun, '/Closest_elect_cylinder_Cond_vs_field_change', name, '.png'))
close all

conds = [];
field = [];
dists = [];
for j = 1:length(ROI_centers)
    conds = [conds; cylindersDiff(closestElectIdx(j),j).condRel];
    field = [field; cylindersDiff(closestElectIdx(j),j).fieldRel];
    dists = [dists; cylindersDiff(closestElectIdx(j),j).avgDist];
end
%scatter(conds, field)
scatter(conds, field, [],dists, "filled")
c = colorbar;
c.Label.String = "Distance to closest electrode (mm)";
c.Label.Rotation = 270;
xlabel('Relative cylinder conductivity change');
ylabel('Relative magnitude of field change');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/Closest_elect_Rel_cylinder_Cond_vs_field_change', name))
saveas(gcf,append(titleRun, '/Closest_elect_Rel_cylinder_Cond_vs_field_change', name, '.png'))
close all


conds = [];
field = [];
dists = [];
for j = 1:length(ROI_centers)
    conds = [conds; cylindersDiff(closestElectIdx(j),j).condRelAbs];
    field = [field; cylindersDiff(closestElectIdx(j),j).fieldRel];
    dists = [dists; closestDist2Elect(j)];
end
%scatter(conds, field)
scatter(conds, field, [],dists, "filled")
c = colorbar;
c.Label.String = "Distance from ROI to closest electrode (mm)";
c.Label.Rotation = 270;
xlabel('Relative cylinder conductivity change');
ylabel('Relative magnitude of field change');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/ROI2elect_dist_Rel_cylinder_Cond_vs_field_change', name))
saveas(gcf,append(titleRun, '/ROI2elect_dist_Rel_cylinder_Cond_vs_field_change', name, '.png'))
close all


conds = [];
field = [];
dists = [];
for j = 1:length(ROI_centers)
    conds = [conds; cylindersDiff(closestElectIdx(j),j).condsAbs];
    field = [field; cylindersDiff(closestElectIdx(j),j).fieldMag];
    dists = [dists; closestDist2Elect(j)];
end
%scatter(conds, field)
scatter(conds, field, [],dists, "filled")
c = colorbar;
c.Label.String = "Distance from ROI to closest electrode (mm)";
c.Label.Rotation = 270;
xlabel('Avg cylinder conductivity change');
ylabel('Avg magnitude of field change');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/ROI2elect_dist_cylinder_Cond_vs_field_change', name))
saveas(gcf,append(titleRun, '/ROI2elect_dist_cylinder_Cond_vs_field_change', name, '.png'))
close all



conds = [];
field = [];
dists = [];
for j = 1:length(ROI_centers)
    conds = [conds; cylindersDiff(closestElectIdx(j),j).condsAbs];
    field = [field; cylindersDiff(closestElectIdx(j),j).fieldMag];
    dists = [dists; closestDist2Elect(j)];
end
%scatter(conds, field)
scatter(dists, field, [],conds, "filled")
c = colorbar;
c.Label.String = "Absolute cond. Change (S/m)";
c.Label.Rotation = 270;
xlabel('Distance from ROI to closest electrode (mm)');
ylabel('Avg magnitude of field change (V/m)');
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/ROI2elect_dist_vs_field_change', name))
saveas(gcf,append(titleRun, '/ROI2elect_dist_vs_field_change', name, '.png'))
close all
end

