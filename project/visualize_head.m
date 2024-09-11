function [] = visualize_head(name,Epost,Epre, Geometry, ROI, Segmentation, ROI_radii, electrodes)

% %Highest level field and voxel analysis program
% Make Run Folder
datetimeRun = strrep(strrep(strrep(datestr(datetime('now')), ' ', '_'), '-', '_'), ':', '_');
titleRun = append('Visualize_', name, datetimeRun);
mkdir(titleRun);
voxWeight3d = Segmentation.final.voxWeight(1,:) * Segmentation.final.voxWeight(2,:) * Segmentation.final.voxWeight(3,:); 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tissues = {'wm','gm','csf','bone','skin','eyes','air'};
CI = Geometry.CI;
Geometry = Geometry.Geometry;
tissues = tissues(1:length(unique(CI)));
% Generate proof of proper centering
element_centers_tot = (Geometry.node(Geometry.cell(:,1),:) + Geometry.node(Geometry.cell(:,2),:) + Geometry.node(Geometry.cell(:,3),:) + Geometry.node(Geometry.cell(:,4),:))/4;
nodeDownSample1 = Geometry.node([1:10:length(Geometry.node)], :);
conversion = ROI(1, 10).TwoDSegmentationD.skin;
EDiff = Epost - Epre;

plot3(conversion(:,1), conversion(:,2), conversion(:,3), 'b.'); %or whatever marker you want
hold on

quiver3(element_centers_tot(:,1),element_centers_tot(:,2),element_centers_tot(:,3),Epost(:, 1),Epost(:, 2),Epost(:, 3), 'r')
quiver3(element_centers_tot(:,1),element_centers_tot(:,2),element_centers_tot(:,3),EDiff(:, 1),EDiff(:, 2),EDiff(:, 3), 'r')

% saveas(gcf,append(titleRun, '/3D_centering_proof')) --> Large file, data
% not very useful to begin with

saveas(gcf,append(titleRun, '/3D_centering_proof.png'))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get voxel information for slice-by-slice viewing, divided by
%concuctivities or tissue type


sumI = zeros(size(Segmentation(1).initial.skin));
sumF = zeros(size(Segmentation(1).initial.skin));
sumD = zeros(size(Segmentation(1).initial.skin));

condsumI = zeros(size(Segmentation(1).initial.skin));
condsumF = zeros(size(Segmentation(1).initial.skin));
condsumD = zeros(size(Segmentation(1).initial.skin));
conductivities = [0.1400, 0.3300, 1.6540, 0.0100, 0.4630, 1.5000, 0];

for t = 1:length(tissues)
    sumI = sumI + Segmentation(1).initial.(tissues{t}) * t;
    sumF = sumF + Segmentation(1).final.(tissues{t}) * t;
    sumD = sumD + Segmentation(1).diff.(tissues{t}) * t;
    
    condsumI = condsumI + Segmentation(1).initial.(tissues{t}) * conductivities(t);
    condsumF = condsumF + Segmentation(1).final.(tissues{t}) * conductivities(t);
end
condsumD = condsumF - condsumI;

% By tissue index
tiledlayout(3,1)
figure('units','normalized', 'OuterPosition', [0 0 1 1])

nexttile
imagesc(sumI(:,:,130))
bottom = min(min(sumI(:)), min(sumF(:)));
top = max(max(sumI(:)), max(sumF(:)));
colorbar;
caxis manual
caxis([bottom top])
title(['Initial mask, min=' num2str(min(sumI(:))) ', max=' num2str(max(sumI(:)))])
set(gca,'FontSize',14)
% c.Ticks=0:12;


nexttile
imagesc(sumF(:,:,130))
colorbar;
caxis manual
caxis([bottom top])
title(['Final mask, min=' num2str(min(sumF(:))) ', max=' num2str(max(sumF(:)))])% c.Ticks=0:12;
set(gca,'FontSize',14)

nexttile
imagesc(sumD(:,:,130))
c=colorbar;
title(['Change in masks, min=' num2str(min(sumD(:))) ', max=' num2str(max(sumD(:)))])
set(gca,'FontSize',14)

% c.Ticks=0:12
% print(['Axial slices of initial, final and difference voxel masks, respectively.'] ) %Add name

close all

%For conductivities:

%Initial
tiledlayout(3,1)
figure('units','normalized', 'OuterPosition', [0 0 1 1])

nexttile
imagesc(condsumI(:,:,130))
bottom = min(min(condsumI(:)), min(condsumF(:)));
top = max(max(condsumI(:)), max(condsumF(:)));
colorbar;
caxis manual 
caxis([bottom top])
title(['Initial mask, min=' num2str(min(condsumI(:))) ', max=' num2str(max(condsumI(:)))])
set(gca,'FontSize',14)


%Final
nexttile
imagesc(condsumF(:,:,130))
colorbar;
caxis manual
caxis([bottom top])
title(['Final mask, min=' num2str(min(condsumF(:))) ', max=' num2str(max(condsumF(:)))])% c.Ticks=0:12;
set(gca,'FontSize',14)

%Difference
nexttile
imagesc(condsumD(:,:,130))
c=colorbar;
title(['Change in masks, min=' num2str(min(condsumD(:))) ', max=' num2str(max(condsumD(:)))])
set(gca,'FontSize',14)

saveas(gcf,append(titleRun, '/Axial_Slices', name))
saveas(gcf,append(titleRun, '/Axial_Slices', name, '.png'))

close all

% %Difference
% nexttile
% ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1) * 
% imagesc(condsumD(:,:,130))
% c=colorbar;
% title(['Change in masks, min=' num2str(min(condsumD(:))) ', max=' num2str(max(condsumD(:)))])
% set(gca,'FontSize',14)
% 
% saveas(gcf,append(titleRun, '/Weighted_Cond', name))
% saveas(gcf,append(titleRun, '/Weighted_Cond', name, '.png'))






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get ROI shells per tissue type with respect to the electric field:

for r = length(ROI_radii(:,1))
    
    %Get average field change in ROI shell by tissue and for all
    EChanged(r).noTissue = []; 
    for t = 1:length(tissues)
        EChanged(r).(tissues{t}) = [];
    end
    
    electCen = electrodes(1,:);
    centers_rep = repmat(electCen, length(Geometry.node), 1);
    oldRad = 0;

    for rad = ROI_radii(r, :)
        flagNode = oldRad < vecnorm(Geometry.node - centers_rep, 2, 2) & vecnorm(Geometry.node - centers_rep, 2, 2) <= rad;
        flag = flagNode(Geometry.cell(:,1)) & flagNode(Geometry.cell(:,2)) & flagNode(Geometry.cell(:,3));
        E = abs(vecnorm(Epost, 2,2) - vecnorm(Epre, 2,2));
        placeholder = E(flag,:);
        EChanged(r).noTissue= [EChanged(r).noTissue, mean(placeholder)];
        for f = 1:length(tissues)
            placeholderTissue = E(flag & CI == f);
            EChanged(r).(tissues{f})= [EChanged(r).(tissues{f}), mean(placeholderTissue)];
        end
        
        oldRad = rad;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quantify voxel change by ROI shell, visualizing with avg change of field

%Calculating total voxel change
for t = 1:length(tissues)

    for r = 1:length(ROI_radii)
        voxChangeCombined(r, t) = voxWeight3d * size(ROI(1,r).TwoDSegmentationD.shell.(tissues{t}),1);
    end
    plot(ROI_radii, voxChangeCombined(:,t))
    hold on
end
legend(tissues);
xlabel('ROI radius (voxels)')
ylabel('Volume of tissue change (mm3)')
title(['Voxel Changes as Function of Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/DistVsVoxChange', name))
saveas(gcf,append(titleRun, '/DistVsVoxChange', name, '.png'))
close all


%Voxel additions only:
for t = 1:length(tissues)
    for r = 1:length(ROI_radii)
        voxAddCombined(r, t) = voxWeight3d * size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1);
    end
    plot(ROI_radii, voxAddCombined(:,t))
    hold on
end
legend(tissues);
xlabel('ROI radius (voxels)')
ylabel('Volume of voxel additions (mm3)')
title(['Voxel Additions as Function of Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/DistVsVoxAdd', name))
saveas(gcf,append(titleRun, '/DistVsVoxAdd', name, '.png'))
close all

%Voxel Deletions only:
for t = 1:length(tissues)
    for r = 1:length(ROI_radii)
        voxDelCombined(r, t) = voxWeight3d * size(ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t}),1);
    end
    plot(ROI_radii, voxDelCombined(:,t))
    hold on
end
legend(tissues);
xlabel('ROI radius (voxels)')
ylabel('Volume of voxel deletions (mm3)')
title(['Voxel Deletions as Function of Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/DistVsVoxDel', name))
saveas(gcf,append(titleRun, '/DistVsVoxDel', name, '.png'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quantify voxel change by conductivity, finding average and visualizing
%with change of field over cumulative ROI spheres (not shells)

%Vox ROI Sphere overall conductivity change
lenROI = length(ROI);
voxChangeCombinedCond = zeros(1, lenROI);
for r = 1:lenROI
    for t = 1:length(tissues)
        voxChangeCombinedCond(r) = voxChangeCombinedCond(r) + size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1) * conductivities(t);
    end
end

% ROI Sphere conductivity change just from vox additions, deemed unneeded
% voxAddCombinedCond = zeros(1, lenROI);
% for r = 1:lenROI
%     for t = 1:length(tissues)
%         voxAddCombinedCond(r) = voxAddCombinedCond(r) + size(ROI(1,r).TwoDSegmentationDAdd.(tissues{t}),1) * conductivities(t);
%     end
% end

% ROI Sphere conductivity change just from vox deletions, deemed unneeded
% voxDelCombinedCond = zeros(1, lenROI);
% for r = 1:lenROI
%     for t = 1:length(tissues)
%         voxDelCombinedCond(r) = voxDelCombinedCond(r) + size(ROI(1,r).TwoDSegmentationDDel.(tissues{t}),1) * conductivities(t);
%     end
% end

EChanged.noTissue(isnan(EChanged.noTissue)) = 0;
condChangeRep = repmat(voxChangeCombinedCond, lenROI, 1);
% condAddRep = repmat(voxAddCombinedCond, lenROI, 1); Add/Delete deemed
% unneccesary
% condDelRep = repmat(voxDelCombinedCond, lenROI, 1);
eChangeRep = repmat(EChanged.noTissue, lenROI, 1);
ratioCondE = condChangeRep' ./ eChangeRep;
imagesc(ratioCondE)
colorbar;


for r = 1:length(ROI)
    cmap = jet(length(ROI)); %20 colors
    scatter(condChangeRep(r,:), eChangeRep(r,:),40, cmap, 'filled');
%     hold on
%     scatter(condAddRep(r,:), eChangeRep(r,:))
%     scatter(condDelRep(r,:), eChangeRep(r,:))
end
c = colorbar;
c.Label.String = 'Distance from electrode to end of sphere (vox)';
caxis manual
caxis([ROI_radii(1), ROI_radii(length(ROI_radii))])
xlabel('ROI Conductivity Change (S/m)');
ylabel('Magnitude of Field Change (V/m)');
title(['Conductivity Change vs Field Change in Spheres'])

% legend({('All Changes'), ('Additions'), ('Deletions')});
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Cond_vs_Field_Change', name))
saveas(gcf,append(titleRun, '/Cond_vs_Field_Change', name, '.png'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get sum of vox change w/in ROI shell and visualize, along with field change


% For 2D:
tiledlayout(3,1)
nexttile
%Visualize radius vs vox changes by tissue
for t = 1:length(tissues)
    for r = 1:length(ROI)
        voxChangeCombined(r, t) = voxWeight3d * size(ROI(1,r).TwoDSegmentationD.shell.(tissues{t}),1);
    end
    plot(ROI_radii, voxChangeCombined(:,t))
    xlabel('ROI Radius')
    ylabel('Number of Voxels Changed')
    title(['Distance vs. Voxels Changed'])

    hold on
end
legend(tissues);
set(gca,'FontSize',7.5)

nexttile

%Visualize field change by tissue type compared to vox changes
for t = 1:length(tissues)
    plot(voxChangeCombined(:,t), EChanged.(tissues{t}))
    xlabel('Volume of Voxels Changed (mm3)')
    ylabel('Magnitude of Field Change (V/m)')
    title(['Voxel Changes vs. Field Magnitude Change'])
    hold on
end
legend(tissues);
set(gca,'FontSize',7.5)

nexttile

%Visualize field change by radii
plot(ROI_radii, EChanged.noTissue)
xlabel('ROI Radius')
ylabel('Magnitude of Field Change (V/m)')
title(['Field Magnitude Change Over Distance'])
set(gca,'FontSize',7.5)

saveas(gcf,append(titleRun, '/ROI_Shell_Vox_Field_Radii_Plots', name))
saveas(gcf,append(titleRun, '/ROI_Shell_Vox_Field_Radii_Plots', name, '.png'))
close all


% For just additions:
tiledlayout(3,1)
nexttile
%Visualize radius vs vox changes by tissue
for t = 1:length(tissues)
    for r = 1:length(ROI)
        voxAddCombined(r, t) = voxWeight3d * size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1);
    end
    plot(ROI_radii, voxAddCombined(:,t))
    xlabel('ROI Radius')
    ylabel('Volume of Voxels Added (mm3)')
    title(['Distance vs. Voxel Additions'])
    hold on
end
legend(tissues);
set(gca,'FontSize',7.5)

nexttile

%Visualize field change by tissue type compared to vox additions
for t = 1:length(tissues)
    plot(voxAddCombined(:,t), EChanged.(tissues{t}))
    xlabel('Volume of Voxels Added (mm3)')
    ylabel('Magnitude of Field Change (V/m)')
    title(['Voxel Additions vs. Field Magnitude Change'])
    hold on
end
legend(tissues);
set(gca,'FontSize',7.5)

nexttile
%Visualize field change by radii
plot(ROI_radii, EChanged.noTissue)
xlabel('ROI Radius')
ylabel('Magnitude of Field Additions (V/m)')
title(['Field Magnitude Change Over Distance'])
set(gca,'FontSize',7.5)

saveas(gcf,append(titleRun, '/ROI_Shell_Vox_Add_Field_Radii_Plots', name))
saveas(gcf,append(titleRun, '/ROI_Shell_Vox_Add_Field_Radii_Plots', name, '.png'))
close all

% For just deletions:
tiledlayout(3,1)
nexttile
%Visualize radius vs vox deletions by tissue
for t = 1:length(tissues)
    for r = 1:length(ROI)
        voxDelCombined(r, t) = voxWeight3d * size(ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t}),1);
    end
    plot(ROI_radii, voxDelCombined(:,t))
    xlabel('ROI Radius')
    ylabel('Volume of Voxels Deleted (mm3)')
    title(['Distance vs. Voxel Deletions'])

    hold on
end
legend(tissues);
set(gca,'FontSize',7.5)

nexttile

%Visualize field change by tissue type compared to vox additions
for t = 1:length(tissues)
    plot(voxDelCombined(:,t), EChanged.(tissues{t}))
    xlabel('Volume of Voxels Deleted (mm3)')
    ylabel('Magnitude of Field Change (V/m)')
    title(['Voxel Additions vs. Field Magnitude Change'])
    hold on
end
legend(tissues);
set(gca,'FontSize',7.5)

nexttile
%Visualize field change by radii
plot(ROI_radii, EChanged.noTissue)
xlabel('ROI Radius')
ylabel('Magnitude of Field Additions (V/m)')
title(['Field Magnitude Change Over Distance'])
set(gca,'FontSize',7.5)

saveas(gcf,append(titleRun, '/ROI_Shell_Vox_Del_Field_Radii_Plots', name))
saveas(gcf,append(titleRun, '/ROI_Shell_Vox_Del_Field_Radii_Plots', name, '.png'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Avg Conductivity Change Over all Voxels in ROI shell

% Addition:
voxAddCombinedCond = zeros(1, length(ROI));
for r = 1:length(ROI)
    totVox(r) = 0;
    for t = 1:length(tissues)
        
        %ROI shells
        newAddSize = size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1);
        newTotSize = size(ROI(1,r).TwoDSegmentationF.shell.(tissues{t}),1);
        
        voxAddCombinedCond(r) = voxAddCombinedCond(r) + newAddSize * conductivities(t);
        totVox(r) = newTotSize + totVox(r);
    end
    voxAddCombinedCond(r) = voxAddCombinedCond(r) / totVox(r);
end

scatter(voxAddCombinedCond, EChanged.noTissue)
xlabel('Avg Conductivity Additions Over All Voxels (S/m)')
ylabel('Magnitude of Field Change (V/m)')
title(['Avg Voxel Conductivity Addition vs. Field Magnitude Change'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Vox_Tot_Add_Cond_vs_Field', name))
saveas(gcf,append(titleRun, '/Vox_Tot_Add_Cond_vs_Field', name, '.png'))
close all

plot(ROI_radii, voxAddCombinedCond)
xlabel('ROI radius')
ylabel('Avg Conductivity Additions Over All Voxels (S/m)')
title(['Avg Voxel Conductivity Addition Over Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Tot_Add_Cond', name))
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Tot_Add_Cond', name, '.png'))
close all

% Deletion:
voxDelCombinedCond = zeros(1, length(ROI));
for r = 1:length(ROI)
    totVox(r) = 0;
    for t = 1:length(tissues)
        
        %ROI shells
        newDelSize = size(ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t}),1);
        newTotSize = size(ROI(1,r).TwoDSegmentationF.shell.(tissues{t}),1);
        
        voxDelCombinedCond(r) = voxDelCombinedCond(r) + newDelSize * conductivities(t);
        totVox(r) = newTotSize + totVox(r);
    end
    voxDelCombinedCond(r) = voxDelCombinedCond(r) / totVox(r);
end

scatter(voxDelCombinedCond, EChanged.noTissue)
xlabel('Avg Conductivity Deletions Over All Voxels (S/m)')
ylabel('Magnitude of Field Change (V/m)')
title(['Avg Voxel Conductivity Deletion vs. Field Magnitude Change'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Vox_Tot_Del_Cond_vs_Field', name))
saveas(gcf,append(titleRun, '/Vox_Tot_Del_Cond_vs_Field', name, '.png'))
close all

plot(ROI_radii, voxDelCombinedCond)
xlabel('ROI radius')
ylabel('Avg Conductivity Deletions Over All Voxels (S/m)')
title(['Avg Voxel Conductivity Deletion Over Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Tot_Del_Cond', name))
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Tot_Del_Cond', name, '.png'))
close all

% Total Change:
voxChangeCombinedCond = zeros(1, length(ROI));

for r = 1:length(ROI)
    totVox(r) = 0;
    for t = 1:length(tissues)
        
        %ROI shells
        newChangeSize = size(ROI(1,r).TwoDSegmentationD.shell.(tissues{t}),1);
        newTotSize = size(ROI(1,r).TwoDSegmentationF.shell.(tissues{t}),1);
        
        voxChangeCombinedCond(r) = voxChangeCombinedCond(r) + newChangeSize * conductivities(t);
        totVox(r) = newTotSize + totVox(r);
    end
    voxChangeCombinedCond(r) = voxChangeCombinedCond(r) / totVox(r);
end

scatter(voxChangeCombinedCond, EChanged.noTissue)
xlabel('Avg Conductivity Changes Over All Voxels in Shell (S/m)')
ylabel('Magnitude of Field Change (V/m)')
title(['Avg Voxel Conductivity Change vs. Field Magnitude Change'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Vox_Tot_Change_Cond_vs_Field', name))
saveas(gcf,append(titleRun, '/Vox_Tot_Change_Cond_vs_Field', name, '.png'))
close all

voxChangeCombinedCond(isnan(voxAddCombinedCond)) = 0;
plot(ROI_radii, voxChangeCombinedCond)
xlabel('ROI radius')
ylabel('Avg Conductivity Changes Over All Voxels (S/m)')
title(['Avg Voxel Conductivity Change Over Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Change_Cond_Tot', name))
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Change_Cond_Tot', name, '.png'))
close all


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Avg Conductivity Change per Voxel Change

voxAddAvgCond = zeros(1, length(ROI));
for r = 1:length(ROI)
    voxSumAdd(r) = 0;
    for t = 1:length(tissues)
        
        %ROI shells

        newSize = size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1);
        
        voxAddAvgCond(r) = voxAddAvgCond(r) + newSize * conductivities(t);
        voxSumAdd(r) = newSize + voxSumAdd(r); % # of vox additions in ROI
    end
    voxAddAvgCond(r) = voxAddAvgCond(r) / voxSumAdd(r);
end

voxAddAvgCond(isnan(voxAddAvgCond)) = 0;
scatter(voxAddAvgCond, EChanged.noTissue)
xlabel('Conductivity Additions (S/m)')
ylabel('Magnitude of Field Change (V/m)')
title(['Avg Voxel Conductivity Addition vs. Field Magnitude Change'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Vox_Add_Cond_vs_Field', name))
saveas(gcf,append(titleRun, '/Vox_Add_Cond_vs_Field', name, '.png'))
close all

plot(ROI_radii, voxAddAvgCond)
xlabel('ROI radius')
ylabel('Conductivity Additions (S/m)')
title(['Avg Voxel Conductivity Addition Over Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Add_Cond', name))
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Add_Cond', name, '.png'))
close all

voxDelAvgCond = zeros(1, length(ROI));
for r = 1:length(ROI)
    voxSumDel(r) = 0;
    for t = 1:length(tissues)
        
        %ROI shells
        newSize = size(ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t}),1);
        
        voxDelAvgCond(r) = voxDelAvgCond(r) - newSize * conductivities(t);
        voxSumDel(r) = newSize + voxSumDel(r); % # of vox deletions in ROI
    end
    voxDelAvgCond(r) = voxDelAvgCond(r) / voxSumDel(r);
end

voxDelAvgCond(isnan(voxDelAvgCond)) = 0;
scatter(voxDelAvgCond, EChanged.noTissue)
xlabel('Conductivity Deletions (S/m)')
ylabel('Magnitude of Field Change (V/m)')
title(['Avg Voxel Conductivity Deletion vs. Field Magnitude Change'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Vox_Del_Cond_vs_Field', name))
saveas(gcf,append(titleRun, '/Vox_Del_Cond_vs_Field', name, '.png'))
close all

plot(ROI_radii, voxDelAvgCond)
xlabel('ROI radius')
ylabel('Conductivity Deletions (S/m)')
title(['Avg Voxel Conductivity Deletion Over Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Del_Cond', name))
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Del_Cond', name, '.png'))
close all

voxChangeAvgCond = zeros(1, length(ROI));
for r = 1:length(ROI)
    addCoord = [];
    delCoord = [];
    voxSumChange(r) = 0;
    numChange(r) = 0;
    
    for t = 1:length(tissues)
        
        %ROI shells
             
         if size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}), 1) > 0
            addCoord = [ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t})];
%             oldSize = size(ROI(1,r - 1).TwoDSegmentationD.(tissues{t}),1);
         end
  
         
         if size(ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t}), 1) > 0
            delCoord = [ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t})];
         end
%         

%         newSize = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1) - oldSize;
        
%         voxChangeCombinedCond(r) = voxChangeCombinedCond(r) + newSize * conductivities(t);
         voxSumChange(r) = conductivities(t) * size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1)...
             + voxSumChange(r);
         voxSumChange(r) = -conductivities(t) * size(ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t}),1)...
             + voxSumChange(r);
%          numChange(r) = size(ROI(1,r).TwoDSegmentationD.shell.(tissues{t}),1) + numChange(r);
         
    end
    
    %Number of additions + not-included deletions = total changes:
    if size(addCoord,1) == 0 || size(delCoord,1) == 0
             numChange(r) = size(addCoord, 1) + size(delCoord, 1);
    else
        numChange(r) = size(addCoord, 1) + sum(~ismember(delCoord, addCoord, 'rows'));
    end
%     voxChangeAvgCond(r) = (voxAddAvgCond(r) * voxSumAdd(r) + voxAddAvgCond(r) * voxSumDel(r)) / voxSumChange(r);
     voxChangeAvgCond(r) = voxSumChange(r) / numChange(r);

end

voxChangeAvgCond(isnan(voxChangeAvgCond)) = 0;
scatter(voxChangeAvgCond, EChanged.noTissue)
xlabel('Conductivity Change Avg per Voxel Changed (S/m)')
ylabel('Magnitude of Field Change (V/m)')
title(['Avg Voxel Conductivity Change vs. Field Magnitude Change'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Vox_Change_Cond_vs_Field', name))
saveas(gcf,append(titleRun, '/Vox_Change_Cond_vs_Field', name, '.png'))
close all

plot(ROI_radii, voxChangeAvgCond)
xlabel('ROI radius')
ylabel('Conductivity Change Avg per Voxel Changed (S/m)')
title(['Avg Voxel Conductivity Change Over Distance'])
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Change_Cond', name))
saveas(gcf,append(titleRun, '/Radius_vs_Vox_Change_Cond', name, '.png'))
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating unweighted & weighted Dice Coefficient and visualzizing field
% for ROI spheres


for t = 1:length(tissues)
    diceCoeff.(tissues{t}) = [];
    for r = 1:length(ROI)
        %If no tissue or negative dice ---> 0
        dcNumerator = max(0, 2*(size(ROI(1,r).TwoDSegmentationF.(tissues{t}),1)...
            - size(ROI(1,r).TwoDSegmentationDAdd.(tissues{t}),1)));
        
        dcDenominator = (size(ROI(1,r).TwoDSegmentationF.(tissues{t}),1)...
            + size(ROI(1,r).TwoDSegmentationI.(tissues{t}),1));
        
        diceCoeff.(tissues{t})(r) = dcNumerator / dcDenominator;    
    end
    diceCoeff.(tissues{t})(isnan(diceCoeff.(tissues{t}))) = 0;
    plot(ROI_radii, diceCoeff.(tissues{t}))
    xlabel('ROI Radius (vox)')
    ylabel('Dice Coefficient of Tissue')
    title(['Dice Coeff of Tissue over Distance in Spheres '])
    set(gca,'FontSize',14)
    hold on
end
legend(tissues);
saveas(gcf,append(titleRun, '/Radius_vs_Dice_Sphere', name))
saveas(gcf,append(titleRun, '/Radius_vs_Dice_Sphere', name, '.png'))


figure;
for t = 1:length(tissues)
    EChanged.(tissues{t})(isnan(EChanged.(tissues{t}))) = 0;
    dcPholder = diceCoeff.(tissues{t});
    [diceCoeff.(tissues{t}), sortIdx] = sort(diceCoeff.(tissues{t}),'ascend');
    EChanged.(tissues{t}) = EChanged.(tissues{t})(sortIdx);
    scatter(diceCoeff.(tissues{t}), EChanged.(tissues{t}))
    xlabel('Dice Coefficient of Tissue in Spheres')
    ylabel('Magnitude of Field Change (V/m)')
    title(['Dice Coeff of Tissue in Spheres vs Field Change'])
    hold on
    diceCoeff.(tissues{t}) = dcPholder;
end
legend(tissues);
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Dice_Spheres_vs_Field_Mag_Change', name))
saveas(gcf,append(titleRun, '/Dice_Spheres_vs_Field_Mag_Change', name, '.png'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Weighted Dice Coefficient for ROI spheres and visualzizing with it
figure;
for t = 1:length(tissues)
    diceCoeffWeighted.(tissues{t}) = [];
    for r = 1:length(ROI)
        diceCoeffWeighted.(tissues{t})(r) = conductivities(t) * diceCoeff.(tissues{t})(r) / (ROI_radii(:,r))^(2);
    end
    diceCoeffWeighted.(tissues{t})(isnan(diceCoeffWeighted.(tissues{t}))) = 0;
    
    plot(diceCoeffWeighted.(tissues{t}), EChanged.(tissues{t}))
    xlabel('Weighted Dice Coefficient of Tissue in ROI Spheres')
    ylabel('Magnitude of Field Change (N/C)')
    title(['Weighted Dice Coeff in ROI Spheres vs. Field Change'])
    hold on
end
legend(tissues);
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Weighted_Dice_Spheres_vs_Field_Mag_Change', name))
saveas(gcf,append(titleRun, '/Weighted_Dice_Spheres_vs_Field_Mag_Change', name, '.png'))
close all
uwPlot = [];
wPlot = [];
for t = 1:length(tissues)
    uwCorr.(tissues{t}) = [];
    weightedCorr.(tissues{t}) = [];
    uwCorr.(tissues{t}) = corrcoef(diceCoeff.(tissues{t}), EChanged.(tissues{t}));
    weightedCorr.(tissues{t}) = corrcoef(diceCoeffWeighted.(tissues{t}), EChanged.(tissues{t}));
    
    uwPlot = [uwPlot uwCorr.(tissues{t})(1,2)];
    wPlot = [wPlot  weightedCorr.(tissues{t})(1,2)];
end
figure;
scatter(uwPlot, wPlot)
xlabel('Correlation of Unweighted Dice Coefficient of Tissue With Field Change')
ylabel('Correlation of Weighted Dice Coefficient of Tissue With Field Change')
title(['UnWeighted vs. Weighted Dice Correlation Coeff With Field Change'])
text(uwPlot(1) + 0.01,wPlot(1)-0.05,'wm')
text(uwPlot(2) + 0.01,wPlot(2)-0.05,'gm')
text(uwPlot(3) + 0.01,wPlot(3)-0.05,'csf')
text(uwPlot(4) + 0.01,wPlot(4)-0.05,'bone')
text(uwPlot(5) + 0.01,wPlot(5)-0.05,'skin')
text(uwPlot(6) + 0.01,wPlot(6)-0.05,'eyes')
if (length(uwPlot) == 7)
    text(uwPlot(7) + 0.01,wPlot(7)-0.05,'air')
end
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/UnWeighted_Dice_vs_Weighted_Sphere', name))
saveas(gcf,append(titleRun, '/UnWeighted_Dice_vs_Weighted_Sphere', name, '.png'))
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating unweighted & weighted Dice Coefficient and visualzizing field
% for ROI shells


for t = 1:length(tissues)
    diceCoeff.(tissues{t}) = [];
    for r = 1:length(ROI)
        %If no tissue or negative dice ---> 0
        dcNumerator = max(0, 2*(size(ROI(1,r).TwoDSegmentationF.shell.(tissues{t}),1)...
            - size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1)));
        
        dcDenominator = (size(ROI(1,r).TwoDSegmentationF.shell.(tissues{t}),1)...
            + size(ROI(1,r).TwoDSegmentationI.shell.(tissues{t}),1));
        
        diceCoeff.(tissues{t})(r) = dcNumerator / dcDenominator;    
    end
    diceCoeff.(tissues{t})(isnan(diceCoeff.(tissues{t}))) = 0;
    plot(ROI_radii, diceCoeff.(tissues{t}))
    xlabel('ROI (sphere)Radius (vox)')
    ylabel('Dice Coefficient of Tissue')
    title(['Dice Coeff of Tissue over Distance'])
    set(gca,'FontSize',14)
    hold on
end
legend(tissues);
saveas(gcf,append(titleRun, '/Radius_vs_Dice_Shell', name))
saveas(gcf,append(titleRun, '/Radius_vs_Dice_Shell', name, '.png'))


for t = 1:length(tissues)
    diceCoeff.shell.(tissues{t}) = [];
    for r = 1:length(ROI)
        %If no tissue or negative dice ---> 0
        
        dcNumerator = max(0, 2*(size(ROI(1,r).TwoDSegmentationF.shell.(tissues{t}),1)...
            - size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1)));
        
        dcDenominator = (size(ROI(1,r).TwoDSegmentationF.shell.(tissues{t}),1)...
            + size(ROI(1,r).TwoDSegmentationI.shell.(tissues{t}),1));
        
        diceCoeff.shell.(tissues{t})(r) = dcNumerator / dcDenominator;
    end
    diceCoeff.shell.(tissues{t})(isnan(diceCoeff.shell.(tissues{t}))) = 0;
    plot(ROI_radii, diceCoeff.shell.(tissues{t}))
    xlabel('ROI Radius (vox)')
    ylabel('Dice Coefficient of Tissue in ROI Shells')
    title(['Dice Coeff of Tissue over Distance (Shells)'])
    set(gca,'FontSize',14)
    hold on
end
legend(tissues);
saveas(gcf,append(titleRun, '/Radius_vs_Dice_Shell', name))
saveas(gcf,append(titleRun, '/Radius_vs_Dice_Shell', name, '.png'))

figure;
for t = 1:length(tissues)
    EChanged.(tissues{t})(isnan(EChanged.(tissues{t}))) = 0;
    dcPholder = diceCoeff.(tissues{t});
    [diceCoeff.(tissues{t}), sortIdx] = sort(diceCoeff.(tissues{t}),'ascend');
    EChanged.(tissues{t}) = EChanged.(tissues{t})(sortIdx);
    scatter(diceCoeff.(tissues{t}), EChanged.(tissues{t}))
    xlabel('Dice Coefficient of Tissue')
    ylabel('Magnitude of Field Change (V/m)')
    title(['Dice Coeff of Tissue in ROI Shells vs Field Change'])
    hold on
    diceCoeff.(tissues{t}) = dcPholder;
end
legend(tissues);
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Dice_Shells_vs_Field_Mag_Change', name))
saveas(gcf,append(titleRun, '/Dice_Shells_vs_Field_Mag_Change', name, '.png'))
close all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Weighted Dice Coefficient and visualzizing with it
figure;
for t = 1:length(tissues)
    diceCoeffWeighted.(tissues{t}) = [];
    for r = 1:length(ROI)
        diceCoeffWeighted.(tissues{t})(r) = conductivities(t) * diceCoeff.(tissues{t})(r) / (ROI_radii(:,r))^(2);
    end
    diceCoeffWeighted.(tissues{t})(isnan(diceCoeffWeighted.(tissues{t}))) = 0;
    
    plot(diceCoeffWeighted.(tissues{t}), EChanged.(tissues{t}))
    xlabel('Weighted Dice Coefficient of Tissue')
    ylabel('Magnitude of Field Change (N/C)')
    title(['Weighted Dice Coeff of ROI Shells vs. Field Change'])
    hold on
end
legend(tissues);
set(gca,'FontSize',14)
saveas(gcf,append(titleRun, '/Weighted_Dice_Shells_vs_Field_Mag_Change', name))
saveas(gcf,append(titleRun, '/Weighted_Dice_Shells_vs_Field_Mag_Change', name, '.png'))
close all
uwPlot = [];
wPlot = [];
for t = 1:length(tissues)
    uwCorr.(tissues{t}) = [];
    weightedCorr.(tissues{t}) = [];
    uwCorr.(tissues{t}) = corrcoef(diceCoeff.(tissues{t}), EChanged.(tissues{t}));
    weightedCorr.(tissues{t}) = corrcoef(diceCoeffWeighted.(tissues{t}), EChanged.(tissues{t}));
    
    uwPlot = [uwPlot uwCorr.(tissues{t})(1,2)];
    wPlot = [wPlot  weightedCorr.(tissues{t})(1,2)];
end
figure;
scatter(uwPlot, wPlot)
xlabel('Correlation of Unweighted Dice Coefficient of Tissue With Field Change')
ylabel('Correlation of Weighted Dice Coefficient of Tissue With Field Change')
title(['UnWeighted vs. Weighted Dice Correlation Coeff With Field Change for Shells'])
text(uwPlot(1) + 0.01,wPlot(1)-0.05,'wm')
text(uwPlot(2) + 0.01,wPlot(2)-0.05,'gm')
text(uwPlot(3) + 0.01,wPlot(3)-0.05,'csf')
text(uwPlot(4) + 0.01,wPlot(4)-0.05,'bone')
text(uwPlot(5) + 0.01,wPlot(5)-0.05,'skin')
text(uwPlot(6) + 0.01,wPlot(6)-0.05,'eyes')
if (length(uwPlot) == 7)
    text(uwPlot(7) + 0.01,wPlot(7)-0.05,'air')
end
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/UnWeighted_Dice_vs_Weighted_Shells', name))
saveas(gcf,append(titleRun, '/UnWeighted_Dice_vs_Weighted_Shells', name, '.png'))
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Distance-Weighted Conductivity avg by ROI and visualzizing with field
%change
figure;
condChangeWeighted = [];
for r = 1:length(ROI)
    condChangeWeighted(r) = 0;
    
    for t = 1:length(tissues)
    condChangeWeighted(r) = condChangeWeighted(r) ...
        + size(ROI(1,r).TwoDSegmentationDAdd.shell.(tissues{t}),1) * conductivities(t)/ (ROI_radii(:,r))^(2)...
        - size(ROI(1,r).TwoDSegmentationDDel.shell.(tissues{t}),1) * conductivities(t)/ (ROI_radii(:,r))^(2);
    end    

end
scatter(condChangeWeighted, EChanged.noTissue)
xlabel('Distance-Weighted Conductivity')
ylabel('Magnitude of Field Change (N/C)')
title(['Distance-Weighted Conductivity of ROI Shells vs. Field Change'])
hold on
set(gca,'FontSize',11)
saveas(gcf,append(titleRun, '/Weighted_Cond_Change_Shells_vs_Field_Mag_Change', name))
saveas(gcf,append(titleRun, '/Weighted_Cond_Change_Shells_vs_Field_Mag_Change', name, '.png'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Distance-Weighted Conductivity and visualzizing with field
%change
vox_weights = [];
all_coords = [];
for t = 1:length(tissues)
    coords_tissue = ROI(1,length(ROI)).TwoDSegmentationD.(tissues{t});
    numVox = length(coords_tissue);
    ROI_centers_repeated = repmat(electrodes(1,:), numVox, 1);
    distances_vox = vecnorm(ROI_centers_repeated - coords_tissue, 2, 2);
    if isempty(distances_vox)
        continue
    else
        new_weight = conductivities(t)./(distances_vox).^(2);
        vox_weights = [vox_weights new_weight'];
        all_coords = [all_coords; coords_tissue];
    end
end



scatter3(all_coords(:,1),all_coords(:,2),all_coords(:,3),[],vox_weights);
colorbar
saveas(gcf,append(titleRun, '/Weighted_Cond_Change_3D', name))
saveas(gcf,append(titleRun, '/Weighted_Cond_Change_3D', name, '.png'))

slice_origin = [0, 0, 38.9879];
slice_select = round(all_coords(:,3), 4) == slice_origin(:,3);
slice_coords = all_coords(slice_select,:);
slice_weights = vox_weights(slice_select);
scatter(slice_coords(:,1), slice_coords(:,2), [], slice_weights, 'filled');
colorbar
saveas(gcf,append(titleRun, '/Weighted_Cond_Change_2D', name))
saveas(gcf,append(titleRun, '/Weighted_Cond_Change_2D', name, '.png'))
end
