% %Highest level field and voxel analysis program
%
% %load voxel differences:
% [headDiff, files] = compare_heads('/Volumes/spiral/BSSLab/MNL/edited/017/run1/m2m_017/mask_prep', '/Volumes/spiral/BSSLab/MNL/edited/017/run2/m2m_017/mask_prep');
%
% %load the head model
% %TO DO: Specify the variable that you're loading into
% load('/Volumes/spiral/BSSLab/MNL/Simulations/017/017pre.mat');
%
% load('/Volumes/spiral/BSSLab/MNL/Simulations/017/017post.mat', 'Geometry');
% load('/Volumes/spiral/BSSLab/MNL/edited/016/m2m_016/skin.mat');
% geometrySTL = Geometry;
% import_msh('/Volumes/spiral/BSSLab/MNL/edited/016/m2m_016/016.msh', 'v');
% element_centers = (Geometry.node(Geometry.cell(:,1),:) + Geometry.node(Geometry.cell(:,2),:) + Geometry.node(Geometry.cell(:,3),:) + Geometry.node(Geometry.cell(:,4),:))/4;
% ROI(1).center = [0 75 60];
% ROI(2).center = [0 -50 -120];
% ROI(3).center = [-83 -28 23];
% ROI(4).center = [83 -28 23];
% ROI_radius = 30;
% ROI_labels = categorical({'Bone','GM', 'WM', 'CSF', 'Eyes', 'Skin'});
% Electrode_labels = categorical({'Front','Back', 'Left', 'Right'});
%
% len_elem = length(element_centers);
% %Fill arrays with zeros for now:
% %Will be filled with field and distance data in loop
% E_ROI_byTissue = zeros(6, 4);
% E_ROI_noTissue = zeros(1, 4);
% E_Percentage = zeros(6, 4);
% dist_elem_ROI = zeros(len_elem, 4);
%
% %Loop through each electrode configuration and collect field data
% for r = 1:4
%
%     %Function to select the voxel data within the ROI
%     %ROI(r).select = ROI_select_heads('/Volumes/spiral/BSSLab/MNL/edited/017/run1/m2m_017/mask_prep', '/Volumes/spiral/BSSLab/MNL/edited/017/run2/m2m_017/mask_prep', Geometry, 40, ROI(r).center);
%    %Load appropriate field data
%     if (r == 1 || r == 2)
%         load('/Volumes/spiral/BSSLab/MNL/simulations/017/017pre_AP_result_mapped.mat')
%         Epre_AP = E;
%         load('/Volumes/spiral/BSSLab/MNL/simulations/017/017post_AP_result_mapped.mat')
%         Epost_AP = E;
%         E = Epost_AP - Epre_AP;
%         Emag_AP = vecnorm(Epost_AP, 2, 2) - vecnorm(Epre_AP, 2, 2);
%         Eangle_AP = [];
%         Eangle_AP = acosd(dot(Epost_AP,Epre_AP, 2)./(vecnorm(Epost_AP, 2, 2).*vecnorm(Epre_AP, 2, 2)));
%
%
%         Etot = Epost_AP;
%     end
%     if (r == 3 || r == 4)
%         load('/Volumes/spiral/BSSLab/MNL/simulations/017/017pre_LR_result_mapped.mat')
%         Epre_LR = E;
%         load('/Volumes/spiral/BSSLab/MNL/simulations/017/017post_LR_result_mapped.mat')
%         Epost_LR = E;
%         E = Epost_LR - Epre_LR;
%         Etot = Epost_LR;
%         Emag_LR = vecnorm(Epost_LR, 2, 2) - vecnorm(Epre_LR, 2, 2);
%         Eangle_LR = [];
%         Eangle_LR = acosd(dot(Epost_LR,Epre_LR, 2)./(vecnorm(Epost_LR, 2, 2).*vecnorm(Epre_LR, 2, 2)));
%     end
%
%     %Collect data that is non-specific to tissue
%     ROI_centers_rep = repmat(ROI(r).center, len_elem, 1);
%     dist_elem_ROI_temp = element_centers - ROI_centers_rep; % vector from each elem's center to the ROI center % vector from each elem's center to the ROI center
%     dist_elem_ROI(:, r) = vecnorm(dist_elem_ROI_temp,2,2); % absolute distance from each elem's center to the ROI center
%     ROI(r).flag = dist_elem_ROI(:, r) < ROI_radius;
%     ROI(r).noTissue = ROI(r).flag;
%     E_ROI_noTissue(r) = mean(E(ROI(r).noTissue));
%
%     %Loop through each tissue type
%     for f = 1:6
%
%         ROI(r).flags(:, f) = ROI(r).flag & Geometry.field == f;
%
%         E_ROI_byTissue(f, r) = mean(E(ROI(r).flags(:, f)));
%         E_ROI_byTissue(isnan(E_ROI_byTissue))=0;
%         E_Percentage(f, r) =  E_ROI_byTissue(f, r)/mean(Etot(ROI(r).flags(:, f)));
%     end
%
% end
%
%
% % E_ROI_byTissue = E_ROI_byTissue';
% noAirDiff = headDiff(2:7)';
% totVox = headDiff(:,2);
% percentChange = noAirDiff ./ totVox(2:7);
%
% %Visualization junkyard:
%
% %  for r = 1:4
% %
% % figure(1);
% %
% % imagesc(ROI(r).select(:, f));
% %
% % hold on;
% % Efield = E(ROI(r).flag);
% %
% % imagesc(Efield);
% %
% % hold on;
% %
% % imagesc(isoMask('/Volumes/spiral/BSSLab/MNL/edited/017/run1/m2m_017/mask_prep/MASK_SKIN.nii.gz'));
% %
% %
% %  end
%
% set(gca,'xticklabel',ROI_labels)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%COMPARE E FIELD ANGLE CHANGE TO CSF VOX CHANGE
avgAngle = mean(EAngleOut(1,:)');
avgMag = mean(EMagOut(1,:)');
load(pre)
Epre = E;
load(post)
Epost = E;
percentChange = mean(vecnorm(Epost,2,2)./vecnorm(Epre,2,2));
bar(Segmentation.voxChangeSum.csf, avgAngle)
xlabel('CSF Voxel Change')
ylabel('Average Angle Change (degrees)')
set(gca,'FontSize',14)



% bar3(noAirDiff, E_ROI_byTissue, 'grouped');
% bar3(percentChange, E_ROI_byTissue, 'stacked');
%
% figure(3)
% bar(Electrode_labels, E_ROI_noTissue)
%
% %FROM ROI_Select code%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
geometry = Geometry;
element_centers = (geometry.node(geometry.cell(:,1),:) + geometry.node(geometry.cell(:,2),:) + geometry.node(geometry.cell(:,3),:) + geometry.node(geometry.cell(:,4),:))/4;
%
% ROI_centers_rep = repmat(electCenter, length(element_centers), 1);
%
% flagROI = vecnorm(element_centers - ROI_centers_rep, 2, 2) < radius;
% elem_trim = element_centers(flagROI, :);
% E = E(flagROI, :);
%
% path1 = '/Volumes/spiral/BSSLab/MNL/edited/017/run1/m2m_017/mask_prep';
% path2 = '/Volumes/spiral/BSSLab/MNL/edited/017/run2/m2m_017/mask_prep';
% radius = 20;
% [isoTotal,filenames1, filenames2] = ROI_select_heads(path1, path2, geometry, electCenter, ROI_radius);
% [conversion20 conversion2] = ROI_select(skinMask1, skinMask2, skinMask, geometry, electCenter, radius);
% conversion50 = ROI_select(skinMask1, skinMask2, skinMask, geometry, 50, electCenter);
%
%
%
%
% changeCSF = length(isoTotal.csf);
% changeFrontROI = E_ROI_noTissue(1);
% bar(changeCSF, changeFrontROI)
%
load(geometry)
nodeDownSample1 = Geometry.node([1:10:length(Geometry.node)], :);
conversion2 = ROI(1, 10).TwoDSegmentationD.skin;
% nodeDownSample2 = geometrySTL.node([1:10:length(geometrySTL.node)], :);
convDownSample = conversion2([1:10:length(conversion2)], :);
% convDownSample20 = conversion20([1:10:length(conversion20)], :);
% convDownSample50 = conversion50([1:10:length(conversion50)], :);
%
Eangle = EAngleOut(1,:)';
Emag = EMagOut(1,:)';
eDownSample = E([1:100:length(E)], :);
eMagDownSample = Emag([1:100:length(Eangle)], :);
eAngleDownSample = Eangle([1:100:length(Eangle)], :);
%
scatter3(E(:, 1), E(:,2), E(:, 3), EAngleOut(1,:)')
% hold on
scatter3(eDownSample(:, 1), eDownSample(:,2), eDownSample(:, 3), eMagDownSample)
scatter3(eDownSample(:, 1), eDownSample(:,2), eDownSample(:, 3), eAngleDownSample)
hold on
%
% [conversion20 conversion2] = ROI_select(skinMask1, skinMask2, skinMask, geometry, electCenter, radius);
plot3(conversion2(:,1), conversion2(:,2), conversion2(:,3), 'b.'); %or whatever marker you want
hold on
% plot3(conversion(:,1), conversion(:,2), conversion(:,3), 'r.'); %or whatever marker you want
% hold on
% plot3(elem_trim(:,1) + E(:,1), elem_trim(:,2) + E(:,2), elem_trim(:,3) + E(:,3), 'b.');
% hold on
% plot3(geometry.node(:,1), geometry.node(:,2), geometry.node(:,3), 'g.');
% plot3(geometrySTL.node(:,1), geometrySTL.node(:,2), geometrySTL.node(:,3), 'g.');
% hold on
% plot3(geo1.Geometry.node(:,1), geo1.Geometry.node(:,2), geo1.Geometry.node(:,3), 'r.');
plot3(nodeDownSample1(:,1), nodeDownSample1(:,2), nodeDownSample1(:,3), 'r.');
% plot3(nodeDownSample2(:,1), nodeDownSample2(:,2), nodeDownSample2(:,3), 'b.');
hold on
pholder = ROI(1,5).TwoDSegmentationF.skin;
plot3(pholder(:,1), pholder(:,2), pholder(:,3))
hold on
plot3(Geometry.node(:,1), Geometry.node(:,2), Geometry.node(:,3), 'g.');
plot3(convDownSample(:,1), convDownSample(:,2), convDownSample(:,3), 'g.');
% plot3(convDownSample20(:,1), convDownSample20(:,2), convDownSample20(:,3), 'y.');
% plot3(convDownSample50(:,1), convDownSample50(:,2), convDownSample50(:,3), 'g.');
quiver3(element_centers(:,1),element_centers(:,2),element_centers(:,3),E(:, 1),E(:, 2),E(:, 3))
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
tissues = {'wm','gm','csf','bone','skin','eyes','air'};
% sumSel = zeros(size(ROI(1).select.skin));
% for t = 1:length(tissues)
%    sumSel = sumSel + ROI(1).select.(tissues{t}) * t;
% end
%
% sumRef = zeros(size(ROI(1).ref.skin));
% for t = 1:length(tissues)
%    tissues{t}
%    sumRef = sumRef + ROI(1).ref.(tissues{t}) * t;
% end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run Analyze function
for r = length(ROI_radii(:,1))
    [EAngleOut(r),EMagOut(r), ROI, Segmentation(r)] = Analyze(electrodes, '/Users/sidsimon/Documents/MATLAB/head17.mat', pre, post, headI, headF, ROI_radii(:,r));
    storedROIs(r) = ROI;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get voxel information for slice-by-slice viewing, divided by
%concuctivities or tissue type
sumI = zeros(size(Segmentation1TMS(1).initial.skin));
sumF = zeros(size(Segmentation1TMS(1).initial.skin));
sumD = zeros(size(Segmentation1TMS(1).initial.skin));

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

% 
% sumAdd = zeros(size(maskAddition));
% sumDel = zeros(size(maskDeletion));
% sumD = zeros(size(maskChange));
% for t = 1:length(tissues)
%     sumAdd = sumAdd + maskAddition * t;
%     sumDel = sumDel + maskDeletion * t;
%     sumD = sumD + maskChange * t;
%   imagesc(sumD(:,:,130))
% end
% for t = 1:length(tissues)
%
%    ROI(1).select.(tissues{t}) = ROI(1).ref.(tissues{t}) - ROI(1).select.(tissues{t});
% end

% c.TickLabels=tissues

tiledlayout(3,1)

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


%For conductivities:
tiledlayout(3,1)

nexttile
imagesc(condsumI(:,:,130))
bottom = min(min(condsumI(:)), min(condsumF(:)));
top = max(max(condsumI(:)), max(condsumF(:)));
colorbar;
caxis manual %Should I change to 0 -> 7
caxis([bottom top])
title(['Initial mask, min=' num2str(min(condsumI(:))) ', max=' num2str(max(condsumI(:)))])
set(gca,'FontSize',14)
% c.Ticks=0:12;


nexttile
imagesc(condsumF(:,:,130))
colorbar;
caxis manual
caxis([bottom top])
title(['Final mask, min=' num2str(min(condsumF(:))) ', max=' num2str(max(condsumF(:)))])% c.Ticks=0:12;
set(gca,'FontSize',14)

nexttile
imagesc(condsumD(:,:,130))
c=colorbar;
title(['Change in masks, min=' num2str(min(condsumD(:))) ', max=' num2str(max(condsumD(:)))])
set(gca,'FontSize',14)

voxChangeCombined = zeros(length(ROI),length(tissues));



for t = 1:length(tissues)
    oldSize = 0;
    for r = 1:length(ROI_radii)
        newSize = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1) - oldSize;
        voxChangeCombined(r, t) = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1);
        oldSize = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1);
    end
    plot(ROI_radii, voxChangeCombined(:,t))
    hold on
end
legend(tissues);
xlabel('ROI radius (voxels)')
ylabel('Number of voxel changes')
title(['Voxel Changes as Function of Distance'])
set(gca,'FontSize',14)
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get E ROIs per tissue type with logspace radii
Epre = load(pre);
Epre = Epre.E;
Epost = load(post);
Epost = Epost.E;

for r = length(ROI_radii(:,1))
    
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
            placeholderTissue = E(flag & Geometry.field == f);
            EChanged(r).(tissues{f})= [EChanged(r).(tissues{f}), mean(placeholderTissue)];
        end
        
        oldRad = rad;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Quantify voxel change by conductivity, finding average and vizualizeing
%with change of field


lenROI = length(ROI);
voxChangeCombinedCond = zeros(1, lenROI);
for r = 1:lenROI
    for t = 1:length(tissues)
        voxChangeCombinedCond(r) = voxChangeCombinedCond(r) + size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1) * conductivities(t);
    end
end



EChanged.noTissue(isnan(EChanged.noTissue)) = 0;
condChangeRep = repmat(voxChangeCombinedCond, lenROI, 1);
eChangeRep = repmat(EChanged.noTissue, lenROI, 1);
ratioCondE = condChangeRep' ./ eChangeRep;
imagesc(ratioCondE)
colorbar;
for r = 1:length(ROI)
    scatter(condChangeRep(r,:), eChangeRep(r,:))
    hold on
end
xlabel('ROI Conductivity Change (S/m)');
ylabel('Magnitude of Field Change (V/m)');
set(gca,'FontSize',14)
%  for r = lenROI
%      condChangeAndEChange(r,:) =
%      length(ROI_radii(:,1))
%  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get sum of vox change and visualize, along with field change

for t = 1:length(tissues)
    for r = 1:length(ROI)
        voxChangeCombined(r, t) = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1);
    end
    %     plot3(ROI_radii(2,:), voxChangeCombined(:,t)', EChanged.noTissue)
    %     hold on
end
xlabel('ROI Radius');
ylabel('Number of Voxels Changed');
zlabel('Magnitude of Field Change (V/m)');
legend(tissues);
set(gca,'FontSize',14)


% For 2D:
tiledlayout(3,1)
nexttile
for t = 1:length(tissues)
    oldSize = 0;
    for r = 1:length(ROI)
        newSize = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1);
        voxChangeCombined(r, t) = newSize - oldSize;
        oldSize = newSize;
    end
end
legend(tissues);
set(gca,'FontSize',14)

nexttile
for t = 1:length(tissues)
    plot(voxChangeCombined(:,t), EChanged.(tissues{t}))
    xlabel('Number of Voxels Changed')
    ylabel('Magnitude of Field Change (V/m)')
    title(['Voxel Changes vs. Field Magnitude Change'])
    hold on
end
legend(tissues);
set(gca,'FontSize',14)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Visualize field change by radii

nexttile
plot(ROI_radii, EChanged.noTissue)
xlabel('ROI Radius')
ylabel('Magnitude of Field Change (V/m)')
title(['Field Magnitude Change Over Distance'])
set(gca,'FontSize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Conductivity Changes from Voxel Changes
condsumD = condsumF - condsumI;
imagesc(condsumD(:, :, 130))
colorbar;
set(gca,'FontSize',14)

voxChangeCombinedCond = zeros(1, length(ROI));
for r = 1:length(ROI)
    voxSum = 0;
    for t = 1:length(tissues)
        
        %ROI shells
        if r - 1 > 0
            oldSize = size(ROI(1,r - 1).TwoDSegmentationF.(tissues{t}),1);
        else
            oldSize = 0;
        end
        
        newSize = size(ROI(1,r).TwoDSegmentationF.(tissues{t}),1) - oldSize;
        
        voxChangeCombinedCond(r) = voxChangeCombinedCond(r) + newSize * conductivities(t);
        voxSum = newSize + voxSum;
    end
    voxChangeCombinedCond(r) = voxChangeCombinedCond(r) / voxSum;
end

scatter(voxChangeCombinedCond, EChanged.noTissue)
xlabel('Conductivity Changes (S/m)')
ylabel('Magnitude of Field Change (V/m)')
title(['Avg Voxel Conductivity Change vs. Field Magnitude Change'])
set(gca,'FontSize',14)

plot(ROI_radii, voxChangeCombinedCond)
xlabel('ROI radius')
ylabel('Conductivity Changes (S/m)')
title(['Avg Voxel Conductivity Change Over Distance'])
set(gca,'FontSize',14)

sumTest = sumD;
sumTest(sumTest == 0) = 8;
sumTest(sumTest ~= 8) = 0;
sumTest(sumTest == 8) = 1;
sum(sumTest, 'all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Dice Coefficient and visualzizing with it


for t = 1:length(tissues)
    diceCoeff.(tissues{t}) = [];
    shellFOld = 0;
    shellIOld = 0;
    shellDOld = 0;
    for r = 1:length(ROI)
        shellFNew = size(ROI(1,r).TwoDSegmentationF.(tissues{t}),1) - shellFOld;
        shellINew = size(ROI(1,r).TwoDSegmentationI.(tissues{t}),1) - shellIOld;
        shellDNew = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1) - shellDOld;
        diceCoeff.(tissues{t})(r) = 2*(shellFNew - shellDNew) / (shellFNew + shellINew);
        shellFOld = size(ROI(1,r).TwoDSegmentationF.(tissues{t}),1);
        shellIOld = size(ROI(1,r).TwoDSegmentationI.(tissues{t}),1);
        shellDOld = size(ROI(1,r).TwoDSegmentationD.(tissues{t}),1);
    end
    diceCoeff.(tissues{t})(isnan(diceCoeff.(tissues{t}))) = 0;
    plot(ROI_radii, diceCoeff.(tissues{t}))
    xlabel('ROI Radius (vox)')
    ylabel('Dice Coefficient of Tissue')
    title(['Dice Coeff of Tissue over Distance'])
    set(gca,'FontSize',14)
    hold on
end
legend(tissues);

figure;
for t = 1:length(tissues)
    EChanged.(tissues{t})(isnan(EChanged.(tissues{t}))) = 0;
    dcPholder = diceCoeff.(tissues{t});
    [diceCoeff.(tissues{t}), sortIdx] = sort(diceCoeff.(tissues{t}),'ascend');
    EChanged.(tissues{t}) = EChanged.(tissues{t})(sortIdx);
    scatter(diceCoeff.(tissues{t}), EChanged.(tissues{t}))
    xlabel('Dice Coefficient of Tissue')
    ylabel('Magnitude of Field Change (V/m)')
    title(['Dice Coeff of Tissue vs Field Change'])
    hold on
    diceCoeff.(tissues{t}) = dcPholder;
end
legend(tissues);
set(gca,'FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating Weighted Dice Coefficient and visualzizing with it
figure;
for t = 1:length(tissues)
    diceCoeffWeighted.(tissues{t}) = [];
    for r = 1:length(ROI)
        diceCoeffWeighted.(tissues{t})(r) = conductivities(t) * diceCoeff.(tissues{t})(r) / (r * 50)^(2);
    end
    diceCoeffWeighted.(tissues{t})(isnan(diceCoeffWeighted.(tissues{t}))) = 0;
    
    plot(diceCoeffWeighted.(tissues{t}), EChanged.(tissues{t}))
    xlabel('Weighted Dice Coefficient of Tissue')
    ylabel('Magnitude of Field Change (N/C)')
    title(['Weighted Dice Coeff vs. Field Change'])
    hold on
end
legend(tissues);
set(gca,'FontSize',14)

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
xlabel('Unweighted Dice Coefficient of Tissue')
ylabel('Weighted Dice Coefficient of Tissue')
title(['UnWeighted vs. Weighted Dice Coeff'])
set(gca,'FontSize',14)


% %Looping for different radii:
% for i = 50:50:250
%    tissues{t}
%    sumRef = sumRef + ROI(1).ref.(tissues{t}) * t;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Try to scale all datapoints according to radius
EChanged.noTissue = [];
for t = 1:length(tissues)
    EChanged.(tissues{t}) = [];
end

electCen = electrodes(1,:);
centers_rep = repmat(electCen, length(Geometry.node), 1);
for rad = 1:1:250
    flagNode = vecnorm(Geometry.node - centers_rep, 2, 2) < rad;
    flag = flagNode(Geometry.cell(:,1)) & flagNode(Geometry.cell(:,2)) & flagNode(Geometry.cell(:,3));
    E = abs(vecnorm(Epost, 2,2) - vecnorm(Epre, 2,2));
    placeholder = E(flag,:);
    EChanged.noTissue= [EChanged.noTissue, mean(placeholder)];
    for f = 1:length(tissues)
        placeholderTissue = E(flag & Geometry.field == f);
        EChanged.(tissues{f})= [EChanged.(tissues{f}), mean(placeholderTissue)];
    end
end

