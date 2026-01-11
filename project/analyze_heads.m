

tissues = {'wm','gm','csf','bone','skin','eyes','air'};

% subjects = append('TIME',{'020','041','045','052', '057','072','081','082','085'});
subjects = append('TIME',{'045','052', '057','072','081','082','085'});
just_masks = 1;
TIME=1;
is_DE_run=1;

if(~TIME)
    directory = '/Volumes/spiral/BSSLab/TMSTestRetest/Models/'; 

    % List contents of the directory
    contents = dir(directory);

    % Extract folder names using logical indexing
    folder_names = {contents([contents.isdir]).name};

    % Remove '.' and '..' from the list
    subjects = folder_names(~ismember(folder_names, {'.', '..'}));

end
    
ROI_radii = logspace(0,250,20);
for i =1:length(subjects)

    if (~just_masks)
        run1 = ['/Volumes/spiral/BSSLab/TIME/Modeling/' subjects{i} '/Model/headreco_T1/run1/m2m_' subjects{i} '/mask_prep/'];
        run2 = ['/Volumes/spiral/BSSLab/TIME/Modeling/' subjects{i} '/Model/headreco_T1/m2m_' subjects{i} '/mask_prep/'];
        post = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/' subjects{i} '_run2/' subjects{i} '_run2_c1_result_map-r2.mat']);
        preFit = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/' subjects{i} '_run1/' subjects{i} '_run1_c1_result_map-r2.mat']);
        pre = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/' subjects{i} '_run1/' subjects{i} '_run1_c1_result.mat']);
        geoPost = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/' subjects{i} '_run2/' subjects{i} '_run2_model.mat']);
        geoPreFit = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/' subjects{i} '_run1/' subjects{i} '_run1_model_map-r2.mat']);
        geoPre = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/' subjects{i} '_run1/' subjects{i} '_run1_model.mat']);
        configPath = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/setup/configs.mat']);
        configs = load(configPath).configs;
        try
            electrodes = configs{1}(1,:);
        catch
            warning("config load didn't work");
        end

    else

    if(~TIME)
        subject_dir = join([directory subjects{i}]);

        % List contents of the directory
        contents = dir(fullfile(subject_dir, 'run*'));

        % Extract folder names
        folder_names = {contents.name};

        % Extract numbers from folder names
        folder_numbers = cellfun(@(x) str2double(regexp(x, '\d+', 'match')), folder_names);

        % Find unique folder numbers
        unique_numbers = unique(folder_numbers);

        if numel(unique_numbers) < 2
            fprintf('There are not enough differently named folders with the pattern "run" followed by a number.\n');
            continue
        else
            % Find lowest and highest numbers
            lowest_number = min(unique_numbers);
            highest_number = max(unique_numbers);

            % Find corresponding folder names
            lowest_folder = sprintf('run%d', lowest_number);
            highest_folder = sprintf('run%d', highest_number);

            % Display the results
            fprintf('Lowest number folder: %s\n', lowest_folder);
            fprintf('Highest number folder: %s\n', highest_folder);
        end

        run1 = join([subject_dir '/' lowest_folder '/m2m_' subjects{i} '/mask_prep/']);
        if is_DE_run
            run2 = append(subjects{i}, '_VoxCount/', 'dilated/');
        else
            run2 = join([subject_dir '/' highest_folder '/m2m_' subjects{i} '/mask_prep/']);
        end
        geoPost = join([subject_dir '/' highest_folder '/m2m_' subjects{i} '/' subjects{i} '.geo']);
        

        if ~exist(run2, 'dir')
            fprintf('No run final dir found, trying w/o m2m \n');
            run2 = join([subject_dir '/' highest_folder '/mask_prep/']);
            geoPost = join([subject_dir '/' highest_folder '/' subjects{i} '.geo']);
        end

        if ~exist(run2, 'dir')
            fprintf('No run final dir found, skipping \n');
            continue
        end

        if ~exist(run1, 'dir')
            fprintf('No run initial dir found, trying w/o m2m \n');
            run1 = join([subject_dir '/' lowest_folder '/mask_prep/']);
        end

        if ~exist(run1, 'dir')
            fprintf('No run initial dir found, skipping \n');
            continue
        end

    else
        
        run1 = ['/Volumes/spiral/BSSLab/TIME/Modeling/' subjects{i} '/Model/headreco_T1/run1/m2m_' subjects{i} '/mask_prep/'];
        if is_DE_run
            run2 = append(subjects{i}, '_VoxCount/', 'dilated/');
        else
            run2 = ['/Volumes/spiral/BSSLab/TIME/Modeling/' subjects{i} '/Model/headreco_T1/m2m_' subjects{i} '/mask_prep/'];
        end
            geoPost = join(['/Volumes/spiral/BSSLab/MEQ/' subjects{i} '/' subjects{i} '_run2/' subjects{i} '_run2_model.mat']);

    end
        electrodes = [];
        preFit = {};
        post = {};
    end

    [EAngle, EMagnitude, ROI, Segmentation] = analyze_head(electrodes, electrodes, geoPost, preFit, post, run1, run2, ROI_radii, string(subjects(i)));
    if (~just_masks)
        EPost = load(post).E_vector;
        EPre = load(pre).E_vector;
        EPreFit = load(preFit).E_vector;
        GeometryPost = load(geoPost);
        GeometryPre = load(geoPre);
        GeometryPreFit = load(geoPreFit);


        brain_blocking(subjects(i), GeometryPost, GeometryPre, GeometryPreFit, ROI, configs, EPost, EPreFit, EPre);
    else
        sumI = zeros(size(Segmentation(1).initial.csf));
        sumF = zeros(size(Segmentation(1).final.csf));
        sumD = zeros(size(Segmentation(1).diff.csf));

        for t = 1:length(tissues)
            if(isempty(Segmentation(1).initial.(tissues{t})))
                 warning(join(['Missing Tissue: ' tissues{t}]));
                continue
            end
            sumI = sumI + Segmentation(1).initial.(tissues{t}) * t;
            sumF = sumF + Segmentation(1).final.(tissues{t}) * t;
            sumD = sumD + Segmentation(1).diff.(tissues{t}) * t;

        end

        % By tissue index
        tiledlayout(3,1)
        figure('units','normalized', 'OuterPosition', [0 0 1 1])

        nexttile
        imagesc(sumI(:,:,160))
        bottom = min(min(sumI(:)), min(sumF(:)));
        top = max(max(sumI(:)), max(sumF(:)));
        colorbar;
        caxis manual
        caxis([bottom top])
        title(['Initial mask, min=' num2str(min(sumI(:))) ', max=' num2str(max(sumI(:)))])
        set(gca,'FontSize',14)


        nexttile
        imagesc(sumF(:,:,160)) %50 for BI
        colorbar;
        caxis manual
        caxis([bottom top])
        title(['Final mask, min=' num2str(min(sumF(:))) ', max=' num2str(max(sumF(:)))])
        set(gca,'FontSize',14)

        nexttile
        imagesc(sumD(:,:,160))
        c=colorbar;
        title(['Change in masks, min=' num2str(min(sumD(:))) ', max=' num2str(max(sumD(:)))])
        set(gca,'FontSize',14)

        voxcount_dir = join([subjects{i} '_VoxCount']);
        if ~exist(voxcount_dir, 'dir')
            mkdir(voxcount_dir);
        end
        
        if is_DE_run
            saveas(gcf,append(voxcount_dir, '/Dilate_Axial_View'))
            saveas(gcf,append(voxcount_dir, '/Dilate_Axial_View', '.png'))
        else
            saveas(gcf,append(voxcount_dir, '/Axial_View'))
            saveas(gcf,append(voxcount_dir, '/Axial_View', '.png'))
        end

        close all

         % By tissue index
        tiledlayout(3,1)
        figure('units','normalized', 'OuterPosition', [0 0 1 1])

        nexttile
        imagesc(rot90(squeeze(sumI(:,120,:))))
        bottom = min(min(sumI(:)), min(sumF(:)));
        top = max(max(sumI(:)), max(sumF(:)));
        colorbar;
        caxis manual
        caxis([bottom top])
        title(['Initial mask, min=' num2str(min(sumI(:))) ', max=' num2str(max(sumI(:)))])
        set(gca,'FontSize',14)


        nexttile
        imagesc(rot90(squeeze(sumF(:,120,:))))
        colorbar;
        caxis manual
        caxis([bottom top])
        title(['Final mask, min=' num2str(min(sumF(:))) ', max=' num2str(max(sumF(:)))])% c.Ticks=0:12;
        set(gca,'FontSize',14)

        nexttile
        imagesc(rot90(squeeze(sumD(:,120,:))))
        c=colorbar;
        title(['Change in masks, min=' num2str(min(sumD(:))) ', max=' num2str(max(sumD(:)))])
        set(gca,'FontSize',14)
        if is_DE_run
            saveas(gcf,append(voxcount_dir, '/Dilate_Coronal_View'))
            saveas(gcf,append(voxcount_dir, '/Dilate_Coronal_View', '.png'))
        else
            saveas(gcf,append(voxcount_dir, '/Coronal_View'))
            saveas(gcf,append(voxcount_dir, '/Coronal_View', '.png'))
        end
        close all
        if is_DE_run
            save(join([voxcount_dir '/Dilate_Segmentation.mat']), 'Segmentation', '-v7.3');
        else
            save(join([voxcount_dir '/Segmentation.mat']), 'Segmentation', '-v7.3');
        end

    end

end
