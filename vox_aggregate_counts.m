% List contents of the directory matching the pattern
contents = dir('*_VoxCount');

% Extract folder names using logical indexing
folder_names = {contents([contents.isdir]).name};

replace_dir = 1;

TIME_folder_names = folder_names(~cellfun(@isempty,strfind(folder_names,'TIME')));
BI_folder_names = folder_names(~cellfun(@isempty,strfind(folder_names,'BI')));

% Loop through each folder name
for i = 1:numel(BI_folder_names)
    folder_name = BI_folder_names{i};

    seg_file_path = fullfile(folder_name, 'Segmentation.mat');

    % Check if the Segmentation.mat file exists in the folder
    if exist(seg_file_path, 'file')
        % Load the Segmentation.mat file
        data = load(seg_file_path, 'Segmentation');
        voxChanges = data.Segmentation.voxChangeSum;

        tissues = fieldnames(voxChanges);

        subplot(4,3,i);

        tissues = fieldnames(voxChanges);

        voxChangeAccumulator = [];
        for j = 1:numel(tissues)
            voxChangeAccumulator = [voxChangeAccumulator; voxChanges.(tissues{j})];
        end

        bar(voxChangeAccumulator);
        xlabel('Tissue Type');
        ylabel('Number of Voxel Changes');
        set(gca,'xticklabel',tissues);
        vox_agg_dir = 'VoxCountAgg';
        if replace_dir && exist(vox_agg_dir, 'dir')
            rmdir(vox_agg_dir, 's');
            mkdir(vox_agg_dir);
        elseif ~exist(vox_agg_dir, 'dir')
            mkdir(vox_agg_dir);
        end

        title(strrep(folder_name, "_", " "));
    else
        disp(['Segmentation.mat not found in folder: ', folder_name]);
    end
end
saveas(gcf,append(vox_agg_dir, '/BI_Counts'))
close all


% Loop through each folder name
for i = 1:numel(TIME_folder_names)
    folder_name = TIME_folder_names{i};

    seg_file_path = fullfile(folder_name, 'Segmentation.mat');

    % Check if the Segmentation.mat file exists in the folder
    if exist(seg_file_path, 'file')
        % Load the Segmentation.mat file
        data = load(seg_file_path, 'Segmentation');
        voxChanges = data.Segmentation.voxChangeSum;
        totalCount.(TIME_folder_names{i}).voxChanges =  voxChanges;
        
        tissues = fieldnames(voxChanges);

        % subplot(3,2,i);
        subplot(3,3,i);

        tissues = fieldnames(voxChanges);

        voxChangeAccumulator = [];
        for j = 1:numel(tissues)
            totalCount.(TIME_folder_names{i}).voxInitial.(tissues{j}) = sum(data.Segmentation.initial.(tissues{j}),"all");
            totalCount.(TIME_folder_names{i}).voxFinal.(tissues{j}) = sum(data.Segmentation.final.(tissues{j}),"all");
            voxChangeAccumulator = [voxChangeAccumulator; voxChanges.(tissues{j})];
        end

        bar(voxChangeAccumulator);
        xlabel('Tissue Type');
        ylabel('Number of Voxel Changes');
        set(gca,'xticklabel',tissues);
        vox_agg_dir = 'VoxCountAgg';


        title(strrep(folder_name, "_", " "));

    else
        disp(['Segmentation.mat not found in folder: ', folder_name]);
    end

end
saveas(gcf,append(vox_agg_dir, '/TIME_Counts'))
save(append(vox_agg_dir, '/VoxCounts.mat'), 'totalCount')
close all


