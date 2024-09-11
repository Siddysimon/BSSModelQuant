
BIdir = '/Volumes/spiral/BSSLab/TMSTestRetest/Models/';
subjectsTIME = append('TIME',{'020','041','045','052','057','072','081','082','085'});
tissues = {'wm','gm','csf','bone','skin','eyes','air'};


% List contents of the directory
contents = dir(BIdir);

% Extract folder names using logical indexing
folder_names = {contents([contents.isdir]).name};

% Remove '.' and '..' from the list
subjectsBI  = folder_names(~ismember(folder_names, {'.', '..'}));

% replace_dir = 1;


% Loop through each folder name
for i = 1:numel(subjectsBI)
    subject_dir = join([BIdir subjectsBI{i}]);

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

        highest_number = max(unique_numbers); %Get highest run (final)

        % Find corresponding folder name
        highest_folder = sprintf('run%d', highest_number);

        runF = join([subject_dir '/' highest_folder '/m2m_' subjectsBI{i} '/mask_prep/']);

        % Check if the Segmentation.mat file exists in the folder
        if exist(runF, 'dir')
            dilatedDir = append(runF, 'dilated/');
            erodedDir = append(runF, 'eroded/');
            if ~exist(dilatedDir, 'dir')
                mkdir(dilatedDir)
            end
            if ~exist(erodedDir, 'dir')
                mkdir(erodedDir)
            end
            for j = 1:numel(tissues)
                if exist(append(dilatedDir, 'dilated_', tissueNifti.name), 'file') && exist(append(erodedDir, 'eroded_', tissueNifti.name), 'file')
                    continue
                end
                tissueNifti = dir(fullfile(runF,append('*', upper((tissues{j})), '*.nii*')));
                if length(tissueNifti) > 1
                    tissueNifti = tissueNifti(1);
                end
                mask = nifti_load(append(runF, tissueNifti.name));
                maskDilate = mask;
                maskErode = mask;
                maskDilate.vol = thickenbinvol(mask.vol, 2);
                maskErode.vol = thinbinvol(mask.vol, 2);
                nifti_save(maskDilate, append(dilatedDir, 'dilated_', tissueNifti.name));
                nifti_save(maskErode, append(erodedDir, 'eroded_', tissueNifti.name));
                
            end

        else
            disp(['Final run folder not found: ', runF]);
        end
    end
end



% Loop through each folder name
for i = 1:numel(subjectsTIME)
    runF = ['/Volumes/spiral/BSSLab/TIME/Modeling/' subjectsTIME{i} '/Model/headreco_T1/m2m_' subjectsTIME{i} '/mask_prep/'];

    % Check if the final run folder exists
    if exist(runF, 'dir')
            dilatedDir = append(runF, 'dilated/');
            erodedDir = append(runF, 'eroded/');
            if ~exist(dilatedDir, 'dir')
                mkdir(dilatedDir)
            end
            if ~exist(erodedDir, 'dir')
                mkdir(erodedDir)
            end
            for j = 1:numel(tissues)
                if exist(append(dilatedDir, 'dilated_', tissueNifti.name), 'file') && exist(append(erodedDir, 'eroded_', tissueNifti.name), 'file')
                    continue
                end

                tissueNifti = dir(fullfile(runF,append('*', upper((tissues{j})), '*.nii*')));
                if length(tissueNifti) > 1
                    tissueNifti = tissueNifti(1)
                end
                mask = nifti_load(append(runF, tissueNifti.name));
                maskDilate = mask;
                maskErode = mask;
                maskDilate.vol = thickenbinvol(mask.vol, 2);
                maskErode.vol = thinbinvol(mask.vol, 2);
                nifti_save(maskDilate, append(dilatedDir, 'dilated_', tissueNifti.name));
                nifti_save(maskErode, append(erodedDir, 'eroded_', tissueNifti.name));
            end

    else
        disp(['Final run folder not found: ', runF]);
    end

end


