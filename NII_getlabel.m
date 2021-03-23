function Results = NII_getlabel(nii_map, threshold)
% reference map: spm12 Neuromorphometrics, resolution: 1.5mm
% add reference map: yeo 7 networks (from conn dir), resolution: 1mm
if ( ~exist('threshold', 'var') || isempty(threshold) )
    threshold = 1; % default threshold set to 1
end
tic;
disp(['threshold is ', num2str(threshold)]);
% initate report to speedup
% [Results(1,1).spmNeuromorphometrics, Results(1,1).MaxRegion, Results(1,1).Yeo7Networks, ...
%     Results(1,1).MaxNet, Results(1,1).Threshold, Results(1,1).Resolution] = deal([]);
% Results(1,1).Resolution = [];

%% reference map: spm12 Neuromorphometrics, resolution: 1.5mm
% update : add other atlas support; reslice nii map to match atlas file
% spmNeuromorphometric resliced to 1*1*1 mm^3 (Nearest Neighborhood; Original: 1.5*1.5*1.5)
atlas_set = {'spmNeuromorphometrics'; 'AAL3v1';'cat12anatomy3'};

% reslice ica map to smp12 label nifti resolution （1*1*1）
% 不能reslice atlas到ica文件空间，会产生原index外的值（其他插值），或减少index（最近邻插值）
disp('Resize nii file to match label file (1mm * 1mm * 1mm) ing...');
resize_img(nii_map,[1 1 1], nan(2,3));

% load resliced ica_map (nii_map) data
[fpath, fname] = fileparts(nii_map);
rnii_map = fullfile(fpath, ['r', fname, '.nii']);

rnii_data = niftiread(rnii_map);

% delete resliced nii map
delete(rnii_map);

for i = 1:numel(atlas_set)
    % prepare

    % unzip .nii.gz
    disp(['Calculationg report for ',atlas_set{i}]);
    label_map = [atlas_set{i}, '.nii.gz'];
%     label_gz = [atlas_set{i}, '.nii.gz'];
%     gunzip(label_gz);
%     label_map = [atlas_set{i}, '.nii'];
    label_names = importdata([atlas_set{i}, '.txt']);
    if i == 1 % if spm neuromorphometrics, load specified label indexs
        load([atlas_set{i}, '_indexs.mat']);
    else
        label_indexs = 1:numel(label_names);
    end
    % read label nifti 
    % label_niistruct = spm_vol_nifti(label_map);
    % label_niimat = label_niistruct.private.dat(); % load 3d mat data
    label_data = niftiread(label_map); % niftiread can load nii.gz directly
    label_info = niftiinfo(label_map);
    label_dims = label_info.PixelDimensions;

    % start to calculate
    if numel(size(rnii_data)) == 3 % if only one volume
        command = ['[Results(1,1).', atlas_set{i}, ', Results(1,1).MaxRegion', atlas_set{i}, '] = getlabel_fromVolume(rnii_data, threshold, label_data, label_indexs, label_names);'];
        eval(command);
    elseif numel(size(rnii_data)) == 4 % if multiple volumes
%         Results = repmat(Results, size(rnii_data, 4), 1); % make row = ics'number
        for n = 1:size(rnii_data, 4)
            % loop for every volume
            command = ['[Results(n,1).', atlas_set{i}, ', Results(n,1).MaxRegion', atlas_set{i}, '] = getlabel_fromVolume(squeeze(rnii_data(:,:,:,n)), threshold, label_data, label_indexs, label_names);'];
            eval(command);
        end
    end
%     delete(label_map); % delete .nii
%     Results(1,1).Resolution = [Results(1,1).Resolution; ['Report based on ', atlas_set{i}, ' Label file Resolution ', num2str(label_dims)]];
end
Results(1,1).atlasResolution = [1, 1, 1];
% --------------------------atlas part end-------------------------------

%% reference map: yeo 7/17 networks
% voxel statistics reported based on original resolution
% conn_roisPath = fullfile(fileparts(which('conn')), 'utils', 'otherrois');
% label_map = fullfile(conn_roisPath, 'Yeo2011.nii');
% label_txt = fullfile(conn_roisPath, 'Yeo2011.txt');
nets_set = {'Yeo_7Nets'; 'Yeo_17Nets'};
for j = 1:numel(nets_set)
    gunzip([nets_set{j}, '.nii.gz']);
    label_map = [nets_set{j}, '.nii'];
    label_txt = [nets_set{j}, '.txt'];
    label_names = importdata(label_txt);
    label_indexs = 1:numel(label_names);

    % reslice yeomap to match nii map
    matlabbatch{1}.spm.spatial.coreg.write.ref = {[which(nii_map), ',1']};
    matlabbatch{1}.spm.spatial.coreg.write.source = {[label_map, ',1']};
    matlabbatch{1}.spm.spatial.coreg.write.roptions.interp = 0; % Nearest neighborhood
    matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.write.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix = 'rcoreg';

    spm('Defaults','fMRI');
    spm_jobman('initcfg');
    spm_jobman('run',matlabbatch);

    % load nifti data
    nii_data = niftiread(nii_map);
    % load resliced yeo map
    rlabel_map = ['rcoreg', label_map];

    label_data = niftiread(rlabel_map);

    % delete resliced label map
    delete(rlabel_map);

    if numel(size(nii_data)) == 3 % if only one volume
        %[Results(1,1).Yeo7Networks, Results(1,1).MaxNet] = getlabel_fromVolume(nii_data, threshold, label_data, label_indexs, label_names);
        command = ['[Results(1,1).', nets_set{j}, ', Results(1,1).MaxRegion', nets_set{j}, '] = getlabel_fromVolume(nii_data, threshold, label_data, label_indexs, label_names);'];
        eval(command);
    elseif numel(size(nii_data)) == 4 % if multiple volumes
        for n = 1:size(nii_data, 4)
            % loop for every volume
            %[Results(n,1).Yeo7Networks, Results(n,1).MaxNet] = getlabel_fromVolume(squeeze(nii_data(:,:,:,n)), threshold, label_data, label_indexs, label_names);
            command = ['[Results(n,1).', nets_set{j}, ', Results(n,1).MaxRegion', nets_set{j}, '] = getlabel_fromVolume(squeeze(nii_data(:,:,:,n)), threshold, label_data, label_indexs, label_names);'];
            eval(command);
        end
    end
    % get nii dimensions info
    nii_info = niftiinfo(nii_map);
    nii_dims = nii_info.PixelDimensions;
    % Results(1,1).Resolution = [Results(1,1).Resolution; ...
    %     {['Yeo 7 networks Report based on: Original Niftis Resolution ', num2str(nii_dims)]}];

    delete(label_map);
end
% End of Yeo 7 networks part
Results(1,1).Threshold = threshold;
Results(1,1).netsResolution = nii_dims(1:3);
% -----------------------End of code-------------------------------%
toc;
end



%% -----------------------------------------------------------------------
% nested function, get report from one volume
function [Report, maxRegion] = getlabel_fromVolume(rnii_data, threshold, label_data, label_indexs, label_names)
    % threshold to binary
    upper_thresh = label_data(rnii_data > threshold);
    upthresh_VoxelSizes = numel(find(rnii_data > threshold)); % Whole Voxel size upper Threshold
    upthresh_VoxelMass = sum(rnii_data(rnii_data > threshold)); % Whole Voxel mass upper Threshold
    % only left with label & > threshold
    upper_thresh = upper_thresh(upper_thresh ~= 0);
    % statistics
    uni_inds = unique(upper_thresh); % indexs of label upper threshold
    
	% to store numbers of voxels; percent of voxel size / whole voxel size (upper threshold)
	% ;mass of voxels (size * statistics value);
    [uni_voxel_nums, uni_voxel_mass, uni_voxel_perc] = deal(zeros(numel(uni_inds), 1)); 

    uni_labels = cell(numel(uni_inds), 1); % to store name of lable
    for i = 1:numel(uni_inds)
        % cal voxel sizes in 1.5 resolution
        uni_voxel_nums(i) = sum(upper_thresh == uni_inds(i) );
        % cal voxel masses in 1.5 resolution
        uni_voxel_mass(i) = sum(rnii_data(rnii_data > threshold & label_data  == uni_inds(i)));
        % find corresponding label
        uni_labels(i) = label_names(label_indexs == uni_inds(i));
        uni_voxel_perc(i) = uni_voxel_nums(i) / upthresh_VoxelSizes * 100;
    end
    %resort by voxel numbers, descend
    [VoxelSize, Ind] = sort(uni_voxel_nums, 'descend');
    anatLabel = uni_labels(Ind);
    VoxelMass = uni_voxel_mass(Ind);
    Voxel_upThresh_Perc = uni_voxel_perc(Ind);
    
    % add to last row, report without label's voxels info
    VoxelSize = [VoxelSize; upthresh_VoxelSizes - sum(VoxelSize)];
    Voxel_upThresh_Perc = [Voxel_upThresh_Perc; 100 - sum(Voxel_upThresh_Perc)];
    anatLabel = [anatLabel; {'Unlabeled upper threshold area'}];
    VoxelMass = [VoxelMass; upthresh_VoxelMass - sum(VoxelMass)];

    Report = table(VoxelSize, Voxel_upThresh_Perc, VoxelMass, anatLabel);
    % sort by descend
    [~, newInd] = sort(VoxelSize, 'descend');
    Report = Report(newInd, :);
    
    if VoxelSize(1) == 0 % if no upper threshold voxel with label
        maxRegion =NaN;
    else
        maxRegion = anatLabel(1);
    end
end
