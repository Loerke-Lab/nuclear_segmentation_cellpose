
function [] = segment_nuclei_cellpose_LS(data, exec_envir)
% Function to segment nuclei using cellpose. Requires  Deep Learning Toolbox™, 
%       Computer Vision Toolbox™, and the Medical Imaging Toolbox™ Interface 
%       for Cellpose Library. Find cellpose library download instructions here: 
%       https://www.mathworks.com/matlabcentral/fileexchange/130629-medical-imaging-toolbox-interface-for-cellpose-library?s_tid=srchtitle_support_results_3_cellpose 
%
% INPUT:    data = data structure minimally containing fields Source, 
%                   ImageFileListNuc, and PixRes. This implementations
%                   assumes the file convention associated with lattice
%                   light sheet data, where each time point is saved as a
%                   volumetric image.
%           exec_envir = execution environment of choice, either 'CPU' or
%                   'GPU' (if available, GPU runs 2-3x faster)
% 
% OUTPUT:   no output variables-- results saved to SegmentationData 
%
% L. Russell (updated 24 June, 2025)
od = cd; % original directory
warning('off','all'); % silence command line warnings

%% CONSTANTS + LOAD DATA
t_in = 1; % initial time point (change as needed)
mn = 1; % movie index in data structure
ds_factor = 2; % down-sampling factor
pix_res = data(mn).PixRes; % # of pixels / micron
imageFileList = data(mn).ImageFileListNuc; % list of raw image file names

% load first image, to extract dimensions
fprintf('Loading First Volume...'); % command line message
cd(data(mn).Source); % navigate to Source
im_temp = tiffreadVolume(imageFileList{1}); % load the first image volume
im_size = size(im_temp); % get size of volumes
% z_scale = round(im_size(3) * pix_res/2); % optional rescaling factor
z_scale = im_size(3);
z_lim = round(im_size(3)) - 13; % segment full volume
% z_lim = round(im_size(3) .* (2/3)); % only segment apical most 2/3 of volume
% z_lim = round(im_size(3) .* (2/3)); % only segment apical most 1/2 of volume
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % clear command line message

% load cellpose model, with preferred environment
fprintf('Loading Cellpose model...'); % command line message
if strcmp(exec_envir, 'GPU')
    cp_nuc = cellpose('ExecutionEnvironment','gpu', 'Model','nuclei');
elseif strcmp(exec_envir, 'CPU')
    cp_nuc = cellpose('ExecutionEnvironment','cpu', 'Model','nuclei');
end
fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % clear command line message

%% MAIN
cd(data(mn).Source); cd('SegmentationData') % navigate to segmentation folder

% loop through all time points:
t = t_in;
cframefoldername = sprintf('frame%04d',t); % naming convention of seg folders
count = t_in; % counting variable, for tracking
while exist(cframefoldername)==7
% for t = t_in  % alternate loop start, for debugging

    % display current progress of processing
    fprintf('extracting @ timepoint %04d\n',t);

    %% (1) load data
    fprintf('1/6 (load data)'); % command line message

    % loop through each Z-plane and filter/ pre-process images
    cd(data(mn).Source); % navigate to Source
    % load image from Source
    imageNuc_full = tiffreadVolume(imageFileList{t});
    % limit analysis to top 2/3 of volume
    imageNuc_full = imageNuc_full(:,:,1:z_lim);
    % % resample volume  so voxels have equal dimensions
    % imageNuc_full = imresize3(imageNuc_full, [im_size(1), im_size(2), z_scale]);
    % downsize volume by factor of 3 in x/y
    resize_x = round(size(imageNuc_full,1) / ds_factor);
    resize_y = round(size(imageNuc_full,2) / ds_factor);
    resize_mat = [resize_x, resize_y, size(imageNuc_full,3)];
    imageNuc_full = imresize3(imageNuc_full,resize_mat);

    % navigate to SegmentationData folder
    cd(data(mn).Source); cd('SegmentationData');
    cd(cframefoldername);

    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % clear command line message
    
    %% (2) Pre-Process
    fprintf('2/6 (pre-process)'); % command line message

    % % background subtraction:
    % sigma = 1.5; % interpolation sigma for pre-processing
    % sigma_background = round(10 .* pix_res); % sigma for background subtraction
    % % Apply a gaussian filter to smooth signal
    % imageNuc_adjust = imgaussfilt3(imageNuc_full, sigma);
    % % subtract a large sigma filtered version, for local noise sub:
    % imageNuc_adjust = imageNuc_adjust - imgaussfilt3(imageNuc_adjust, sigma_background);

    % Apply a median filter to reduce noise
    imageNuc_adjust = medfilt3(imageNuc_full); % median filter
    imageNuc_adjust = mat2gray(imageNuc_adjust); % normalize to 1
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % clear command line message
    
    %% (3) cellpose

    fprintf('3/4 (cellpose)'); % command line message
    
    % segment volume based on exec_envir input
    if strcmp(exec_envir, 'GPU')
        % % call segmentCells3D function, from cellpose plugin library
        % nuclei = segmentCells3D(cp_nuc, imageNuc_adjust); % default batch size
        % smaller batch size, if computational resources are limited:
        nuclei = segmentCells3D(cp_nuc, imageNuc_adjust, 'GPUBatchSize', 4);
    elseif strcmp(exec_envir, 'CPU')
        % call segmentCells3D function, from cellpose plugin library
        nuclei = segmentCells3D(cp_nuc, imageNuc_adjust);
    end

    % upsample segmentation to match original image
    nuclei = imresize3(nuclei, [im_size(1), im_size(2), z_lim], 'box');
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % clear command line message

       %% (4) Match tracking #'s to previous timepoint
    fprintf('4/4 (tracking)'); % command line message

    % mask out regions without matching cell segmentation
    % navigate to current SegmentationData folder
    cd(data.Source); cd('SegmentationData'); cd(cframefoldername);
    % load mask and erode slightly
    mask = load('mask.mat').mask; mask = imerode(mask, strel('disk',3));
    % find nuclei that overlap with mask and delete their masks before tracking
    ind = nuclei(mask); ind = unique(ind); ind(ind==0) = [];
    for i = 1:length(ind); nuclei(nuclei==ind(i)) = 0; end

    if count > t_in
        % if starting at frame > t=1, load previous frame for tracking
        cframefoldername_last = sprintf('frame%04d',t-1); % naming convention of seg folders
        cd(data.Source); cd('SegmentationData'); cd(cframefoldername_last);
        nuclei_last = load('nuclei_cp.mat').nuclei_cp;
        nuclei_last = max(nuclei_last,[],3); % max projection to get ID list

        nuclei_last = nuclei_last .* 10000; % re-scale tracking from previous frame
        nuc_max_proj = max(nuclei, [], 3); % max projection to get list of tracking IDs
        ind_list = unique(nuc_max_proj(:)); ind_list(ind_list==0) = [];
        
        % on first tracked time point, initialize counter variable as max + 1. 
        % on subsequent loops, keep counting up from there.
        if count == t_in+1
            counter = max(ind_list) + 1; % updating variable for tracking IDs
        end

        % loop through each nucleus in current frame:
        for n = 1:length(ind_list)
            fprintf(' nuc #%04d/%04d', n, length(ind_list)); % command line message
            % find pixels in last frame that overlap with nuc n in current frame
             track_list = nuclei_last(nuc_max_proj == ind_list(n)); 
             % find IDs with highest overlap between frames
             ID = mode(track_list(:));
             if ID == 0
                 ID = counter .* 10000;
                 counter = counter + 1;
             end
             % re-label nuclei ID #'s to match
             nuclei(nuclei==ind_list(n)) = ID;
             fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % clear command line message
        end
        nuclei = nuclei ./ 10000; % re-scale tracking IDs
    end

    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % clear command line message

    %% SAVING + PLOTTING
    fprintf('Saving Segmentation...'); % command line progress message

    nuclei_cp = nuclei;
    % % OPTIONAL: plot results after each time point is segmented
    % for z = 1:size(nuclei,3)
    %     imshowpair(mat2gray(imageNuc_full(:,:,z)), imbinarize(nuclei(:,:,z)));
    %     title_stg = sprintf('t = %04d', t);
    %     title(title_stg);
    %     pause(0.01);
    % end

    % write result to file
    cd(data.Source); cd('SegmentationData'); cd(cframefoldername);
    % % save under usual variable name, to match ELSA and NSNO methods
    % save('nuclei.mat', 'nuclei', '-v7.3'); 
    % save under alternate variable name, for performance testing/ validation
    save('nuclei_cp.mat', 'nuclei_cp', '-v7.3');

    % update t and folder name
    count = count + 1;
    t=t+1;
    cframefoldername = sprintf('frame%04d',t);
    % return to upper directory
    cd(data.Source); cd('SegmentationData');

    % delete command line message (saving)
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');

    % delete message 'extracting @ timepoint %04d' before next loop:
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); 
end % time loop

cd(od); % return to original directory
end

