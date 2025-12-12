function [] = segment_nuclei_cellpose(data, exec_envir, varargin)
% segment_nuclei_cellpose(data, exec_envir, crop_corner, crop_size)
%
% Function to segment nuclei using cellpose. Requires  Deep Learning Toolbox™, 
%       Computer Vision Toolbox™, and the Medical Imaging Toolbox™ Interface 
%       for Cellpose Library. Find cellpose library download instructions here: 
%       https://www.mathworks.com/matlabcentral/fileexchange/130629-medical-imaging-toolbox-interface-for-cellpose-library?s_tid=srchtitle_support_results_3_cellpose 
%
% INPUT:    data = data structure minimally containing fields Source, 
%                   ImageFileListNuc, and PixRes.
%           exec_envir = execution environment of choice, either 'CPU' or
%                   'GPU' (if available, GPU runs 2-3x faster)
% OPTIONAL INPUTS:
%           crop_corner (varargin{1}) = upper right position to crop input image from
%           crop_size (varargin{2}) = size of XY crop window
%           cellseg (varargin{3}) = logical input to turn on optional cell segmentation 
%                   restriction parameter. 1 to limit nuclear segmenataion to regions 
%                   with cell segmentation, 0 to segment all objects in frame.
% 
% OUTPUT:   no output variables-- results saved to SegmentationData 
%
% L. Russell (updated 24 June, 2025)
od = cd; % original directory

%% CONSTANTS + LOAD DATA
t_in = 1; % initial time point (change as needed)
mn = 1; % movie index in data structure
% pix_res = data(mn).PixRes; % # of pixels / micron
% z_step = data(mn).Zstep; % # of microns per z-slice (usually 1)

% optional input to restrict results to regions with cell segmentation
if nargin == 5; cellseg = varargin{3}; else; cellseg = 0; end

% get list of raw image file names from data
imageFileList = data(mn).ImageFileListNuc;

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

    tic; % start clock, for performance testing

    % display current progress of processing
    fprintf('extracting @ timepoint %04d\n',t);
    % navigate to SegmentationData folder
    cd(data(mn).Source); cd('SegmentationData');
    cd(cframefoldername);

    %% (1) load data
    fprintf('1/6 (load data)'); % command line message

    imageNuc_full = [];
    % loop through each Z-plane and filter/ pre-process images
    % imageNuc_full = zeros(im_size); % pre-allocate storage for raw image loading
    for z = 1:size(imageFileList,2) % loop through Z
        imageNuc_temp = imread(imageFileList{t,z}); % read raw image
        imageNuc_full(:,:,z) = double(imageNuc_temp(:,:,1)); % convert from RGB to double
    end

    % if exactly 3 inputs, assume 512 x 512 cropping window:
    if nargin == 3
        crop_corner = varargin{1};
        imageNuc_full = imageNuc_full(crop_corner(1):(crop_corner(1)+512), ...
            crop_corner(2):(crop_corner(2)+512),:);
    elseif nargin > 3 % otherwise, use crop_size input:
        crop_corner = varargin{1}; crop_size = varargin{2};
        imageNuc_full = imageNuc_full(crop_corner(1):(crop_corner(1)+crop_size(1)), ...
            crop_corner(2):(crop_corner(2)+crop_size(2)),:);
    end

    % rescaling factor, based on z-spacing
    % z_scale = round(size(imageFileList,2) .* pix_res .* z_step); 
    z_scale = round(size(imageFileList,2));
    % resample volume  so voxels have equal dimensions:
    imageNuc_full = imresize3(imageNuc_full, ...
        [size(imageNuc_full,1), size(imageNuc_full,2), z_scale]);

    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % clear command line message
    
    %% (2) Pre-Process
    fprintf('2/6 (pre-process)'); % command line message

    % Apply a median filter to reduce noise
    imageNuc_adjust = imageNuc_full; % pre-allocate space for adjusted image
    imageNuc_adjust = medfilt3(imageNuc_adjust); % median filter
    imageNuc_adjust = imgaussfilt3(imageNuc_adjust, [2 2 1]); % gaussian filter
    imageNuc_adjust = mat2gray(imageNuc_adjust); % normalize to 1
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % clear command line message
    
    %% (3) cellpose
    fprintf('3/4 (cellpose)'); % command line message

    % call segmentCells3D function, from cellpose plugin library
    nuclei = segmentCells3D(cp_nuc, imageNuc_adjust); % default batch size

    % % % smaller batch size, if computational resources are limited:
    % nuclei = segmentCells3D(cp_nuc, imageNuc_adjust, 'GPUBatchSize', 4);

    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b') % clear command line message

    %% (4) Match tracking #'s to previous timepoint
    fprintf('4/4 (tracking)'); % command line message

     % navigate to current SegmentationData folder
    cd(data.Source); cd('SegmentationData'); cd(cframefoldername);
    % optionally mask out regions without matching cell segmentation:
    if cellseg == 1
        % load mask and erode slightly
        mask = load('mask.mat').mask; mask = imerode(mask, strel('disk',3));
        % find nuclei that overlap with mask and delete their masks before tracking
        ind = nuclei(mask); ind = unique(ind); ind(ind==0) = [];
        for i = 1:length(ind); nuclei(nuclei==ind(i)) = 0; end
    end

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

    toc; % print clock result, for performance testing
end % time loop

cd(od); % return to original directory
end
