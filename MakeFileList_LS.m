
function [data] = MakeFileList_LS(data)
% This function generates an image file list cell array from the TIF files
% in the directory data.Source or in the current directory (default).
%
% INPUT:    data: data structure with fields ImageFileList and
%               ImageFileListNuc pre-allocated to empty variables.
%
% OUTPUT:   data: results are written to input data structure, and empty
%               segmentation folders are saved
od = cd; % store original directory

%% CONSTANTS
mn = 1; % index of movie in data structure
path = data(mn).Source; % get Source

%% MAIN
% check if path was inputted and change directory to that path
if nargin == 1
    cd(path)
end

% get list of all tif files in Source folder
Mylist = dir('*.tif');
MylistNames = {Mylist.name}';
Nfiles = numel(MylistNames);
idxListCh0 = 1;
FrameVecCh0 = [];
idxListCh1 = 1;
FrameVecCh1 = [];

for ii=1:Nfiles
    NamePartsCells = strsplit(MylistNames{ii},'_');
    % naming conventions for MOSAIC light sheet images, 2025 visit
    Channel = NamePartsCells{6};
    Frame   = NamePartsCells{12};
    
    Frame(1:5) = [];
    FrameInt = str2num(Frame);
    
    % if strcmpi(Channel,'ch0')
    % if strcmpi(Channel,'GFP')   
    if strcmpi(Channel,'488nm') 
        ListCh0{idxListCh0} = MylistNames{ii};
        idxListCh0 = idxListCh0 + 1;
        FrameVecCh0 = [FrameVecCh0;FrameInt];
    end
    % if strcmpi(Channel,'ch1')
    % if strcmpi(Channel,'mCherry')   
    if strcmpi(Channel,'560nm') 
        ListCh1{idxListCh1} = MylistNames{ii};
        idxListCh1 = idxListCh1 + 1;
        FrameVecCh1 = [FrameVecCh1;FrameInt];
    end
end

if any(diff(FrameVecCh0)~=1)
    error('Check frame order or missing frames in Channel 0')
end
if any(diff(FrameVecCh0)~=1)
    error('Check frame order or missing frames in Channel 1')
end
ListCh0 = ListCh0';
ListCh1 = ListCh1';

% write result to data structure
data(mn).ImageFileList = ListCh1;
data(mn).ImageFileListNuc = ListCh0;

% pre-allocate segmentation data folders
cd(data(mn).Source);
% if a segmentation results folder doesn't already exist, create one
if  ~ (exist('SegmentationData')==7) 
    [~,~,~] = mkdir(path,'SegmentationData');
end
cd('SegmentationData');
path2 = cd;
% loop over desired time points and create subfolders
for t=1:length(ListCh0)
    cframefoldername = sprintf('frame%04d',t);
    if  ~ (exist(cframefoldername)==7) 
        [~,~,~] = mkdir(path2,cframefoldername);
    end
end

cd(od); % return to original directory
end

%% HELPER FUNCTIONS
%
% function [] = makeSegmentationDataFiles(Nframes)
%     for t=1:Nframes
%         mkdir(sprintf('frame%04d',t));
%     end
% end