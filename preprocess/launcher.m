function launcher(alignmentType)
% clc; clear; close all
addpath(genpath('thirdparty'));

%%
if ~exist('alignmentType','var')
    alignmentType = 1; %0 for nowarp, 1 for optical flow, 2 for homography, 3 for similarity
end
path2dataset = '../dataset/qualitative_datasets';   

%%
datasets = dir(path2dataset);
dirFlags = [datasets.isdir];
datasets(~dirFlags) = [];
datasets = {datasets(3:end).name};

% permute and divide into training and testing sets
num_datasets = numel(datasets);
seed = 12;
rng(seed);
ds_idx = randperm(num_datasets);

alignment = '';
if alignmentType == 0
    alignment = '_nowarp';
elseif alignmentType == 1
    alignment = '_OF';
elseif alignmentType == 2
    alignment = '_homography';
end

path2output = @(stage,d) sprintf('../data/testing_real_all_nostab%s/%s/image_%s',alignment,stage,d);   

fn_in = @(dsi,type,fri) [path2dataset '/' datasets{ds_idx(dsi)} '/' type '/' sprintf('%05d.jpg',fri)];
fn_out = @(stage,fri,frid) sprintf('%s/%05d.jpg',path2output(stage,frid),fri);

ds_range = 1:num_datasets;

for l = -2:2
    for i = 1:num_datasets
        checkDir(path2output(datasets{i},num2str(l)));
    end
end

%%
for ii = 1:length(ds_range)
    fr_cnt = 0;
    i = ds_range(ii);
    % get the frame range
    datasets{ds_idx(i)}
    files = dir([path2dataset '/' datasets{ds_idx(i)} '/input/*.jpg']);
    if isempty(files)
        files = dir([path2dataset '/' datasets{ds_idx(i)} '/input/*.png']);
    end
    if ~isempty(files)
        [~,ststr,~] = fileparts(files(1).name);
        [~,enstr,~] = fileparts(files(end).name);
        start_frame = str2num(ststr);
        end_frame = str2num(enstr);
        frame_range = start_frame:min(start_frame+99,end_frame);
        num_frame = numel(frame_range);

        fr_idx = floor(linspace(frame_range(1),frame_range(end),num_frame));
        for j = 1:num_frame
            fr_cnt = fr_cnt+1;
            % save image_1 to image_5
            v0 = im2double(imread(fn_in(i,'input',fr_idx(j)+0)));
            v0g = single(rgb2gray(v0));
            [h,w,~] = size(v0);
            
            for l = -2:2
                if l ~= 0
                    vi = im2double(imread(fn_in(i,'input',max(min(fr_idx(j)+l,frame_range(end)),frame_range(1)))));
                    vig = single(rgb2gray(vi));
                    if alignmentType == 0
                        v_i0 = vi;
                    elseif alignmentType == 1
                        flo_i0 = genFlow(v0g, vig);
                        [v_i0, ~] = warpToRef(v0, vi, flo_i0);
                    elseif alignmentType == 2
                        v_i0 = homographyAlignment(v0,vi,0);
                    elseif alignmentType == 3
                        v_i0 = similarityAlignment(v0,vi,0);
                    end
                else
                    v_i0 = v0;
                end
                imwrite(v_i0, fn_out(datasets{ds_idx(i)},fr_cnt,num2str(l)));
            end
        end
    end        
end
