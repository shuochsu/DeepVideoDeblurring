%% Find Image Rotation and Scale Using Automated Feature Matching
% This example shows how to automatically align two images that differ by a
% rotation and a scale change. It closely parallels another example titled
% <matlab:showdemo('RotationFitgeotransExample') Find Image Rotation and Scale>. 
% Instead of using a manual approach to register the two images, it
% utilizes feature-based techniques found in the Computer Vision System
% Toolbox(TM) to automate the registration process.
%
% In this example, you will use |detectSURFFeatures| and 
% |vision.GeometricTransformEstimator| System object to recover rotation 
% angle and scale factor of a distorted image. You will then transform the 
% distorted image to recover the original image.

% Copyright 1993-2014 The MathWorks, Inc. 

function recovered = homographyAlignment(original, distorted, verbose)
% clc; clear; close all

%% Step 1: Read Image
% Bring an image into the workspace.
if isempty(original)
    original = imread('00001.jpg');
end
original_c = original;
original = rgb2gray(original);
if verbose
    imshow(original);
    text(size(original,2),size(original,1)+15, ...
        'Image courtesy of Massachusetts Institute of Technology', ...
        'FontSize',7,'HorizontalAlignment','right');
end

%% Step 2: Resize and Rotate the Image
if isempty(distorted)
    % scale = 0.7;
    % J = imresize(original, scale); % Try varying the scale factor.
    % 
    % theta = 30;
    % distorted = imrotate(J,theta); % Try varying the angle, theta.
    distorted = rgb2gray(imread('00003.jpg'));
end
distorted_c = distorted;
distorted = rgb2gray(distorted);
if verbose
    figure, imshow(distorted)
end

%%
% You can experiment by varying the scale and rotation of the input image.
% However, note that there is a limit to the amount you can vary the scale
% before the feature detector fails to find enough features.

%% Step 3: Find Matching Features Between Images
% Detect features in both images.
ptsOriginal  = detectSURFFeatures(original);
ptsDistorted = detectSURFFeatures(distorted);

%%
% Extract feature descriptors.
[featuresOriginal,   validPtsOriginal]  = extractFeatures(original,  ptsOriginal);
[featuresDistorted, validPtsDistorted]  = extractFeatures(distorted, ptsDistorted);

%%
% Match features by using their descriptors.
indexPairs = matchFeatures(featuresOriginal, featuresDistorted);

%%
% Retrieve locations of corresponding points for each image.
matchedOriginal  = validPtsOriginal(indexPairs(:,1));
matchedDistorted = validPtsDistorted(indexPairs(:,2));

%%
% Show point matches. Notice the presence of outliers.
if verbose
    figure;
    showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
    title('Putatively matched points (including outliers)');
end

%% Step 4: Estimate Transformation
% Find a transformation corresponding to the matching point pairs using the
% statistically robust M-estimator SAmple Consensus (MSAC) algorithm, which
% is a variant of the RANSAC algorithm. It removes outliers while computing
% the transformation matrix. You may see varying results of the
% transformation computation because of the random sampling employed by the
% MSAC algorithm.
[tform, inlierDistorted, inlierOriginal, status] = estimateGeometricTransform(...
    matchedDistorted, matchedOriginal, 'projective');
if status
    recovered = original_c;
    return;
end

%%
% Display matching point pairs used in the computation of the
% transformation matrix.
if verbose
    figure;
    showMatchedFeatures(original,distorted, inlierOriginal, inlierDistorted);
    title('Matching points (inliers only)');
    legend('ptsOriginal','ptsDistorted');
end

%% Step 5: Solve for Scale and Angle
% Use the geometric transform, TFORM, to recover 
% the scale and angle. Since we computed the transformation from the
% distorted to the original image, we need to compute its inverse to 
% recover the distortion.
%
%  Let sc = scale*cos(theta)
%  Let ss = scale*sin(theta)
%
%  Then, Tinv = [sc -ss  0;
%                ss  sc  0;
%                tx  ty  1]
%
%  where tx and ty are x and y translations, respectively.
%

%%
% Compute the inverse transformation matrix.
Tinv  = tform.invert.T;

ss = Tinv(2,1);
sc = Tinv(1,1);
scale_recovered = sqrt(ss*ss + sc*sc);
theta_recovered = atan2(ss,sc)*180/pi;

%%
% The recovered values should match your scale and angle values selected in
% *Step 2: Resize and Rotate the Image*.

%% Step 6: Recover the Original Image
% Recover the original image by transforming the distorted image.
outputView = imref2d(size(original_c));
recovered = zeros(size(original_c));

[h,w,l] = size(original_c);
[mx,my] = meshgrid(1:w,1:h);
tmp = [mx(:)';my(:)';ones(1,h*w)];
tmp = (Tinv)'*tmp;
X = reshape(tmp(1,:)./tmp(3,:),h,w);
Y = reshape(tmp(2,:)./tmp(3,:),h,w);
recovered = ba_interp2(single(distorted_c),X,Y,'linear');

% for i=1:l
%     recovered(:,:,i) = interp2(mx,my,distorted_c(:,:,i),X,Y,'linear',0);
% end

% for ic = 1:size(original_c,3)
%     recovered(:,:,ic)  = imwarp(squeeze(distorted_c(:,:,ic)),tform,'OutputView',outputView,'FillValues',0,'SmoothEdges',true);
% end

%%
% Compare |recovered| to |original| by looking at them side-by-side in a montage.
if verbose
    figure, imshowpair(original_c,recovered,'montage')
%%
% The |recovered| (right) image quality does not match the |original| (left)
% image because of the distortion and recovery process. In particular, the 
% image shrinking causes loss of information. The artifacts around the edges are 
% due to the limited accuracy of the transformation. If you were to detect 
% more points in *Step 4: Find Matching Features Between Images*, 
% the transformation would be more accurate. For example, we could have
% used a corner detector, |detectFASTFeatures|, to complement the SURF 
% feature detector which finds blobs. Image content and image size also 
% impact the number of detected features.

    displayEndOfDemoMessage(mfilename)
end
