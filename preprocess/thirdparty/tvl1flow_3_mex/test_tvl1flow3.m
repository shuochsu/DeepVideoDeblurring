%% Script to test the optical flow estimation
clear all
close all

img1f = '/datos/research/deblurring/data/jue/videoDeblur/snow/input/0001.png';
img2f = '/datos/research/deblurring/data/jue/videoDeblur/snow/input/0002.png';

%img1 = single(rgb2gray(imread(img1f)));
%img2 = single(rgb2gray(imread(img2f)));

img1 = double((imread(img1f)));
img2 = double((imread(img2f)));


img1 = .299*img1(:,:,1) + .587*img1(:,:,2) + .114*img1(:,:,3);
img2 = .299*img2(:,:,1) + .587*img2(:,:,2) + .114*img2(:,:,3);

img1 = single(floor(img1));
img2 = single(floor(img2));


figure(1)
imagesc(img2)


%%
tau=0.25;
lambda = 0.1;
theta=0.3;
nscales=4;
zfactor=0.5;
nwarps=5;
epsilon=0.01;
verbose=1;

[m,n,~] = size(img1);
tagstruct.ImageLength = m;
tagstruct.ImageWidth = n;
tagstruct.Compression = Tiff.Compression.None;
tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample =  32;                        % float data
tagstruct.SamplesPerPixel = 2; %u and v (2 channels)
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;


flowF = tvl1flow(img1,img2,tau,lambda,theta,nscales,zfactor,nwarps, epsilon,verbose);

%Save into a TIFF file
t = Tiff('aux1.tiff', 'w');
t.setTag(tagstruct);
t.write(single(flowF));
t.close();

%% OLD

setenv('DYLD_LIBRARY_PATH', '/usr/local/lib'),
tvl1flowS='/datos/research/deblurring/code/heat/videodeblur/tvl1flow_3_mex/tvl1flow';

flowB = 'aux.tiff';
system(sprintf('%s %s %s %s 0 0.25 0.1 0.3 4 0.5 5 0.01 1',tvl1flowS,img1f,img2f,flowB)); % flowB.tiff


%%
aux1 = imread('aux1.tiff');
aux = imread('aux.tiff');

figure(1)
imagesc(flowF(:,:,1)-aux(:,:,1))

figure(2)
imagesc(flowF(:,:,2)-aux(:,:,2))

