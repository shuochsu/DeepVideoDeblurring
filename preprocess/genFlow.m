function flo_i0 = genFlow(v0g, vig)

[h,w,~] = size(v0g);

minimum_flow_size = 200;         % length of the smaller dimension in the optical flow
DWNSMP = minimum_flow_size/min(h,w);  
TAU = 0.2;                       % 0.25; (higher = faster, less accurate)
LAMBDA = 0.15;                   % 0.15; (smaller = smoother)
THETA = 0.3;
NSCALES = 6;                     % 5;
ZFACTOR = 0.5;
NWARPS = 7;                      % 5;
EPSILON = 0.01;
VERBOSE = 0;

flo_i0 = tvl1flow(v0g, vig, DWNSMP, TAU, LAMBDA, THETA,...
                    NSCALES, ZFACTOR, NWARPS, EPSILON, VERBOSE);
                
flo_i0_original_scale = zeros(h,w,2);
flo_i0_original_scale(:,:,1) = imresize(flo_i0(:,:,1), [h, w], 'bicubic');
flo_i0_original_scale(:,:,2) = imresize(flo_i0(:,:,2), [h, w], 'bicubic');
flo_i0 = flo_i0_original_scale;

flo_i0 = flo_i0/DWNSMP;
