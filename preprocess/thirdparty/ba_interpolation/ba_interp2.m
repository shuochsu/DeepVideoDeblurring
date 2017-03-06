function Z = ba_interp2(F, X, Y, method)
% Fast nearest, bi-linear and bi-cubic interpolation for image data
%
% Usage:
% ------
%     Z = ba_interp2(F, X, Y, [method])
%
% where method is one off nearest, linear, or cubic.
%
% F    is a WxHxD Image with an arbitray number of layers D.
% X, Y are I_1 x ... x I_n matrices with the x and y coordinates to
%      interpolate.
% Z    is a I_1 x ... x I_n x D matrix, which contains the interpolated image channels.
%
% Notes:
% ------
% This method handles the border by repeating the closest values to the point accessed. 
% This is different from matlabs border handling.
%
% Example
% ------
%
%    %% The mandrills eye
%    clear
%    IMG=load('mandrill');
%    IMG = ind2rgb(IMG.X, IMG.map);
%    [Dx Dy] = meshgrid(130:0.1:250, -150:0.1:-50);
%    
%    R = [cos(pi/4) sin(pi/4); -sin(pi/4) cos(pi/4)];
%    RD = R * [Dx(:)'; Dy(:)'] + 250;
%    RDx = reshape(RD(1,:), size(Dx));
%    RDy = reshape(RD(2,:), size(Dy));
%    
%    methods = {'nearest', 'linear', 'cubic'};
%    la=nan(1,3);
%    for i=1:3
%      la(i) = subplot(2,2,i);
%      tic;
%      IMG_R = ba_interp2(IMG, RDx, RDy, methods{i});
%      elapsed=toc;
%      imshow(IMG_R);
%      title(sprintf('Rotation and zoom using %s interpolation took %gs', methods{i}, elapsed));
%    end
%    linkaxes(la);
%
% Licence:
% --------
% GPL
% (c) 2008 Brian Amberg
% http://www.brian-amberg.de/

  error('ERROR: The mex file was not compiled. Use  $ mex -O ba_interp2.cpp      to compile it');
end
