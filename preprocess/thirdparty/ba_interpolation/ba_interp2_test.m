%% compile
%mex -O ba_interp2.cpp

%% The mandrills eye
clear
IMG=load('mandrill');
IMG = ind2rgb(IMG.X, IMG.map);
[Dx Dy] = meshgrid(-25:0.1:25, -25:0.1:25);

R = [cos(pi/18) sin(pi/18); -sin(pi/18) cos(pi/18)];
RD = R * [Dx(:)'; Dy(:)'];
RDx = reshape(RD(1,:), size(Dx))+320;
RDy = reshape(RD(2,:), size(Dy))+42;

methods = {'nearest', 'linear', 'cubic'};
la=nan(1,3);
for i=1:3
  tic;
  IMG_R = ba_interp2(IMG, RDx, RDy, methods{i});
  elapsed=toc;

  IMG_R_M = nan(size(IMG_R));
  tic;
  for c=1:size(IMG,3)
    IMG_R_M(:,:,c) = interp2(IMG(:,:,c), RDx, RDy, ['*' methods{i}]);
  end
  elapsed_matlab=toc;
    
  la(i) = subplot(2,2,i);
  imshow(IMG_R);
  difference = max(max(max(abs(IMG_R-IMG_R_M))));
  title({...
    sprintf('Rotation and zoom using %s interpolation took %gs.', methods{i}, elapsed), ...
    sprintf('Speedup against inbuilt method by a factor of %g', elapsed_matlab/elapsed), ...
    sprintf('The maximum difference with the built in method was %g', difference)});
  end
linkaxes(la(:));

%% The mandrills eye single
clear
IMG=load('mandrill');
IMG = single(ind2rgb(IMG.X, IMG.map));
[Dx Dy] = meshgrid(single(-25:0.1:25), single(-25:0.1:25));

R = single([cos(pi/18) sin(pi/18); -sin(pi/18) cos(pi/18)]);
RD = R * [Dx(:)'; Dy(:)'];
RDx = reshape(RD(1,:), size(Dx))+320;
RDy = reshape(RD(2,:), size(Dy))+42;

methods = {'nearest', 'linear', 'cubic'};
la=nan(1,3);
for i=1:3
  tic;
  IMG_R = ba_interp2(IMG, RDx, RDy, methods{i});
  elapsed=toc;

  IMG_R_M = nan(size(IMG_R));
  tic;
  for c=1:size(IMG,3)
    IMG_R_M(:,:,c) = interp2(IMG(:,:,c), RDx, RDy, ['*' methods{i}]);
  end
  elapsed_matlab=toc;
    
  la(i) = subplot(2,2,i);
  imshow(IMG_R);
  difference = max(max(max(abs(IMG_R-IMG_R_M))));
  title({...
    sprintf('Rotation and zoom using %s interpolation took %gs.', methods{i}, elapsed), ...
    sprintf('Speedup against inbuilt method by a factor of %g', elapsed_matlab/elapsed), ...
    sprintf('The maximum difference with the built in method was %g', difference)});
  end
linkaxes(la(:));


%% Cumulative effect of resampling.
clear
IMG=load('mandrill');
IMG = ind2rgb(IMG.X, IMG.map);
[Dx Dy] = meshgrid(-349.5:1:349.5, -349.5:1:349.5);

R = [cos(pi/18) sin(pi/18); -sin(pi/18) cos(pi/18)];
RD = R * [Dx(:)'; Dy(:)'];
RDx = reshape(RD(1,:), size(Dx))+350;
RDy = reshape(RD(2,:), size(Dy))+350;

methods = {'nearest', 'linear', 'cubic'};
la=nan(1,numel(methods));
for i=1:numel(methods)
  IMG_R{i} = IMG;
  IMG_R{i}(end+1:700,:,:) = 0;
  IMG_R{i}(:,end+1:700,:) = 0;
  IMG_R{i} = circshift(IMG_R{i}, ceil([[700-size(IMG,1) 700-size(IMG,2)]/2, 0]));
  imshow(IMG_R{i});
end
[X Y] = meshgrid(1:700, 1:700);
mask = ((X-350).^2+(Y-350).^2)<350^2;
for j=1:36
  for i=1:numel(methods)
    tic;
    IMG_R{i} = ba_interp2(IMG_R{i}, RDx, RDy, methods{i});
    IMG_R{i}(repmat(~mask, [1 1 3])) = 0;
    elapsed=toc;

    la(i) = subplot(1,3,i);
    imshow(IMG_R{i});
    title(sprintf('Rotating a mandrill %g times by 10degree using %s interpolation.', j, methods{i}));
    axis([100 600 100 600]);
  end
  drawnow
end
linkaxes(la(:));

for i=1:numel(methods)
  imwrite(IMG_R{i}, sprintf('mandril_%s.png', methods{i}), 'Alpha', double(mask));
end
