function [out, nan_map] = backwardsWarp(im,flow,is_nan2zero)
[h,w,l] = size(im);
[mx,my] = meshgrid(1:w,1:h);

X = mx + flow(:,:,1);
Y = my + flow(:,:,2);

% old approach 
if 0
    out = zeros(h,w,l);
    nan_map = false(h,w,l);

    for i=1:l
        out(:,:,i) = interp2(mx,my,im(:,:,i),X,Y,'linear');

        nan_map(:,:,i) = isnan(out(:,:,i));
        if is_nan2zero
            out(:,:,i) = interp2(mx,my,im(:,:,i),X,Y,'linear',0);
        end
    end
end

% new faster method (no nans)
out = ba_interp2(im,X,Y,'linear');
nan_map = isnan(out);
