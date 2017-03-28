function imageout = reduce_rgb( imagein )
% imageout = reduce_rgb( imagein )
%
% same as 'reduce()' except "imagein" & "imageout" contain RGB colormaps

r = reduce( imagein(:,:,1) ); % 'reduce' red colormap
g = reduce( imagein(:,:,2) ); % 'reduce' green colormap
b = reduce( imagein(:,:,3) ); % 'reduce' blue colormap

imageout = zeros( size(r,1), size(r,2), 3 ); % allocate space for "imageout"

% recombine into a single image
imageout(:,:,1) = r;
imageout(:,:,2) = g;
imageout(:,:,3) = b;