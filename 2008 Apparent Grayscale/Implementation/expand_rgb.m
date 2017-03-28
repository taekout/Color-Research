function imageout = expand_rgb( imagein )
% imageout = expand_rgb( imagein )
%
% same as 'expand()' except "imagein" & "imageout" contain RGB colormaps

if isempty( imagein )
   imageout = zeros(1,1,3);
   return
end

r = expand( imagein(:,:,1) ); % 'expand' red colormap
g = expand( imagein(:,:,2) ); % 'expand' green colormap
b = expand( imagein(:,:,3) ); % 'expand' blue colormap

imageout = zeros( size(r,1), size(r,2), 3 ); % allocate space for "imageout"

% recombine into a single image
imageout(:,:,1) = r;
imageout(:,:,2) = g;
imageout(:,:,3) = b;