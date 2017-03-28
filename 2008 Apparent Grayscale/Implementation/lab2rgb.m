%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert a lab image to a rgb image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rgb = lab2rgb(lab)

kappa = 903.3;
epsilon = 0.008856;
[h, w, d] = size(lab);

xr = zeros(h, w);
yr = zeros(h, w);
zr = zeros(h, w);

fy = zeros(h, w);

% yr
cc = lab(:, :, 1);
inds = find(cc > kappa * epsilon);
yr(inds) = ( (cc(inds)+16)/116 ) .^ 3;
inds = find(cc < kappa * epsilon);
yr(inds) = cc(inds) ./ kappa;

% fy
inds = find(yr > epsilon);
fy(inds) = (cc(inds) + 16) / 116;
inds = find(yr <= epsilon);
fy(inds) = (kappa.*yr(inds) + 16) / 116;

% fx, fz
fx = ( lab(:, :, 2)./500 ) + fy;
fz = fy - ( lab(:, :, 3)./200 );

% xr, zr
inds = find(fx.^3 > epsilon);
xr(inds) = fx(inds).^3;
inds = find(fx.^3 <= epsilon);
xr(inds) = (116.*fx(inds) - 16) / kappa;

inds = find(fz.^3 > epsilon);
zr(inds) = fz(inds).^3;
inds = find(fz.^3 <= epsilon);
zr(inds) = (116.*fz(inds) - 16) / kappa;

xyz(:, :, 1) = xr * 95.047;
xyz(:, :, 2) = yr * 100.000;
xyz(:, :, 3) = zr * 108.883;

rgb = xyz2rgb(xyz);

end % function lab2rgb


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convert a xyz image to a rgb image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rgb = xyz2rgb(xyz)

xyz = xyz / 100.0;

M = [0.412424 0.212656 0.0193324; 0.357579 0.715158 0.119193; 0.180464 0.0721856 0.950444];
M = inv(M);

rgb(:, :, 1) = xyz(:, :, 1) * M(1, 1) + xyz(:, :, 2) * M(2, 1) + xyz(:, :, 3) * M(3, 1);
rgb(:, :, 2) = xyz(:, :, 1) * M(1, 2) + xyz(:, :, 2) * M(2, 2) + xyz(:, :, 3) * M(3, 2);
rgb(:, :, 3) = xyz(:, :, 1) * M(1, 3) + xyz(:, :, 2) * M(2, 3) + xyz(:, :, 3) * M(3, 3);

% for i = 1:3
%     cc = rgb(:, :, i);
%     
%     inds = find(cc > 0.0031308);
%     cc(inds) = 1.055 * ( cc(inds) .^ (1.0 / 2.4) ) - 0.055;
%     
%     inds = find(cc <= 0.0031308);
%     cc(inds) = 12.92 * cc(inds);
%     
%     rgb(:, :, i) = cc;
% end

end % function xyz2rgb