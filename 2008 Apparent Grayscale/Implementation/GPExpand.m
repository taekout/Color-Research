%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expand an image as per the Gaussian Pyramid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IResult = GPExpand(I, originalSize)
% originalSize is supposed to have 2 dimensions at all times

Wt = [0.0500    0.2500    0.4000    0.2500    0.0500];

dim = size(I);
%newdim = dim*2;
% care for the original Size
newdim = originalSize;
if length(dim) == 3
    newdim(3) = 3;
end
IResult = zeros(newdim,class(I)); % Initialize the array in the beginning ..
m = [-2:2];n=m;

switch length(dim)
	case 1
		%% Pad the boundaries.
		I = [ I(1) ;  I ;  I(dim(1))];
		for i = 0 : newdim(1) - 1
			pixeli = (i - m)/2 + 2;  idxi = find(floor(pixeli)==pixeli);
			A = I(pixeli(idxi)) .* Wt(m(idxi)+3);
			IResult(i + 1)= 2 * sum(A(:));
		end
	case 2
		%% Pad the boundaries.
		I = [ I(1,:) ;  I ;  I(dim(1),:) ];  % Pad the top and bottom rows.
		I = [ I(:,1)    I    I(:,dim(2)) ];  % Pad the left and right columns.
		% Wt2 = Wt'*Wt;
        Wt2 = [0.0039 0.0156 0.0234 0.0156 0.0039 ; 0.0156 0.0625 0.0938 0.0625 0.0156 ; 0.0234 0.0938 0.1406 0.0938 0.0234 ; 0.0156 0.0625 0.0938 0.0625 0.0156 ; 0.0039 0.0156 0.0234 0.0156 0.0039 ];

		for i = 0 : newdim(1) - 1
			for j = 0 : newdim(2) - 1
				pixeli = (i - m)/2 + 2;  idxi = find(floor(pixeli)==pixeli);
				pixelj = (j - m)/2 + 2;  idxj = find(floor(pixelj)==pixelj);
				A = I(pixeli(idxi),pixelj(idxj)) .* Wt2(m(idxi)+3,m(idxj)+3);
				IResult(i + 1, j + 1)= 4 * sum(A(:));
			end
		end
	case 3
        IResult = Expand_rgb(I, originalSize);
%		Wt3 = ones(5,5,5);
%		for i = 1:5
%			Wt3(i,:,:) = Wt3(i,:,:) * Wt(i);
%			Wt3(:,i,:) = Wt3(:,i,:) * Wt(i);
%			Wt3(:,:,i) = Wt3(:,:,i) * Wt(i);
%		end
		
%		%% Pad the boundaries
%		I2 = zeros(dim+2,class(I));
%		I2(2:1+dim(1),2:1+dim(2),2:1+dim(3)) = I;
%		I2(1,:,:)=I2(2,:,:);I2(end,:,:)=I2(end-1,:,:);
%		I2(:,1,:)=I2(:,2,:);I2(:,end,:)=I2(:,end-1,:);
%		I2(:,:,1)=I2(:,:,2);I2(:,:,end)=I2(:,:,end-1);
%		I=I2; clear I2;
		
%		for i = 0 : newdim(1) - 1
%			for j = 0 : newdim(2) - 1
%				for k = 0 : newdim(3) - 1
%					pixeli = (i - m)/2 + 2;  idxi = find(floor(pixeli)==pixeli);
%					pixelj = (j - m)/2 + 2;  idxj = find(floor(pixelj)==pixelj);
%					A = I(pixeli(idxi),pixelj(idxj),pixelk(idxk)) .* Wt3(m(idxi)+3,m(idxj)+3,m(idxk)+3);
%					IResult(i + 1, j + 1,k+1)= 8 * sum(A(:));
%				end
%			end
%		end
end

function imageout = Expand_rgb( imagein, originalSize )
% imageout = expand_rgb( imagein )
%
% same as 'expand()' except "imagein" & "imageout" contain RGB colormaps

if isempty( imagein )
   imageout = zeros(1,1,3);
   return
end

r = GPExpand( imagein(:,:,1), originalSize ); % 'expand' red colormap
g = GPExpand( imagein(:,:,2), originalSize ); % 'expand' green colormap
b = GPExpand( imagein(:,:,3), originalSize ); % 'expand' blue colormap

imageout = zeros( size(r,1), size(r,2), 3 ); % allocate space for "imageout"

% recombine into a single image
imageout(:,:,1) = r;
imageout(:,:,2) = g;
imageout(:,:,3) = b;