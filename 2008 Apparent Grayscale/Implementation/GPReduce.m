%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reduce an image applying Gaussian Pyramid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IResult = GPReduce(I)

Wt = [0.0500    0.2500    0.4000    0.2500    0.0500];

dim = size(I);
if length(dim) == 3
    newdim = ceil(dim*0.5);
    newdim(3) = 3;
else
    newdim = ceil(dim*0.5);
end
IResult = zeros(newdim,class(I)); % Initialize the array in the beginning ..
m = -2 : 2 ;% n = m ;

switch length(dim)
	case 1
		%% Pad the boundaries.
		I = [ I(1) ; I(1) ;  I ; I(dim(1));  I(dim(1)) ];  % Add two rows towards the beginning and the end.
		for i = 0 : newdim(1) -1
			A = I(2*i+m+3).*Wt;
			IResult(i + 1)= sum(A(:));
		end
	case 2
		%% Pad the boundaries.
		I = [ I(1,:) ; I(1,:) ;  I ; I(dim(1),:);  I(dim(1),:) ];  % Add two rows towards the beginning and the end.
		I = [ I(:,1)  I(:,1)     I   I(:,dim(2))  I(:,dim(2)) ];  % Add two columns towards the beginning and the end.
		
		%Wt2 = Wt'*Wt;
        Wt2 = [0.0039 0.0156 0.0234 0.0156 0.0039 ; 0.0156 0.0625 0.0938 0.0625 0.0156 ; 0.0234 0.0938 0.1406 0.0938 0.0234 ; 0.0156 0.0625 0.0938 0.0625 0.0156 ; 0.0039 0.0156 0.0234 0.0156 0.0039 ];
		
		for i = 0 : newdim(1) -1
			for j = 0 : newdim(2) -1
				A = I(2*i+m+3,2*j+m+3).*Wt2;
				IResult(i + 1, j + 1) = sum(A(:));
			end
		end

	case 3
        IResult = Reduce_rgb(I);
        
%		Wt3 = ones(5,5,5);
%		for i = 1:5
%			Wt3(i,:,:) = Wt3(i,:,:) * Wt(i);
%			Wt3(:,i,:) = Wt3(:,i,:) * Wt(i);
%			Wt3(:,:,i) = Wt3(:,:,i) * Wt(i);
%		end
		
		%% Pad the boundaries.
%		I2 = zeros(dim+4,class(I));
%		I2(3:2+dim(1),3:2+dim(2),3:2+dim(3)) = I;
%		I2(1,:,:)=I2(3,:,:);I2(1,:,:)=I2(3,:,:);I2(end,:,:)=I2(end-2,:,:);I2(end-1,:,:)=I2(end-2,:,:);
%		I2(:,1,:)=I2(:,3,:);I2(:,2,:)=I2(:,3,:);I2(:,end,:)=I2(:,end-2,:);I2(:,end-1,:)=I2(:,end-2,:);
%		I2(:,:,1)=I2(:,:,3);I2(:,:,2)=I2(:,:,3);I2(:,:,end)=I2(:,:,end-2);I2(:,:,end-1)=I2(:,:,end-2);
%		I=I2; clear I2;

%		for i = 0 : newdim(1) -1
%			for j = 0 : newdim(2) -1
%				for k = 0 : newdim(3) - 1
%					A = I(2*i+m+3,2*j+m+3,2*k+m+3).*Wt3;
%					IResult(i+1,j+1,k+1) = sum(A(:));
%				end
%			end
%		end
end



function imageout = Reduce_rgb( imagein )
% imageout = reduce_rgb( imagein )
%
% same as 'reduce()' except "imagein" & "imageout" contain RGB colormaps

r = GPReduce( imagein(:,:,1) ); % 'reduce' red colormap
g = GPReduce( imagein(:,:,2) ); % 'reduce' green colormap
b = GPReduce( imagein(:,:,3) ); % 'reduce' blue colormap

imageout = zeros( size(r,1), size(r,2), 3 ); % allocate space for "imageout"

% recombine into a single image
imageout(:,:,1) = r;
imageout(:,:,2) = g;
imageout(:,:,3) = b;
