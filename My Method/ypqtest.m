function ypqtest()
FileName = input('Hello! Please enter the file name you like to convert! : ', 's');
originalImage = imread(FileName);  % Import color data as uint8 in the range [0,255]
%originalImage = double(originalImage)/255;        % Cast to double in the range [0,1]
originalImage = im2double(originalImage); % Also, cast to double in the range [0,1]
c=makecform('srgb2lab');
LAB=applycform(originalImage, c);
WhitePointX = 0.950456; WhitePointY = 1; WhitePointZ = 1.088754;
disp('Now processing ...');
originalImage = imread(FileName);  % Import color data as uint8 in the range [0,255]
%originalImage = double(originalImage)/255;        % Cast to double in the range [0,1]
originalImage = im2double(originalImage); % Also, cast to double in the range [0,1]

% apply gamma function to convert to sRGB
R = originalImage(:,:,1);
G = originalImage(:,:,2);
B = originalImage(:,:,3);
grey_from_rgb=sqrt(R.^2 + G.^2 + B.^2);
%originalImage_rgb = originalImage;
%R(R <= 0.03928) = R(R <=0.03928) ./ 12.92;
%R(R > 0.03928) = ((R(R > 0.03928) + 0.055) ./ 1.055) .^ 2.4;
%G(G <= 0.03928) = G(G <= 0.03928) ./ 12.92;
%G(G > 0.03928) = ((G(G > 0.03928) + 0.055) ./ 1.055 ) .^ 2.4;
%B(B <= 0.03928) = B(B <= 0.03928) ./ 12.92;
%B(B > 0.03928) = ((B(B > 0.03928) + 0.055) ./ 1.055) .^ 2.4;
%originalImage(:,:,1) = R;
%originalImage(:,:,2) = G;
%originalImage(:,:,3) = B;
% now SRGB.

c = makecform('srgb2xyz');
XYZ = applycform(originalImage, c);

LUV=xyz2luv(XYZ);

%%%%%%%%%%%%%%%%
% I got it!! LUV matrix has a lot of  NaNs. it's because of when X and
% (X + 15Y + 3Z) are 0 or Y and (X + 15Y + 3Z) are 0.
% As a result, the computation becomes 0 / 0 = NaN.
% Eventually, it makes u, v arrays have lots of NaNs.
% tweak the values very bit later.
% also ask ppl why it would happen. i mean 'ask if 0 value is allowed.
% and what i should do in the situation.
%%%%%%%%%%%%%%%%
%VAC = ConvertLuvToLvac(LUV);
%% NAyatani
L = LAB(:,:,1);
A = LAB(:,:,2);
B = LAB(:,:,3);
X = zeros(size(L));
Y = zeros(size(L));
Z = zeros(size(L));
HA = zeros(size(L));
HB = zeros(size(L));

L = L ./ 2.55;
A = (A - 128.0) ./ 1.27;
B = (B - 128.0) ./ 1.27;

P = (L + 16) / 116;
YYN = P .^ 3;

Y(YYN > 0.008856) = WhitePointY * YYN(YYN > 0.008856);
HA(YYN > 0.008856) = P(YYN > 0.008856) +  A(YYN > 0.008856) ./ 500;
X(YYN > 0.008856) = WhitePointX .* HA(YYN > 0.008856) .* HA(YYN > 0.008856) .* HA(YYN > 0.008856);
HB(YYN > 0.008856) = P(YYN > 0.008856) - B(YYN > 0.008856) / 200;
Z(YYN > 0.008856) = WhitePointZ .* HB(YYN > 0.008856) .* HB(YYN > 0.008856) .* HB(YYN > 0.008856);

Y(YYN > 0.008856) = WhitePointY .* L(YYN > 0.008856) ./ 903.3;
SQYYN = zeros(size(L));
SQYYN(YYN <= 0.008856) = (L(YYN <= 0.008856) ./ 903.3) .^ (1/3);
HA(YYN <= 0.008856) = A(YYN <= 0.008856) ./ 500 ./ 7.787 + SQYYN(YYN <= 0.008856) ;
X(YYN <= 0.008856) = WhitePointX .* HA(YYN <= 0.008856) .* HA(YYN <= 0.008856) .* HA(YYN <= 0.008856);
HB(YYN <= 0.008856) = SQYYN(YYN <= 0.008856) - B(YYN <= 0.008856) ./ 200 ./ 7.787;
Z(YYN <= 0.008856) = WhitePointZ .* HB(YYN <= 0.008856) .* HB(YYN <= 0.008856) .* HB(YYN <= 0.008856);
XYZ(:,:,1) = X;
XYZ(:,:,2) = Y;
XYZ(:,:,3) = Z;

Uc = 0.20917;
Vc = 0.48810;
U = zeros(size(X));
V = zeros(size(X));

U = (4 * X) ./ (X + 15 * Y + 3 * Z) - Uc;
V = (9 * Y) ./ (X + 15 * Y + 3 * Z) - Vc;
HUE = atan2(V, U);
QHUE = -0.01585 - 0.03017 * cos(HUE) - 0.04556 * cos(2*HUE) - 0.02667 * cos(3*HUE) - 0.00295 * cos(4*HUE) + 0.14592 * sin(HUE) + 0.05084 * sin(2 * HUE) - 0.01900 * sin(3*HUE) - 0.00764 * sin(4*HUE);

% Adapting Luminance Dependency
K_br = 0.2717 * ((6.469 + 6.362 * 20)^0.4495) / ((6.469 + 20)^0.4495);
S_uv = 13 * (   (    ( (U).^2) + ((V).^2)     ).^0.5  );
gamma = 1 + (-0.1340 .* QHUE + 0.0872 .* K_br) .* S_uv;
L = gamma .* L;

L = L .* 2.5599;
L(L>255) = 255; L(L<0) = 0;

% I am here!!!!!!!!!!!!!!!!!!!!!!!!!
% according to the source, after blurring, convert lab(G) -> RGB

%GreyImage = MappingLnToGreyscale(VAC); % Not necessary - this function is
%definitely wrong. I don't know why... but Not mapping it to Lab seems
%perfectly the same as the image of the paper. and in GIMP C source, it
%also just uses it as L of Lab

% for the computation, calculate Lab values in advance.

% Invoke the second version of a remade recursive f.

L = L ./ 100; % Nayatani = L
% Define the YPQ color space
colorconvert=[0.2989360212937753847527155, 0.5870430744511212909351327,  0.1140209042551033243121518;
              0.5,                         0.5,                         -1;
              1,                          -1,                            0]';
colorrevert= [1,                           0.1140209042551033243121518,  0.6440535265786729530912086;
              1,                           0.1140209042551033243121518, -0.3559464734213270469087914;
              1,                          -0.8859790957448966756878482,  0.1440535265786729530912086]';


frame=[size(originalImage,1), size(originalImage,2)];
pixels=frame(1) * frame(2);
originalImage=reshape(originalImage,[pixels,3]);
image=originalImage*colorconvert; % now YPQ
YPQ=image(:,1);
YPQ=reshape(L,frame);
figure(1);
imshow(L); %Nayatani
figure(2);
imshow(YPQ(:,:,1)); %YPQ
figure(3);
LAB=im2double(LAB); %LAB
LAB=LAB./100;
imshow(LAB(:,:,1));
figure(4);
imshow(grey_from_rgb); %RGB AVG


end



function Image = xyz2luv(Image)
% Convert to CIE L*u*v* (CIELUV)
WhitePoint = [0.950456,1,1.088754];
WhitePointU = 0.20917; %(4*WhitePoint(1))./(WhitePoint(1) + 15*WhitePoint(2) + 3*WhitePoint(3));
WhitePointV = 0.48810; %(9*WhitePoint(2))./(WhitePoint(1) + 15*WhitePoint(2) + 3*WhitePoint(3));

X = Image(:,:,1);
Y = Image(:,:,2);
Z = Image(:,:,3);

tmp = X + (15 * Y) + (3 * Z);

% compute U, V
% == >if (tmp == 0 && X == 0) || (tmp == 0 && Y == 0) 
%    error('wowowowowowowoowoowwowoowowowofidasojflksd;flsajdfl;sadfj;ajsf');
%end    
U = (4 * X) ./ tmp;
V = (9 * Y) ./ tmp;

% compute L
L = zeros(size(Y));
bState = Y > (6/29)^3;
L(bState) = 116 .* Y(bState) .^ (1/3) - 16;
L(~bState) = (29/3)^3 .* Y(~bState);

% Sometimes in the event that X and Y are 0, U and V are also computed as
% Nans. to prevent it, set Nan values to 0(default);
tmp2 = isnan(U);
tmp3 = isnan(V);
U(tmp2) = 0;
V(tmp3) = 0;

Image(:,:,1) = L;   % L*
Image(:,:,2) = U;   % u*
Image(:,:,3) = V;   % v*
end