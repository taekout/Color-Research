function main2()
% Assume that the number of the dimensions of an image is an even number
% there is the need to fix the conversion from XYZ to LAB
global nMaximumDepth ;
WhitePointX = 0.950456; WhitePointY = 1; WhitePointZ = 1.088754;


FileName = input('Hello! Please enter the file name you like to convert! : ', 's');
Depth = input('What depth(How many levels of L pyramid) would you like ?: ');
K = zeros(Depth);
for  i = 1 : Depth
    disp(sprintf('%s %d %s', 'Please enter K array value - index #', i, ' : '));
    K(i) = input('Value : ');
end
P_GAMMA = input('Please enter P value : ');
disp('Now processing ...');

%% IIR_INIT()
%B , B1, B2,B3,B0,R ,Q
% P



%static struct iir_param
%
% gdouble B,b1,b2,b3,b0,r,q;
% gdouble *p;
%} iir; 
% double r value is from the arguement which may be from c2g_params.radius
% which may be from  param[3].data.d_float (param is a GimpParam, data is a
% union(int, char, string , etc..
%static void iir_init(double r)
%{
%}
  
%  iir.r = r;
%  gdouble q;
%  if ( r >= 2.5) q = 0.98711 * r - 0.96330;
%  else q = 3.97156 - 4.14554 * sqrt(1.0-0.26891 * r);
  
%  iir.q = q;
%  iir.b0 = 1.57825 + ((0.422205 * q  + 1.4281) * q + 2.44413) *  q;
%  iir.b1 = ((1.26661 * q +2.85619) * q + 2.44413) * q / iir.b0;
%  iir.b2 = - ((1.26661*q +1.4281) * q * q ) / iir.b0;
%  iir.b3 = 0.422205 * q * q * q / iir.b0;
%  iir.B = 1.0 - (iir.b1 + iir.b2 + iir.b3);



%% Extract_lab

nMaximumDepth = Depth;
originalImage = imread(FileName);  % Import color data as uint8 in the range [0,255]
originalImage = double(originalImage)/255;        % Cast to double in the range [0,1]
RED = originalImage(:,:,1);
GREEN = originalImage(:,:,2);
BLUE = originalImage(:,:,3);

%{
%   apply gamma function to convert to sRGB
RED(RED<=0.03928) = RED(RED<=0.03928) ./ 12.92;
RED(RED> 0.03928) = ( (RED(RED>0.03928) + 0.055) ./ 1.055 ) .^ (2.4);
GREEN(GREEN<=0.03928) = GREEN(GREEN<=0.03928) ./ 12.92;
GREEN(GREEN> 0.03928) = ( (GREEN(GREEN>0.03928) + 0.055) ./ 1.055 ) .^ (2.4);
BLUE(BLUE<=0.03928) = BLUE(BLUE<=0.03928) ./ 12.92;
BLUE(BLUE> 0.03928) = ( (BLUE(BLUE>0.03928) + 0.055) ./ 1.055 ) .^ 2.4;
%}

%      //apply gamma function to convert to sRGB
%      if ( red <=0.03928 ) red = red/12.92;
%      else red = pow((red+0.055)/1.055,2.4);
%      if ( green <=0.03928 ) green = green/12.92;
%      else green = pow((green+0.055)/1.055,2.4);
%      if ( blue <=0.03928 ) blue = blue/12.92;
%      else blue = pow((blue+0.055)/1.055,2.4);
L = zeros(size(RED));
FTX = zeros(size(RED));
FTY = zeros(size(RED));
FTZ = zeros(size(RED));
X = 0.431 * RED + 0.342 * GREEN + 0.178 * BLUE;
Y = 0.222 * RED + 0.707 * GREEN + 0.071 * BLUE;
Z = 0.020 * RED + 0.130 * GREEN + 0.939 * BLUE;
%      x = 0.431 * red + 0.342 * green + 0.178 * blue;
%      y = 0.222 * red + 0.707 * green + 0.071 * blue;
%      z = 0.020 * red + 0.130 * green + 0.939 * blue;


TY = Y ./ WhitePointY;
L(TY > 0.008856) = 116 * TY(TY > 0.008856) .^ (1/3) - 16;
FTY(Y > 0.008856) = TY(Y > 0.008856) .^ 1/3;

L(Y <= 0.008856) = 903.3 .* TY(Y <= 0.008856);
FTY(Y <= 0.008856) = 7.78 .* TY(Y <= 0.008856) + 16.0 / 116.0;
%      if (( ty = y/Yn ) > 0.008856)
%        {
%          l   = 116.0 * cbrt( ty ) - 16.0;
%          fty = cbrt( ty );
%        }
%      else
%        {
%          l   = 903.3 * ty;
%          fty = 7.78*ty + sixteenth;
%        }

TX = X ./ WhitePointX;
FTX(TX > 0.008856) = TX(TX > 0.008856) .^ 1/3;
FTX(TX <=0.008856) = 7.78 * TX(TX <=0.008856) + 16.0 / 116.0;
TZ = X ./ WhitePointX;
FTZ(TZ > 0.008856) = TZ(TZ > 0.008856) .^ 1/3;
FTZ(TZ <=0.008856) = 7.78 * TZ(TZ <=0.008856) + 16.0 / 116.0;
%      ftx = ((tx = x/Xn) > 0.008856) ? cbrt (tx) : 7.78 * tx + sixteenth;
%      ftz = ((tz = z/Zn) > 0.008856) ? cbrt (tz) : 7.78 * tz + sixteenth;

%      lab_dst[0] = (guchar) (l * 2.5599);
%      lab_dst[1] = (guchar) (128.0 + (ftx - fty) * 635 );
%      lab_dst[2] = (guchar) (128.0 + (fty - ftz) * 254 );

%      rgb_src += offset;
%      lab_dst += offset;
%    }

LAB_DST = zeros(size(originalImage));
LAB_DST(:,:,1) = L * 2.5599;
LAB_DST(:,:,2) = 128 + (FTX - FTY) * 635;
LAB_DST(:,:,3) = 128 + (FTY - FTZ) * 254;

LAB = LAB_DST;
L = LAB(:,:,1);
A = LAB(:,:,2);
B = LAB(:,:,3);


%% NAyatani
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

LAB_G = zeros(size(LAB));
L = L .* 2.5599;
L(L>255) = 255; L(L<0) = 0;

LAB_G(:,:,1) = L;
LAB_G(:,:,2) = 128;
LAB_G(:,:,3) = 128;
%%

% I am here!!!!!!!!!!!!!!!!!!!!!!!!!
% according to the source, after blurring, convert lab(G) -> RGB

%GreyImage = MappingLnToGreyscale(VAC); % Not necessary - this function is
%definitely wrong. I don't know why... but Not mapping it to Lab seems
%perfectly the same as the image of the paper. and in GIMP C source, it
%also just uses it as L of Lab
GreyImage = LAB_G(:,:,1);
% convert it
%GreyImage(:,:,1) = GreyImage(:,:,1) ./ 255;
%GreyImage(:,:,2) = (GreyImage(:,:,2) - 128) ./255;
%GreyImage(:,:,3) = (GreyImage(:,:,3) - 128) ./255;

GreyImage = (GreyImage - 128) ./ 127;

%c = makecform('lab2srgb');
%DimensionalGreyImage = applycform(GreyImage, c);
% for the computation, calculate Lab values in advance.

% Invoke the second version of a remade recursive f.
Grey_Final = GreyImage + Recursive_f(originalImage, GreyImage, Depth, K, P_GAMMA);

figure;
imshow(originalImage);
axis image


figure;
imshow(GreyImage);%DimensionalGreyImage);
axis image
figure;
imshow(Grey_Final);
axis image

tmp_filename = strcat(FileName, '_Grey_', num2str(P_GAMMA), '_');
for i = 1 : Depth
    tmp_filename = strcat(tmp_filename, num2str(K(i)));
end

tmp_filename = strcat(tmp_filename, '_MAIN2.bmp');
imwrite(GreyImage, tmp_filename, 'bmp');
tmp_filename = strcat(FileName, '_Int_', num2str(P_GAMMA), '_');
for i = 1 : Depth
    tmp_filename = strcat(tmp_filename, num2str(K(i)));
end
tmp_filename = strcat(tmp_filename, '_MAIN2.bmp');
imwrite(Grey_Final, tmp_filename, 'bmp');






% The rest of the job is here.
% Now we are supposed to get a new value matrix.. actually the final value G'
% Try to test the values...

function Grey_adjusted = Recursive_f(Image, GreyImage, N, K, P)
% every recursive invocation, Image and GreyImage become half the size of
% them.
N = N - 1;
if N == -1
    Grey_adjusted = zeros(size(GreyImage), class(GreyImage));
    return
end

% Compute Hi(Laplacian bandpass), Kn(a constant per each depth), delta
tmp = size(GreyImage);
tmp = [tmp(1), tmp(2)];
reduced_G = GPReduce(GreyImage); expanded_G = GPExpand(reduced_G, tmp);
Hn_G = GreyImage - expanded_G;
% compute delta
tmp = size(Image); tmp = [tmp(1), tmp(2)];
reduced_Image = GPReduce(Image); expanded_Image = GPExpand(reduced_Image, tmp); clear tmp;
Hn_Image = Image - expanded_Image;

%%%% We gotta check the bandpass of G and I... their absolute values are
%%%% too different, which caused this code not work.

% check abs later. abs(Hn_Image);
c = makecform('srgb2lab');
lab_Laplacian = applycform(abs(Hn_Image), c);
%lab_Laplacian = colorspace('rgb->lab', abs(Hn_Image));

% Just in case, some modification is applied to Lab values.
% It seems like the ranges of Lab here are 0 < [L,a,b] < 1 even though they
% are practically 0 ~ 100, -xx ~ +xx , -xx ~ +xx.
% Let us set it to 0 ~ Lab ~ 1.
%lab_Laplacian = lab_Laplacian ;   %./ 100;
numerator = (lab_Laplacian(:,:,1) .^ 2 + lab_Laplacian(:,:,2) .^ 2 + lab_Laplacian(:,:,3) .^ 2) .^ 0.5;
denominator = abs(Hn_G);
% (denominator X= 0) : the denominator can't be zero because is is a
% denominator. We concluded that we will put regularization values such as
% almost zero and a similar value to the numerator.
denominator(denominator == 0) = 0.00001;
% there has got to be something wrong with the numerator
delta = (numerator ./ denominator) .^ P;
% ! Program ! = sometimes delta(i) value is extremely high. it's because
% sometimes a denominator and a numerator is almost the same. but not
% exactly the same... think about how to take care of this.

global nMaximumDepth;

tmp =  (Hn_G .* K(nMaximumDepth - N) .* delta) ;

Grey_adjusted = tmp + GPExpand(Recursive_f(reduced_Image, reduced_G, N, K, P), size(GreyImage));



function Image = luv2xyz(Image)
% Convert to CIE XYZ from 'SrcSpace'
WhitePoint = [0.950456,1,1.088754];  


% Convert CIE L*uv to XYZ
WhitePointU = 0.20917; %(4*WhitePoint(1))./(WhitePoint(1) + 15*WhitePoint(2) + 3*WhitePoint(3));
WhitePointV = 0.48810; %(9*WhitePoint(2))./(WhitePoint(1) + 15*WhitePoint(2) + 3*WhitePoint(3));

L = Image(:,:,1);
Y(L>8) = ((L(L>8) + 16)/116) .^ 3;
Y(L<=8) = L(L<=8) * ((3/29) .^ 3);
X = Y * 9 * WhitePointU / (4 * WhitePointV);
Y = Y * ( ( 12 - 3 * WhitePointU - 20 * WhitePointV) / (4 * WhitePointV) );


Image(:,:,1) = -(9*Y.*U)./((U-4).*V - U.*V);                  % X
Image(:,:,2) = Y;                                             % Y
Image(:,:,3) = (9*Y - (15*V.*Y) - (V.*Image(:,:,1)))./(3*V);  % Z


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
return;  

function AdjustedGrey = Recursive_adjustment(image, GreyImage, N, K, P)
% I stopped developing this, so this is still incomplete.
% Local adjustment function
% image = the original image
% GreyImage = the grey image in VAC format calculated from the original I
% N = how many laplacian pyramid depths.
% K = user defined values
% P = how many times contrast is squared. (should be 0 <= P <= 1)

N = N - 1;
if N == -1
    AdjustedGrey = zeros(size(GreyImage), class(GreyImage));
    return
end

% Compute Hi(G_L*)
reducedG = GPReduce(GreyImage);
expandedG = GPExpand(reducedG, size(GreyImage));
Hi_G = GreyImage - expandedG;

% Compute Hi(Image)
reducedImage = GPReduce(image);
tmp = size(image); tmp = [tmp(1), tmp(2)];
expandedImage = GPExpand(reducedImage, tmp); clear tmp;
deltaE = image - expandedImage; % the bandpass Image == deltaE (depth - N)

% a (-) value causes some value out of the range of Lab
deltaE = abs(deltaE);
deltaE = colorspace('rgb->lab', deltaE);
image_L = deltaE(:,:, 1); % extract L
image_a = deltaE(:,:, 2); % extract a
image_b = deltaE(:,:, 3); % extract b

% It is time to compute delta_N
numerator = (image_L .^ 2 + image_a .^ 2 + image_b .^ 2) .^ (0.5) ;
denominator = abs(Hi_G);
% (denominator X= 0) : the denominator can't be zero because is is a
% denominator. We concluded that we will put regularization values such as
% almost zero and a similar value to the numerator.
denominator(denominator == 0) = numerator(denominator == 0); % ./ 0.5;

disparate_N = numerator ./ denominator;
disparate_N = disparate_N .^ P;


% Compute the return value
maximum_Ndepth = 4;
AdjustedGrey = K(maximum_Ndepth - N) * disparate_N .* Hi_G;
AdjustedGrey = AdjustedGrey + Recursive_adjustment(image, AdjustedGrey, N, K, P);





