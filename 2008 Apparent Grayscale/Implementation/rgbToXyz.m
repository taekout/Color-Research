function imageXYZ = rgbToXyz(imageRGB)

image(:,:, 1) = imageRGB(:,:, 1);
image(:,:, 2) = imageRGB(:,:, 2);
image(:,:, 3) = imageRGB(:,:, 3);

image = double(image) / 255;

var_R = image(:,:,1);
var_G = image(:,:,2);
var_B = image(:,:,3);
%var_R = ( image / 255 )        //R from 0 to 255
%var_G = ( G / 255 )        //G from 0 to 255
%var_B = ( B / 255 )        //B from 0 to 255

if ( var_R > 0.04045 ) var_R = ( ( var_R + 0.055 ) / 1.055 ) ^ 2.4
else                   var_R = var_R / 12.92
if ( var_G > 0.04045 ) var_G = ( ( var_G + 0.055 ) / 1.055 ) ^ 2.4
else                   var_G = var_G / 12.92
if ( var_B > 0.04045 ) var_B = ( ( var_B + 0.055 ) / 1.055 ) ^ 2.4
else                   var_B = var_B / 12.92
end

var_R = var_R * 100
var_G = var_G * 100
var_B = var_B * 100

%Observer. = 2��, Illuminant = D65
X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805
Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722
Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505