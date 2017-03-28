originalImage=imread('pascu.jpg');
c=makecform('srgb2lab');
LAB=applycform(originalImage,c);
L=LAB(:,:,1);
imwrite(L, 'pascu_grey.jpg', 'jpg');