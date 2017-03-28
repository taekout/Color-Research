function main()

oPicture = imread('colororder.jpg');
oPicture = im2double(oPicture);
tones = decolorize(oPicture, 0.3, 25, 0.001); % The author uses 0.5 but he recommended 0.3 for normal color2gray usage.
imwrite(tones, 'result.jpg', 'jpg');
end
