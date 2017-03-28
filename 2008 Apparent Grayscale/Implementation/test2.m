function test2()

testImage = imread('sunrise.bmp');  % Import color data as uint8 in the range [0,255]
testImage = double(testImage)/255;        % Cast to double in the range [0,1]


XYZ = colorspace('rgb->xyz', testImage);
LUV = colorspace('rgb->luv', testImage);

VAC = ConvertLuvToLvac(LUV);

TheNewY = MappingLnToGreyscale(VAC);

TheNewYReduced =  GPReduce(TheNewY);
%TheNewYExpanded1 = GPExpand(TheNewYReduced1);

%for nI = 1 : 10
%    TheNewYReduced2 =  GPReduce(TheNewYExpanded2);
%    TheNewYExpanded2 = GPExpand(TheNewYReduced2);
%end

