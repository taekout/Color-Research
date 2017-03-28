function  [RGBlinear] = InverseGammaEncoding(RGBnonlinear)
% Return a value after the inverse process of gamma encoding.

RGBlinear = eye(3,1);

RGBlinear(1) = (RGBnonlinear(1)^2.150) * 0.9521;
RGBlinear(2) = (RGBnonlinear(2)^2.0750) * 0.9583;
RGBlinear(3) = (RGBnonlinear(3)^2.5) * 0.9554;




