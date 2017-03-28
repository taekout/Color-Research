
function [tones,tones2] = ColorReordering3(preprocessedY1,picture,effect,scale,noise,perceptuality)


% Examine inputs
frame=[size(picture,1), size(picture,2)];
pixels=frame(1)*frame(2);
if nargin<2 || isempty(effect)
    effect=0.5;
end;
if nargin<3 || isempty(scale)
    scale=sqrt(2*min(frame));
end;
if nargin<4 || isempty(noise)
    noise=0.001;
end;

% Reset the random number generator
randn('state',0);
tolerance=100*eps;

% Define the YPQ color space
colorconvert=[0.2989360212937753847527155, 0.5870430744511212909351327,  0.1140209042551033243121518;
              0.5,                         0.5,                         -1;
              1,                          -1,                            0]';
colorrevert= [1,                           0.1140209042551033243121518,  0.6440535265786729530912086;
              1,                           0.1140209042551033243121518, -0.3559464734213270469087914;
              1,                          -0.8859790957448966756878482,  0.1440535265786729530912086]';
colorspan=[ 0,  1;
           -1,  1;
           -1,  1];
maxluminance=1;
scaleluminance=0.66856793424088827189; % 1 / sqrt(1^2 + 1^2 + 1^2) ==> Maximum value in RGB -> The distance 1.7xxx ==> 1/ 1.7xxxx
maxsaturation=1.1180339887498948482; 
alter=effect*(maxluminance/maxsaturation);


% Covert picture to the YPQ color space 
picture=reshape(picture,[pixels,3]); % 1024 X 768 X 3 ==> (786432) X 3
image=picture*colorconvert; % now YPQ
original=image;
chroma=sqrt(image(:,2).*image(:,2)+image(:,3).*image(:,3));

% Pair each pixel with a randomly chosen sample site
mesh=reshape(cat(3,repmat((1:frame(1))',[1,frame(2)]),repmat((1:frame(2)),[frame(1),1])),[pixels,2]);
% Basically, I think that the reason why it gets mesh is because it needs
% to pick a random displacement value.
% Below is to get random displacement index.
% the look values are not just going to be similar to the x, y index
% sometimes because y values are multiplied by the width.
displace=(scale*sqrt(2/pi))*randn(pixels,2);
look=round(mesh+displace);
redo=find((look(:,1)<1));
look(redo,1)=2-rem(look(redo,1),frame(1)-1);
redo=find((look(:,2)<2));
look(redo,2)=2-rem(look(redo,2),frame(2)-1);
redo=find((look(:,1)>frame(1)));
look(redo,1)=frame(1)-1-rem(look(redo,1)-2,frame(1)-1);
redo=find((look(:,2)>frame(2)));
look(redo,2)=frame(2)-1-rem(look(redo,2)-2,frame(2)-1);
look=look(:,1)+frame(1)*(look(:,2)-1);

%Prepare CieLAB
c=makecform('srgb2lab');
LAB=applycform(picture,c);
% Calculate the color differences between the paired pixels
preprocessedY1=reshape(preprocessedY1,[pixels,1]);
delta=image-image(look,:); % delta_YPQ
delta_perceptual=LAB-LAB(look,:);
delta_perceptual=reshape(delta_perceptual, [frame, 3]);
delta_perceptual=delta_perceptual(:,:,1).^3+delta_perceptual(:,:,2).^3+delta_perceptual(:,:,3);
delta_perceptual=reshape(delta_perceptual, size(delta(:,1)));
contrastchange=abs(delta(:,1));
contrastdirection=sign(delta(:,1).*(1-perceptuality)+delta_perceptual.*perceptuality); % Oi = sign(delta(YPQ))
% I can change the direction by defining sign(delta(:,1))  % delta == LAB - LAB(look, :);
% delta = LAB(:,:,1) .^ 3 + LAB(:,:,2) .^ 3 + LAB(:,:,3) .^ 3;
colordifference_rgb=picture-picture(look,:); % delta D in paper       %%% Taekyu Shin #1
%%%%%%%%%%%%%%%%% HERE You used to use just LAB! Please walk this through again!!!!!!!!
colordifference_lab=LAB-LAB(look,:); % delta D (LAB)                  %%% Taekyu Shin #2
colordifference_rgb=sqrt(sum(colordifference_rgb.*colordifference_rgb,2))+eps; % Taekyu Shin #3
colordifference_lab=sqrt(sum(colordifference_lab.*colordifference_lab,2))+eps; % Taekyu Shin #4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% So far so good.

% Derive a chromatic axis from the weighted sum of chromatic differences between paired pixels
weight_fromrgb=1-((contrastchange/scaleluminance)./colordifference_rgb); % Ci in paper
weight_fromlab=1-((contrastchange/scaleluminance)./colordifference_lab); % Ci from LAB  TAEKYU Shin #5
weight_fromrgb(find(colordifference_rgb<tolerance))=0;
weight_fromlab(find(colordifference_lab<tolerance))=0; % Taekyu Shin #6
axis_fromrgb=weight_fromrgb.*contrastdirection;
axis_fromlab=weight_fromlab.*contrastdirection;  % Taekyu Shin #7
axis_fromrgb=delta(:,2:3).*[axis_fromrgb, axis_fromrgb]; %P,Q component
axis_fromlab=delta(:,2:3).*[axis_fromlab, axis_fromlab]; %P,Q component - Taekyu Shin #8
axis_fromrgb=sum(axis_fromrgb,1);           % Pi*delta_p + Qi*delta_q
axis_fromlab=sum(axis_fromlab,1);           % Taekyu Shin #9

% Project the chromatic content of the picture onto the chromatic axis
projection_fromrgb=image(:,2)*axis_fromrgb(1)+image(:,3)*axis_fromrgb(2); % The final Ci numerator
projection_fromlab=image(:,2)*axis_fromlab(1)+image(:,3)*axis_fromlab(2); % The final Ci numerator - Taekyu Shin #10
projection_fromrgb=projection_fromrgb/(quantiles(abs(projection_fromrgb),1-noise)+tolerance); % The final Ci divides by denominator.
projection_fromlab=projection_fromlab/(quantiles(abs(projection_fromlab),1-noise)+tolerance); % Taekyu Shin #11

% Combine the achromatic tones with the projected chromatic colors and adjust the dynamic range
% Quantile function cuts out the outlying values beyong proper dynamic range.
preprocessedY_fromrgb=preprocessedY1+effect*projection_fromrgb;
preprocessedY_fromlab=preprocessedY1+effect*projection_fromlab; % Taekyu Shin #12
%%%%%%%%%%% From here, you might want to redo. First of all, quantile(image(:,1),..) might be a good idea because Nayatani
%%%%%%%%%%% range exceeds 1 already. You might want to ignore the outliers.
imagerange_fromrgb=quantiles(preprocessedY_fromrgb,[noise; 1-noise]);
imagerange_fromlab=quantiles(preprocessedY_fromrgb,[noise; 1-noise]); % Taekyu Shin #13
preprocessedY_fromrgb=(preprocessedY_fromrgb-imagerange_fromrgb(1))/(imagerange_fromrgb(2)-imagerange_fromrgb(1)+tolerance);
preprocessedY_fromlab=(preprocessedY_fromlab-imagerange_fromlab(1))/(imagerange_fromlab(2)-imagerange_fromlab(1)+tolerance); % Taekyu Shin #14
targetrange=effect*[0; maxluminance]+(1-effect)*quantiles(original(:,1),[noise; 1-noise]);
preprocessedY_fromrgb=targetrange(1)+(preprocessedY_fromrgb*(targetrange(2)-targetrange(1)+tolerance));
preprocessedY_fromlab=targetrange(1)+(preprocessedY_fromlab*(targetrange(2)-targetrange(1)+tolerance)); % Taekyu Shin #15
preprocessedY_fromrgb=min(max(preprocessedY_fromrgb,original(:,1)-alter.*chroma),original(:,1)+alter.*chroma);
preprocessedY_fromlab=min(max(preprocessedY_fromlab,original(:,1)-alter.*chroma),original(:,1)+alter.*chroma); % Taekyu Shin #16
preprocessedY_fromrgb=min(max(preprocessedY_fromrgb,0),maxluminance);
preprocessedY_fromlab=min(max(preprocessedY_fromlab,0),maxluminance); % Taekyu Shin #17


% Return the results
tones_fromrgb=preprocessedY_fromrgb/maxluminance; % Taekyu Shin #18
tones_fromlab=preprocessedY_fromlab/maxluminance; % Taekyu Shin #19
tones=reshape(tones_fromlab,frame);
tones2=reshape(tones_fromrgb,frame); % Taekyu Shin #20

if nargout>1
    recolor=image*colorrevert;
    recolor=cat(3,reshape(recolor(:,1),frame),reshape(recolor(:,2),frame),reshape(recolor(:,3),frame));
    recolor=min(max(recolor,0),1);
end;








function r = quantiles(x,q)
% Usage:   Finds quantiles q of data x
% Assumes: quantiles 0 <= q <= 1
% Example: r=quantiles(x,[0; 0.5; 1]);
%           xmin=r(1); xmedian=r(2); xmax=r(3)

tolerance=100*eps;
q=q(:);
x=sort(x);
n=size(x,1);
k=size(x,2);
e=1/(2*n);
q=max(e,min(1-e,q));
q=n*q+0.5;
p=min(n-1,floor(q));
q=q-p;
q(find(tolerance>q))=0;
q(find(q>(1-tolerance)))=1; 
q=repmat(q,[1,k]);
r=(1-q).*x(p,:)+q.*x(p+1,:);