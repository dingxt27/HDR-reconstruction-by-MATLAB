clc
clear all;
I_right = imread('C:\Users\dingx\Vision_HW2\hw2q2\uttower_right.jpg');
I_left = imread('C:\Users\dingx\Vision_HW2\hw2q2\uttower_left.jpg');

I_R = rgb2gray(I_right);
I_L = rgb2gray(I_left);

 I_R = double(I_R);
 I_L = double(I_L);

% edge detection with Harris edge detection 
[cimR, rR, cR] = harris(I_R, 3, 1000, 3, 1);
[cimL, rL, cL] = harris(I_L, 3, 1000, 3, 1);

windowsize = 35;
% sub = zeros(windowsize^2,1);
n = round((windowsize)/2);

rc_R = [rR cR];

for i = 1:size(cR)
    if rc_R(i,1)>n-1 && rc_R(i,1)<= size(I_R,1)-n+1 && rc_R(i,2)>n-1 && rc_R(i,2)<= size(I_R,2)-n+1
        a = I_R((rc_R(i,1)-(n-1):(rc_R(i,1)+(n-1))),(rc_R(i,2)-(n-1):(rc_R(i,2)+(n-1))));
         subR(i,:) = reshape(a,windowsize^2,1); %storing window
%         a = a (:);
%         subR(i,:) = a;
    end 
    
end

rc_L = [rL cL];

for i = 1:size(cL)
    if rc_L(i,1)>n-1 && rc_L(i,1)<= size(I_L,1)-n+1 && rc_L(i,2)>n-1 && rc_L(i,2)<= size(I_L,2)-n+1
        b = I_L((rc_L(i,1)-(n-1):(rc_L(i,1)+(n-1))),(rc_L(i,2)-(n-1):(rc_L(i,2)+(n-1))));
        subL(i,:) = reshape(b,windowsize^2,1); %storing window
    end
    
end

%normalize to zero mean and unit standard deviation
for i = 1:size(subR,1)
    x = subR(i,:)-mean(subR(i,:));
    x1 = double(x);
    x2(i,:) = x1/std(double(subR(i,:)));
end

for j = 1:size(subL,1)
    y = subL(j,:)-mean(subL(j,:));
    y1 = double(y);
    y2(j,:) = y1/std(double(subL(j,:)));
end

%distance between 2 descriptor 
D = dist2(y2,x2);


for i = 1:size(D,1)
    p = min(D(i,:));
    p1(i)= p;
end


for i = 1:size(p1,2)
     [aa , bb] = min (p1);
%      bb;
     row(i) = bb;
     aa;
     distance(i) = aa;
     p1(bb) = max(p1);
end
%smallest 200 points
for i = 1:200
    [locationx(i) , locationy(i)] = find( D == distance(i));
end

for i = 1:size(locationx,2)
    match_points_L(i,:) = rc_L(locationx(i),:);
    match_points_R(i,:) = rc_R(locationy(i),:);
end
% 
Match_L = [match_points_L(:,2) match_points_L(:,1)];
Match_R = [match_points_R(:,2) match_points_R(:,1)];

  
figure; ax = axes;
showMatchedFeatures(I_L,I_R,Match_L,Match_R,'montage','Parent',ax);
title(ax, 'Candidate point matches_before');
legend(ax, 'Matched points 1','Matched points 2');%   

% m = RANSAC (Match_L,Match_R)

[M_opti,inlierPercentMax,inlierR,inlierL,average_residual] = getM(1000,Match_L,Match_R,10);

figure; ax = axes;
showMatchedFeatures(I_L,I_R,inlierL,inlierR,'montage','Parent',ax);
title(ax, 'Candidate point matches_inlier');
legend(ax, 'Matched points 1','Matched points 2');% 


% I_LSize = size(I_L);
% I_RSize = size(I_R);

% 

% for k =1:size(inLierL(1))

num1 = find(inlierL(:,1)~=0);
inlierL_t = zeros(size(num1));
inlierR_t = zeros(size(num1));
inlierL_t = inlierL(num1,:);
inlierR_t = inlierR(num1,:);

t = randsample(size(inlierL_t,1),4)
    
randR = zeros(4,2);
randL = zeros(4,2);

randR(1,:) = inlierR_t(t(1),:);
randR(2,:) = inlierR_t(t(2),:);
randR(3,:) = inlierR_t(t(3),:);
randR(4,:) = inlierR_t(t(4),:);

randL(1,:) = inlierL_t(t(1),:);
randL(2,:) = inlierL_t(t(2),:);
randL(3,:) = inlierL_t(t(3),:);
randL(4,:) = inlierL_t(t(4),:);


T = maketform('projective',[randL(1:4,1) randL(1:4,2)],[randR(1:4,1) randR(1:4,2)])
T.tdata.T

[im2t,xdataim2t,ydataim2t]=imtransform(I_R,T);
% now xdataim2t and ydataim2t store the bounds of the transformed im2
xdataout=[min(1,xdataim2t(1)) max(size(I_right,2),xdataim2t(2))];
ydataout=[min(1,ydataim2t(1)) max(size(I_right,1),ydataim2t(2))];
% let's transform both images with the computed xdata and ydata
im2t=imtransform(I_left,T,'XData',xdataout,'YData',ydataout);
im1t=imtransform(I_right,maketform('affine',eye(3)),'XData',xdataout,'YData',ydataout);

%averaging image
ims=im1t/2+im2t/2;
figure, imshow(ims)













