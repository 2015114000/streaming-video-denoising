% demo

%% clear variables 
%
clear; clc; 
warning('off', 'all');
addpath(genpath('toolbox'));
addpath(genpath('external'));
%% read image and load parameters
imName1 = 'E:\博三第二学期\代码\PTSTRC-main\data\VISO dataset\mot\plane\043\img\000082.jpg';
imName2 = 'E:\博三第二学期\代码\PTSTRC-main\data\VISO dataset\mot\plane\043\img\000102.jpg';
im1 = imresize(double(imread(imName1)), 1);
im2 = imresize(double(imread(imName2)), 1);
aa=im1(:,:,2);
bb=im2(:,:,2);

sigma=20;
aa = aa+sigma * randn(size(aa));
bb = bb+sigma * randn(size(bb));


% compute optical flow

prm={'smooth',1,'radius',10,'alpha',20,'nIter',250,'type'};

% tic, [Vx,Vy]=opticalFlow(aa,bb,'smooth',1,'radius',10,'type','LK'); toc
tic, [Vx,Vy,reliab]=opticalFlow(aa,bb,prm{:},'HS'); toc
% tic, [Vx,Vy]=opticalFlow(aa,bb,prm{:},'SD','minScale',1);  toc


I1=imtransform2(aa,[],'vs',-Vx,'us',-Vy,'pad','replicate');


subplot(1,2,1);   % 1 行 3 列，第 1 个子图
imshowpair(uint8(aa), uint8(bb)); title('Image Overlay');
subplot(1,2,2);   % 1 行 3 列，第 1 个子图
imshowpair(uint8(I1), uint8(bb)); title('Image Overlay');


%  figure(2); im(I1); figure(2); im(I2);
% figure(5); imshowpair(uint8(I1), uint8(bb)); title('Image Overlay');

reliab_norm = (reliab - min(reliab(:))) / (max(reliab(:)) - min(reliab(:)));

% 显示图像
% figure(2); imshow(reliab_norm); title('Image reliab');


%compare
% addpath('E:\博二第一学期\视频去噪\代码\flow-code-matlab\flow-code-matlab');
% xy_gt=readFlowFile('E:\博二第一学期\视频去噪\dataset\other-gt-flow\other-gt-flow\Venus\flow10.flo');
% xx=xy_gt(:,:,1);
% yy=xy_gt(:,:,2);
% diff_xx=abs(u-xx);
% diff_yy=abs(v-yy);
% mean(mean(diff_xx))
% mean(mean(diff_yy))

%% store flow file and visualize flow 
