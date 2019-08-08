clear all;
close all;
clc

%%
% colormap setting
red = gray(2^16);
red(:,2) = 0;
red(:,3) = 0;

green = gray(2^16);
green(:,1) = 0;
green(:,3) = 0;

blue = gray(2^16);
blue(:,1) = 0;
blue(:,2) = 0;

magenta = gray(2^16);
magenta(:,2) = 0;

%%
Foldername = pwd;

% Load ROI masks for three tissues: 
% Melanin, Red Blood Cell, and Sebaceous Gland
MASK_Me = imresize(imread([Foldername '\demodata\Me.tif']),0.5);
MASK_RBC = imresize(imread([Foldername '\demodata\RBC.tif']),0.5);
MASK_SG = imresize(imread([Foldername '\demodata\SG.tif']),0.5);

% Load THG image with 3 channel
THGimage = double(imread([Foldername '\demodata\Tmice5_170_512_avg15_4CH_a7z4-1_650670590.tif']));
At = THGimage(:,:,1);
Bt = THGimage(:,:,2);
Ct = THGimage(:,:,3);

% Background values
BgA = 342.4416;
BgB = 355.84;
BgC = 339.6952;

% Substract background
A = At-BgA;
B = Bt-BgB;
C = Ct-BgC;

%%
% Conventional THG intensity image
S = A+B+C;
Sb = ind2rgb(uint16(S*3),magenta);

% denoise process
Sb(:,:,1) = conv2(Sb(:,:,1),ones(2,2),'same')/4;
Sb(:,:,2) = conv2(Sb(:,:,2),ones(2,2),'same')/4;
Sb(:,:,3) = conv2(Sb(:,:,3),ones(2,2),'same')/4;
figure('Name','Conventional THG intensity image')
imshow(Sb);

% Pseudo-colored THG image
Tim = ind2rgb(uint16(A*12),blue)+ind2rgb(uint16(B*12),green)+ind2rgb(uint16(C*12),red);

% denoise process
Timb(:,:,1) = conv2(Tim(:,:,1),ones(2,2),'same')/4;
Timb(:,:,2) = conv2(Tim(:,:,2),ones(2,2),'same')/4;
Timb(:,:,3) = conv2(Tim(:,:,3),ones(2,2),'same')/4;

figure('Name','Pseudo-colored THG image')
imshow(Timb);

% Pseudo-colored THG images of separate channels
TimA = ind2rgb(uint16(A*2^3),blue);
TimB = ind2rgb(uint16(B*2^3),green);
TimC = ind2rgb(uint16(C*2^3),red);

% denoise process
TimAb(:,:,1) = conv2(TimA(:,:,1),ones(2,2),'same')/4;
TimAb(:,:,2) = conv2(TimA(:,:,2),ones(2,2),'same')/4;
TimAb(:,:,3) = conv2(TimA(:,:,3),ones(2,2),'same')/4;
TimBb(:,:,1) = conv2(TimB(:,:,1),ones(2,2),'same')/4;
TimBb(:,:,2) = conv2(TimB(:,:,2),ones(2,2),'same')/4;
TimBb(:,:,3) = conv2(TimB(:,:,3),ones(2,2),'same')/4;
TimCb(:,:,1) = conv2(TimC(:,:,1),ones(2,2),'same')/4;
TimCb(:,:,2) = conv2(TimC(:,:,2),ones(2,2),'same')/4;
TimCb(:,:,3) = conv2(TimC(:,:,3),ones(2,2),'same')/4;

figure('Name','Pseudo-colored THG image for Channel A')
imshow(TimAb);
figure('Name','Pseudo-colored THG image for Channel B')
imshow(TimBb);
figure('Name','Pseudo-colored THG image for Channel C')
imshow(TimCb);

%%
% THG strength distribution for 
% Load ROI of Melanin, Red Blood Cell, and Sebaceous Gland
roiMe = double(MASK_Me==0);
roiRBC = double(MASK_RBC==0);
roiSG = double(MASK_SG==0);

% Select SNR > 5
roiA =  (A>5*BgA);
roiB =  (B>5*BgB);
roiC =  (C>5*BgC);
roi = roiA & roiC;

% Select both SNR > 5 and ROI
irMe = find((roiMe & roi)>0);
irRBC = find((roiRBC & roi)>0);
irSG = find((roiSG & roi)>0);

% THG strength distribution calculation
R = zeros(3,3);

ratio = 1;
ratio1 = 1;
ratio2 = 1;
R(1,1) = ratio*mean(A(irMe));
R(1,2) = ratio2*mean(B(irMe));
R(1,3) = ratio1*mean(C(irMe));
R(2,1) = ratio*mean(A(irRBC));
R(2,2) = ratio2*mean(B(irRBC));
R(2,3) = ratio1*mean(C(irRBC));
R(3,1) = ratio*mean(A(irSG));
R(3,2) = ratio2*mean(B(irSG));
R(3,3) = ratio1*mean(C(irSG));
R = R/6000; % Normalized

RS = zeros(1,9);
RS(1) = sem(A(irMe));
RS(2) = sem(B(irMe));
RS(3) = sem(C(irMe));
RS(4) = sem(A(irRBC));
RS(5) = sem(B(irRBC));
RS(6) = sem(C(irRBC));
RS(7) = sem(A(irSG));
RS(8) = sem(B(irSG));
RS(9) = sem(C(irSG));
RS = RS/6000; % Normalized

figure('Name','THG strength distribution for Mice')
hold on
b = bar(R);
errorbar([0.78,1,1.22,    1.78,2,2.22   2.78,3,3.22],R([1,4,7,2,5,8,3,6,9]),RS,'k.','linewidth',3)

set(b,'Facecolor',[1,1,1])
set(b(1),'Edgecolor',[0.2,0.35,1])
set(b(2),'Edgecolor',[0.35,0.85,0.55])
set(b(3),'Edgecolor',[1,0.15,0.1])
set(b,'linewidth',3)
ylabel('Relative intensity (a.u.)','fontsize',20,'FontWeight','bold');
legend('A','B','C','FontWeight','bold')
legend boxoff
box off;
xlim([0.5 3.8])
ylim([0 1.5])
names = {'Melanocyte'; 'Erythrocyte'; 'Sebaceous gland'};
set(gca,'XTick',1:3,'XTicklabel',names)
set(gca,'YTick',0:0.2:1.2,'FontWeight','bold')
set(gca,'fontsize',15)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])

display(['SG: ' num2str(mean(A(irSG)),'%2.0f') '  ' num2str(mean(B(irSG)),'%2.0f') '  ' num2str(mean(C(irSG)),'%2.0f')]);
display(['RBC:' num2str(mean(A(irRBC)),'%2.0f') '  ' num2str(mean(B(irRBC)),'%2.0f') '  ' num2str(mean(C(irRBC)),'%2.0f')]);
display(['Me: ' num2str(mean(A(irMe)),'%2.0f') '  ' num2str(mean(B(irMe)),'%2.0f') '  ' num2str(mean(C(irMe)),'%2.0f')]);
display(['num of point selected: ' num2str(numel(irSG),'%2.0f') '  ' num2str(numel(irRBC),'%2.0f') '  ' num2str(numel(irMe),'%2.0f')]);