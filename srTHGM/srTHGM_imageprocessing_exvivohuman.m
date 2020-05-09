% Foldername = pwd;
% subFolder = 'exvivo_human';
clearvars -except Foldername subFolder;
close all;
warning('off', 'Images:initSize:adjustingMag');

Data = {};
load([Foldername '\demodata\' subFolder '\Loc.mat']);
folder = ["01", "02", "03", "04", "05", "06", "07", "08", "09"];
bgfile = ["bg01", "bg01", "bg01", "bg01", "bg01", "bg02", "bg02", "bg02", "bg02"];
name = ["V1 Oxy1","V1 Oxy2","V1 Deoxy1","V1 Deoxy2","V1 Deoxy 3","V2-DM Oxy1","V2-DM Oxy2","V3-DM Oxy1","V3-DM Oxy2"];

red = gray(2^16);
red(:,2) = 0;
red(:,3) = 0;

green = gray(2^16);
green(:,1) = 0;
green(:,3) = 0;

seagreen = gray(2^16);
seagreen(:,1) = seagreen(:,1)*0.1;
seagreen(:,3) = seagreen(:,3)*0.6;

blue = gray(2^16);
blue(:,1) = 0;
blue(:,2) = 0;

cyan = gray(2^16);
cyan(:,1) = 0;

magenta = gray(2^16);
magenta(:,2) = 0;

white = gray(2^16);

cmap2 = jet(2^16);
cmap2(1,:) = 0;

for N = 1:9
    Folder = [Foldername '\demodata\' subFolder '\' char(folder(N))];
    bgFolder = [Foldername '\demodata\' subFolder '\' char(bgfile(N))];
    
    ypoint=1024;
    bgypoint=1024;
    
    % background image
    BgAimg = LoadBin(bgFolder,1,bgypoint,1,2^16,5);
    BgBimg = LoadBin(bgFolder,2,bgypoint,1,2^16,5);
    BgCimg = LoadBin(bgFolder,3,bgypoint,1,2^16,5);
    
    BgA = mean(mean(BgAimg));
    BgB = mean(mean(BgBimg));
    BgC = mean(mean(BgCimg));
    
    A = LoadBin_avg(Folder,1,ypoint,1,2^16);
    B = LoadBin_avg(Folder,2,ypoint,1,2^16);
    C = LoadBin_avg(Folder,3,ypoint,1,2^16);
    
    % sub background
    Ab = A - BgA;
    Bb = B - BgB;
    Cb = C - BgC;
    
    Ar = reshape_Bin(Ab,1,3);    
    Br = reshape_Bin(Bb,1,3);
    Cr = reshape_Bin(Cb,1,3);
    
    % load image
    m = 1008;
    i = 1:m;
    mapping = round(m/2*(1-cos(pi/m*i)));
    mapAr = zeros(size(Ar));
    mapBr = zeros(size(Br));
    mapCr = zeros(size(Cr));

    for j = 1:m
        loc = find(mapping==j);
        if ~isempty(loc)
            mapAr(:,j) = mean(Ar(:,loc),2);
            mapBr(:,j) = mean(Br(:,loc),2);
            mapCr(:,j) = mean(Cr(:,loc),2);
        else
        end
    end

    for j = 1:m
        loc = find(mapping==j);
        if isempty(loc)
            mapAr(:,j) = (mapAr(:,j-1)+mapAr(:,j+1))/2;
            mapBr(:,j) = (mapBr(:,j-1)+mapBr(:,j+1))/2;
            mapCr(:,j) = (mapCr(:,j-1)+mapCr(:,j+1))/2;
        else
        end
    end
    Ar = mapAr; Br = mapBr; Cr = mapCr;
    
    % denoise
    mask = zeros(size(Ar));
    mask(200:800,200:800,:)=1;
    fAr =  abs(ifft2(ifftshift(fftshift(fft2(Ar)).*mask)));
    fBr =  abs(ifft2(ifftshift(fftshift(fft2(Br)).*mask)));
    fCr =  abs(ifft2(ifftshift(fftshift(fft2(Cr)).*mask)));
    
    imAr = ind2rgb(uint16((fAr*2^7*1.5*2.1)),blue);
    imBr = ind2rgb(uint16((fBr*2^7*1.5*2.1)),green);
    imCr = ind2rgb(uint16((fCr*2^7*1.5*2.1)),red);
    imTHG = imAr+imBr+imCr;
    figure('Name',name(N));
    imshow(imTHG);
    
    % Otsu's method to select RBC
    intensity = ind2rgb(uint16(Ar*2^5),red) + ind2rgb(uint16(Br*2^5),red) + ind2rgb(uint16(Cr*2^5),red);
    level = graythresh(intensity);
    BW = imbinarize(intensity(:,:,1),level);
    
    % select 20 RBCs for each image
    rbc = reshape(LOC(N,:,:),[20,4]);
    oriCA = Cr./Ar;
    oriBA = Br./Ar;
    oriCB = Cr./Br;
    data = [];
    for i = 1:size(rbc,1)
        roiA = zeros(size(Ar));
        roiB = zeros(size(Br));
        roiC = zeros(size(Cr));
   
        a = rbc(i, 1); b = rbc(i,2); c = rbc(i,3); d = rbc(i,4);

        roiA(a:b,c:d) = Ar(a:b,c:d)>0;
        roiB(a:b,c:d) = Br(a:b,c:d)>0;
        roiC(a:b,c:d) = Cr(a:b,c:d)>0;

        rois(a-1:b+1,c-1:d+1) = 1;
        rois(a+1:b-1,c+1:d-1) = 0;

        roiCA = roiA & roiC;
        roiBA = roiA & roiB;
        iotsuCA = find(BW>0 & roiCA>0);
        iotsuBA = find(BW>0 & roiBA>0);
        
        % calculate C/A and B/A value
        CCAAm = mean(reshape(oriCA(iotsuCA), 1, []));
        BBAAm = mean(reshape(oriBA(iotsuBA), 1, []));

        data = [data;CCAAm,BBAAm];
    end
    Data(N)={data};
end

%% In Vivo Part
% Ex Vivo Calibration Collagen C/A and B/A value
collagen_file = dir(fullfile([Foldername '\demodata\human\calibration'],'C*.tif'));
disp(' ');
disp('Calibration processing...');
    
sumCA = 0; numCA = 0;
sumBA = 0; numBA = 0;
sumA = 0; sumB = 0; sumC = 0;
numA = 0; numB = 0; numC = 0;
for Exp = 1:length(collagen_file)

    filename=[Foldername '\demodata\human\calibration\' collagen_file(Exp).name];

    At = double(imread(filename,3));
    Bt = double(imread(filename,2));
    Ct = double(imread(filename,1));

    BgA = 342.4416;
    BgB = 355.84;
    BgC = 339.6952;

    A = (At-BgA);
    B = (Bt-BgB);
    C = (Ct-BgC);
    
    CA = C./A;
    BA = B./A;
    ica = find(C>0 & A>0);
    iba = find(B>0 & A>0);
    ia = find(A>0);
    ib = find(B>0);
    ic = find(C>0);
    
    sumCA = sumCA + sum(CA(ica));
    numCA = numCA + numel(ica);
    sumBA = sumBA + sum(BA(iba));
    numBA = numBA + numel(iba);
    sumA = sumA + sum(A(ia));
    numA = numA + numel(ia);
    sumB = sumB + sum(B(ib));
    numB = numB + numel(ib);
    sumC = sumC + sum(C(ic));
    numC = numC + numel(ic);
end

ecCA = sumCA/numCA;
ecBA = sumBA/numBA;
ecA = sumA/numA;
ecB = sumB/numB;
ecC = sumC/numC;

%% In Vivo Human RBC value
invivo = [];
img = ["img1","img2"];
for i = 1:2
    all = dir(fullfile([Foldername '\demodata\exvivo_human\invivo\' char(img(i))],'*.tif'));
    filename = [Foldername '\demodata\exvivo_human\invivo\' char(img(i)) '\' all(4).name];
    name = all(4).name;
    name = name(1:end-4);
    
    At = double(imread(filename,3));
    Bt = double(imread(filename,2));
    Ct = double(imread(filename,1));
    Dt = double(imread(filename,4));
    par = [416.53,417.7,389.03,100,4,1.1,1];
    
    BgA = par(1);
    BgB = par(2);
    BgC = par(3);
    BgD = mean([BgA,BgB,BgC]);
    
    A = (At-BgA);
    B = (Bt-BgB);
    C = (Ct-BgC);
    D = (Dt-BgD);
    
    % ROI for RBCs and collagen region in image

    Mroi = double(imread([Foldername '\demodata\exvivo_human\invivo\' char(img(i)) '\' all(2).name]));
    Mroi = double(Mroi==0);
    
    Croi = double(imread([Foldername '\demodata\exvivo_human\invivo\' char(img(i)) '\' all(1).name]));
    Croi = double(Croi==0);
    
    CA = C./A;
    BA = B./A;
    CB = C./B;
    ca = find(C>0 & A>0 & Croi==1);
    ba = find(B>0 & A>0 & Croi==1);
    aa = find(A>0 & Croi==1);
    bb = find(B>0 & Croi==1);
    cc = find(C>0 & Croi==1);
    
    mC_CA = mean(CA(ca));
    mC_BA = mean(BA(ba));
    mCC = mean(C(cc));
    mBB = mean(B(bb));
    mAA = mean(A(aa));
    
    % calculate C/A and B/A value for each RBC
    numF = numel(imfinfo([Foldername '\demodata\exvivo_human\invivo\' char(img(i)) '\' all(3).name]));
    for jF = 1:numF
        F = double(imread([Foldername '\demodata\exvivo_human\invivo\' char(img(i)) '\' all(3).name],jF));
        F = double(F==0);
        iF_ca = find((F & A>0 & C>0)==1);
        iF_ba = find((F & A>0 & B>0)==1);

        mF_CA = mean(CA(iF_ca))/mC_CA*ecCA;
        mF_BA = mean(BA(iF_ba))/mC_BA*ecBA;
        invivo = [invivo; mF_CA, mF_BA];
    end
    Pixel = length(At);
    iP = Pixel^2;
    pp = 1024/Pixel;
    rpp = 8/pp;
    
    roi = Mroi;
    fxr = find(max(roi));
    fyr = find(max(roi,[],2));
    
    lr = max([fxr(numel(fxr))-fxr(1),fyr(numel(fyr))-fyr(1)]);
    cxr = fix((fxr(1)+fxr(numel(fxr)))/2);
    cyr = fix((fyr(1)+fyr(numel(fyr)))/2);
        
    if((cyr-fix(lr/2)-20-rpp)<1)
        fyr(1) = 21+rpp;
    end
    if((cyr+fix(lr/2)+20+rpp)>Pixel)
        fyr(numel(fyr)) = Pixel-20-rpp;
    end
    if((cxr-fix(lr/2)-20-rpp)<1)
        fxr(1) = 21+rpp;
    end
    if((cxr+fix(lr/2)+20+rpp)>Pixel)
        fxr(numel(fxr)) = Pixel-20-rpp;
    end
    
    lr = max([fxr(numel(fxr))-fxr(1),fyr(numel(fyr))-fyr(1)]);
    cxr = fix((fxr(1)+fxr(numel(fxr)))/2);
    cyr = fix((fyr(1)+fyr(numel(fyr)))/2);
    
    rb = 100;
    rr1 = cyr-fix(lr/2)-rb:cyr+fix(lr/2)+rb;
    rr2 = cxr-fix(lr/2)-rb:cxr+fix(lr/2)+rb;
    
    MTHG = ind2rgb(uint16(A*ecA/mAA),blue)+ind2rgb(uint16(B*ecB/mBB),green)+ind2rgb(uint16(C*ecC/mCC),red);
    MTHGb = MTHG;
    MTHGb(:,:,1) = conv2(MTHG(:,:,1),ones(2,2),'same')/4;
    MTHGb(:,:,2) = conv2(MTHG(:,:,2),ones(2,2),'same')/4;
    MTHGb(:,:,3) = conv2(MTHG(:,:,3),ones(2,2),'same')/4;
    
    figure('Name',['in vivo human RBC calibrated image' num2str(i)]);
    imshow(MTHGb(rr1,rr2,:)*5);
    
    figure('Name',['in vivo human RBC calibrated image' num2str(i) ' (only hemoglobin)']);
    MTHGb_selected = MTHGb.*Mroi;
    imshow(MTHGb_selected(rr1,rr2,:)*5);
end
Data(10) = {invivo};
%% P-Value and Box-Plot
data_CA = zeros(9,20);
data_BA = zeros(9,20);
for i=1:9
    data_CA(i,:) = Data{1,i}(:,1);
    data_BA(i,:) = Data{1,i}(:,2);
end
new_Data_CA = {};
new_Data_BA = {};

new_Data_CA(1) = {cat(2,data_CA(1,:),data_CA(2,:))};
new_Data_CA(2) = {cat(2,data_CA(3,:),data_CA(4,:),data_CA(5,:))};
new_Data_CA(3) = {cat(2,data_CA(6,:),data_CA(7,:))};
new_Data_CA(4) = {cat(2,data_CA(8,:),data_CA(9,:))};
new_Data_CA(5) = {reshape(Data{1,10}(:,1),[1,4])};

new_Data_BA(1) = {cat(2,data_BA(1,:),data_BA(2,:))};
new_Data_BA(2) = {cat(2,data_BA(3,:),data_BA(4,:),data_BA(5,:))};
new_Data_BA(3) = {cat(2,data_BA(6,:),data_BA(7,:))};
new_Data_BA(4) = {cat(2,data_BA(8,:),data_BA(9,:))};
new_Data_BA(5) = {reshape(Data{1,10}(:,2),[1,4])};

P = zeros(2,5,5);
for i=1:5
    for j=1:5
        [h,p] = ttest2(new_Data_CA{1,i},new_Data_CA{1,j}); % CA
        P(1,i,j) = p;
        [h,p] = ttest2(new_Data_BA{1,i},new_Data_BA{1,j}); % BA
        P(2,i,j) = p;
    end
end

figure('Name','C/A Value Box-Plot');
boxplot([new_Data_CA{1,1} new_Data_CA{1,2} new_Data_CA{1,3} new_Data_CA{1,4} new_Data_CA{1,5}],...
    [zeros(1,40),ones(1,60),2*ones(1,40),3*ones(1,40),4*ones(1,4)],...
    'Labels',{'V1 Oxy','V1 Deoxy','V2-DM Oxy','V3-DM Oxy','In Vivo'});
title('')
ylabel('C/A Value')
hold on
plot(5, new_Data_CA{1,5},'*','color','m','linewidth', 0.5)


figure('Name','B/A Value Box-Plot');
boxplot([new_Data_BA{1,1} new_Data_BA{1,2} new_Data_BA{1,3} new_Data_BA{1,4} new_Data_BA{1,5}],...
    [zeros(1,40),ones(1,60),2*ones(1,40),3*ones(1,40),4*ones(1,4)],...
    'Labels',{'V1 Oxy','V1 Deoxy','V2-DM Oxy','V3-DM Oxy','In Vivo'});
title('')
ylabel('B/A Value')
hold on
plot(5, new_Data_BA{1,5},'*','color','m','linewidth', 0.5)
