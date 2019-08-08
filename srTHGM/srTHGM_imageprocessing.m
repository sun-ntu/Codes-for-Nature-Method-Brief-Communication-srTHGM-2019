clearvars -except Foldername subFolder calValue;
warning('off', 'Images:initSize:adjustingMag');

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

yellow = gray(2^16);
yellow(:,3) = 0;

white = gray(2^16);

%%

allExp=dir(fullfile([Foldername '\demodata\' subFolder '\test'],'*.tif'));

disp(' ');
disp('Imaging processing...');
for i = 1:numel(allExp)
    disp([allExp(i).name ' _ ' num2str(i)]) 
end
disp(' ');

%%

cmap = jet(256);
cmap(1,:) = 0;
cmap2 = jet(2^16);
cmap2(1,:) = 0;

if(strcmp(subFolder,'human')==1)
    par = [416.53,417.7,389.03,100,4,1.1,1];
else
    par = [342.4416,355.84,339.6952,20,5,0.9,4];
end

for Exp = [4]
    
    filename = [Foldername '\demodata\' subFolder '\test\' allExp(Exp).name];
    name = allExp(Exp).name;
    name = name(1:end-4);
    
    At = double(imread(filename,3));
    Bt = double(imread(filename,2));
    Ct = double(imread(filename,1));
    Dt = double(imread(filename,4));

    Pixel = length(At);
    iP = Pixel^2;
    
    pp = 1024/Pixel;
    
    K = ones(Pixel,Pixel);

    BgA = par(1);
    BgB = par(2);
    BgC = par(3);
    BgD = mean([BgA,BgB,BgC]);

    A = (At-BgA);
    B = (Bt-BgB);
    C = (Ct-BgC);
    D = (Dt-BgD);
    
    
    CA = C./A.*(A>0 & C>0);
    
    %% ROI for RBCs and collagen

    Mroi = double(imread([Foldername '\demodata\' subFolder '\test\' allExp(Exp-2).name]));
    Mroi = double(Mroi==0);
    
    Croi = double(imread([Foldername '\demodata\' subFolder '\test\' allExp(Exp-3).name]));
    Croi = double(Croi==0);
    
    rpp = 8/pp;
    
    roi = Mroi;
    fxr = find(max(roi));
    fyr = find(max(roi,[],2));
    
    lr = max([fxr(numel(fxr))-fxr(1),fyr(numel(fyr))-fyr(1)]);
    cxr = fix((fxr(1)+fxr(numel(fxr)))/2);
    cyr = fix((fyr(1)+fyr(numel(fyr)))/2);
        
    if(lr>Pixel-40-2*rpp)
        fyr = 21+rpp:Pixel-20-rpp;
        fxr = 21+rpp:Pixel-20-rpp;
    else
        if((cyr-fix(lr/2)-20-rpp)<1)
        fyr = fyr+21+rpp+fix(lr/2)-cyr;
        end
        if((cyr+fix(lr/2)+20+rpp)>Pixel)
            fyr = fyr+Pixel-21-rpp-fix(lr/2)-cyr;
        end
        if((cxr-fix(lr/2)-20-rpp)<1)
            fxr = fxr+21+rpp+fix(lr/2)-cxr;
        end
        if((cxr+fix(lr/2)+20+rpp)>Pixel)
            fxr = fxr+Pixel-21-rpp-fix(lr/2)-cxr;
        end
    end
    
    
    lr = max([fxr(numel(fxr))-fxr(1),fyr(numel(fyr))-fyr(1)]);
    cxr = fix((fxr(1)+fxr(numel(fxr)))/2);
    cyr = fix((fyr(1)+fyr(numel(fyr)))/2);
    
    rois = zeros(size(A));
    
    rois(cyr-fix(lr/2)-par(4)-rpp:cyr+fix(lr/2)+par(4)+rpp,cxr-fix(lr/2)-par(4)-rpp:cxr+fix(lr/2)+par(4)+rpp) = 1;
    rois(cyr-fix(lr/2)-par(4):cyr+fix(lr/2)+par(4),cxr-fix(lr/2)-par(4):cxr+fix(lr/2)+par(4)) = 0;

    rr1 = cyr-fix(lr/2)-par(4):cyr+fix(lr/2)+par(4);
    rr2 = cxr-fix(lr/2)-par(4):cxr+fix(lr/2)+par(4);

    Rois = ind2rgb(uint16(rois*2^16),white);

    InvRois = repmat((1-rois),1,1,3);
    
    %%  points selectd based on histogram cluster 
    
    Mvalue = 16383;
    internal = 40;
    peakint = 60;
    

    [ROI_A,ca] = hist(At(1:iP),0:internal:fix(Mvalue));
    [ROI_C,cc] = hist(Ct(1:iP),0:internal:fix(Mvalue));
    
    ca = ca-BgA;
    cc = cc-BgC;

    load([Foldername '\demodata\' subFolder '\test\' name '.mat']);     % icap, icbp, iccp
    
    ict = 5;
    
    if1a = max([icap(2)-ict,1]):icap(2)+ict;
    if1c = max([iccp(2)-ict,1]):iccp(2)+ict;
    
    f1at = fit(ca(if1a)',ROI_A(if1a)','gauss1','startpoint',[ROI_A(icap(2)),ca(icap(2)),50]);
    f1ct = fit(cc(if1c)',ROI_C(if1c)','gauss1','startpoint',[ROI_C(iccp(2)),cc(iccp(2)),50]);
    
    roa = (A>f1at.b1-2*f1at.c1) & (A>5*BgA);
    roc = (C>f1ct.b1-2*f1ct.c1) & (C>5*BgC);
    
    ro = roa & roc;
    
    Croi_T = Croi.*ro;
    fC = find(Croi_T);

    ip = iMp;
    rMA = (A>ip(1)-par(5)*f1at.c1) & (A<ip(1)+par(5)*f1at.c1);
    rMC = (C>ip(2)-par(5)*f1ct.c1) & (C<ip(2)+par(5)*f1ct.c1);

    ip = iCp;
    rCA = (A>ip(1)-par(5)*f1at.c1) & (A<ip(1)+par(5)*f1at.c1);
    rCC = (C>ip(2)-par(5)*f1ct.c1) & (C<ip(2)+par(5)*f1ct.c1);

    
    rM = double(ro & rMA & rMC & Mroi);
    irM = find((ro & rMA & rMC & Mroi)>0);
    
    %% mean CA ratio of Collagen for callibration
    
    ML = 5;
    rCAC = ro.*rCA.*rCC.*Croi;
    irCAC = find(rCAC>0 & CA<ML);
    mC = mean(CA(irCAC));
    
    %% THG and SHG image
    mmI = par(6)*mean(CA(fC));
    MTHG = ind2rgb(uint16(A*sqrt(mmI)),blue)+ind2rgb(uint16(B),green)+ind2rgb(uint16(C/sqrt(mmI)),red);
    MSHG = ind2rgb(uint16(D/par(7)),seagreen);
    
    MTHGb = MTHG;
    MTHGb(:,:,1) = conv2(MTHG(:,:,1),ones(2,2),'same')/4;
    MTHGb(:,:,2) = conv2(MTHG(:,:,2),ones(2,2),'same')/4;
    MTHGb(:,:,3) = conv2(MTHG(:,:,3),ones(2,2),'same')/4;
    
    MSHGb = MSHG;
    MSHGb(:,:,1) = conv2(MSHG(:,:,1),ones(2,2),'same')/4;
    MSHGb(:,:,2) = conv2(MSHG(:,:,2),ones(2,2),'same')/4;
    MSHGb(:,:,3) = conv2(MSHG(:,:,3),ones(2,2),'same')/4;
    
    figure('Name','THG image with ROI')
    imshow(MTHGb*3+Rois)
    
    figure('Name','SHG image')
    imshow(MSHGb*3)

    %% Oxygenation Saturation Calculation
    
    sT = (calValue(1)/calValue(3)-CA/mC)/(calValue(1)/calValue(3)-calValue(2)/calValue(3));
    sT = sT.*rM;
    
    d = 7;
    [sTimb,sTim] = blurimN(sT,rM,Mroi,d,cmap2);
    
    figure('Name','Subcellular RBC Oxygenation Saturation')
    imshow(sTimb(rr1,rr2,:),cmap2);
    cbh = colorbar;
    set(cbh,'YTick',0:0.2:1)
    set(cbh,'YTickLabel',{'0%','20%','40%','60%','80%','100%'})
    
    numF = numel(imfinfo([Foldername '\demodata\' subFolder '\test\' allExp(Exp-1).name]));
    T = zeros(Pixel);
    niF = zeros(1,numF);
    nmF = zeros(1,numF);
    NNsT = zeros(Pixel);
    
    for jF = 1:numF
        F = double(imread([Foldername '\demodata\' subFolder '\test\' allExp(Exp-1).name],jF));
        F = double(F==0);
        iF = find((F & ro & rMA & rMC)==1);
        mF = mean(CA(iF));
        oF = (calValue(1)/calValue(3)-mF/mC)/(calValue(1)/calValue(3)-calValue(2)/calValue(3));
        T = T.*(F==0)+F.*oF;
    end   

    Tim = ind2rgb(uint16(T*2^16),cmap2);

    figure('Name','Cellular RBC Oxygenation Saturation')
    imshow(Tim(rr1,rr2,:),cmap2);
    cbh = colorbar;
    set(cbh,'YTick',0:0.2:1)
    set(cbh,'YTickLabel',{'0%','20%','40%','60%','80%','100%'})

end
