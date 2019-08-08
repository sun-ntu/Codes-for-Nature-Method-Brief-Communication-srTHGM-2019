clearvars -except Foldername subFolder;
warning('off', 'Images:initSize:adjustingMag');

red = gray(2^16);
red(:,2) = 0;
red(:,3) = 0;

green = gray(2^16);
green(:,1) = 0;
green(:,3) = 0;

blue = gray(2^16);
blue(:,1) = 0;
blue(:,2) = 0;

cyan = gray(2^16);
cyan(:,1) = 0;

magenta = gray(2^16);
magenta(:,2) = 0;

yellow = gray(2^16);
yellow(:,3) = 0;

%%
allExp = {};
tissue = {'Deoxy-hemo','Oxy-hemo','Collagen'};
allExp{1}=dir(fullfile([Foldername '\demodata\' subFolder '\calibration'],'D*.tif'));
allExp{2}=dir(fullfile([Foldername '\demodata\' subFolder '\calibration'],'O*.tif'));
allExp{3}=dir(fullfile([Foldername '\demodata\' subFolder '\calibration'],'C*.tif'));

disp(' ');
disp('Calibration processing...');
for i = 1:3
    for j = 1:length(allExp{i})
        disp(allExp{i}(j).name) 
    end
end
disp(' ');

%%

cmap = jet(256);
cmap(1,:) = 0;

RecS = 0;
iRS = 1;

CC=1;
calValue = zeros(3,1);

for j = 1:3
    
    sumCA = 0; numCA = 0;
    
    for Exp = 1:length(allExp{j})

        filename=[Foldername '\demodata\' subFolder '\calibration\' allExp{j}(Exp).name];
        name = allExp{j}(Exp).name;
        name = name(1:end-4);

        At = double(imread(filename,3));
        Bt = double(imread(filename,2));
        Ct = double(imread(filename,1));

        Pixel = length(At);
        iP = Pixel^2;

        pp = 1024/Pixel;

        K = ones(Pixel,Pixel);

        BgA = 342.4416;
        BgB = 355.84;
        BgC = 339.6952;

        A = (At-BgA);
        B = (Bt-BgB);
        C = (Ct-BgC);

        %%  points selectd based on histogram cluster

        Mvalue = 16383;
        internal = 10;

        [ROI_A,ca] = hist(At(1:iP),0:internal:fix(Mvalue));
        [ROI_B,cb] = hist(Bt(1:iP),0:internal:fix(Mvalue));
        [ROI_C,cc] = hist(Ct(1:iP),0:internal:fix(Mvalue));

        ca = ca-BgA;
        cb = cb-BgB;
        cc = cc-BgC;

        load([Foldername '\demodata\' subFolder '\calibration\' name '.mat']);     % icap, icbp, iccp

        ict = 10;

        f3at = fit(ca(icap(3)-ict:icap(3)+ict)',ROI_A(icap(3)-ict:icap(3)+ict)','gauss1','startpoint',[ROI_A(icap(3)),ca(icap(3)),10]);
        if icbp~=0
            f3bt = fit(cb(icbp(3)-ict:icbp(3)+ict)',ROI_B(icbp(3)-ict:icbp(3)+ict)','gauss1','startpoint',[ROI_B(icbp(3)),cb(icbp(3)),10]);
        end
        f3ct = fit(cc(iccp(3)-ict:iccp(3)+ict)',ROI_C(iccp(3)-ict:iccp(3)+ict)','gauss1','startpoint',[ROI_C(iccp(3)),cc(iccp(3)),10]);  

        c1b = 2;

        rTh = (A>5*BgA) & (C>5*BgC);

        rGa3 = (A>f3at.b1-c1b*f3at.c1 & A<f3at.b1+c1b*f3at.c1);
        if icbp~=0
            rGb3 = (B>f3bt.b1-c1b*f3bt.c1 & B<f3bt.b1+c1b*f3bt.c1);
        end
        rGc3 = (C>f3ct.b1-c1b*f3ct.c1 & C<f3ct.b1+c1b*f3ct.c1);

        rG3 = rGa3.*rGc3;


        %% plot and save point CA value

        kk = 2*sqrt(0.75);
        Tim = ind2rgb(uint16((A-BgA)*kk*2^4),blue)+ind2rgb(uint16((B-BgB)*2^4),green)+ind2rgb(uint16((C-BgC)/kk*2^4),red);

        Timb = Tim;
        Timb(:,:,1) = conv2(Tim(:,:,1),ones(2,2),'same')/4;
        Timb(:,:,2) = conv2(Tim(:,:,2),ones(2,2),'same')/4;
        Timb(:,:,3) = conv2(Tim(:,:,3),ones(2,2),'same')/4;

        figure('Name',[tissue{1,j},num2str(Exp)])
        imshow(Timb*1.8)
        
        CA = C./A.*(C>0).*(A>0);
        iF = find((rTh & rG3)==1);

        sumCA = sumCA + sum(CA(iF));
        numCA = numCA + numel(iF);
        
    end
    
    %% average CA ratio of all point selected
    
    calValue(j) = sumCA/numCA;
    disp([tissue{1,j},' CA ratio = ',num2str((sumCA/numCA),'%5.3f')])
    
end
