clearvars -except Foldername;
clc;

format long;
%%

load([Foldername '\AC.mat']);

% z-location of intensity peak (or selected) occured interface for different samples,
% value 0 for automatically detected by program
Ip  = [ 0,  0,    0,  0,   36, 35,    0,  0,    0,  0,    0,  0,    0,  0];
    %   1,  2,    3,  4,    5,  6,    7,  8,    9, 10,   11, 12,   13, 14,   15, 16,   17, 18

% File number set according to sample orders
spec = [ 1,  2,    7,  8,    3,  4,    5,  6,   11, 12,    9, 10,   13,14];
     %   1,  2,    3,  4,    5,  6,    7,  8,    9, 10,   11, 12,   13, 14,   15, 16,   17, 18

% Load sample file attributes
load([Foldername '\demodata\Attribute_org.mat']); 
isp =length(spec); 

%%

iT = isp;
Tag = {'Fused silica_1','Fused silica_2','Water_1','Water_2','Hemoglobin_1','Hemoglobin_2','Hemoglobin/O_2_1','Hemoglobin/O_2_2','Melanin_1','Melanin_2','Linoleic acid_1','Linoleic acid_2','Oleic acid_1','Oleic acid_2'};
Tag2 = {'Fused silica','Water','Deoxy-hemoglobin','Oxy-hemoglobin','Melanin','Linoleic acid','Oleic acid'};
Tag2a = Tag2;

chsv = hsv(8);
ctmap = zeros(isp/2,3);
ctmap(1:isp/2,:) = chsv([1,2,3,5,6,7,8],:);
ctmap(:,2) = 0.7*ctmap(:,2);

%%
load([Foldername '\lpfilter.mat']);
load([Foldername '\frequencybase.mat']);  %   c0, de=10, dde=2, s=1e10, f(Hz),  ln(nm),   ef(Iw datapoint),
lu = ln/1000; %um
i = length(f);

x = Attribute(1).x;
fexp = c0*10^9./x/dde;
fexpc = roundn(fexp,de)*dde;

iy = length(Attribute(1).y);

%%

it = max([Attribute.test]);     %indber of test

yInf1_lp = zeros(it,i);

%%

check = 1;

PT=1; % plot spectra or not
ST=1; % save spectra data or not

disp('THG spectra processing...');
for q1 = 1:14
    
    % find peak (interface)
    stack = find([Attribute.test]==q1);
    name = char(Attribute(stack(q1)).sample);
    disp(['(' num2str(q1) ')  ' name]);
    is = length(stack);
    is2 = round(is/2);
    [~,IndexSorted] = sort([Attribute(stack).zDist]); 
    stack = stack(IndexSorted);
    
    y = zeros(is,iy);
    y_lp = y;
    y_nn = y;
    y_lpnn = y;
    z = 25-[Attribute(stack).zDist]/1000;
    
    for q2 = 1:is
        tAC = Attribute(stack(q2)).t/100;
        y(q2,:) = Attribute(stack(q2)).y;
        y(q2,:) = (y(q2,:)-min(y(q2,:)))./AC(tAC).scale;
        y_lp(q2,:) = Attribute(stack(q2)).lpy;
        y_lp(q2,:) = (y_lp(q2,:)-min(y_lp(q2,:)))/AC(tAC).scale_m;
        y_nn(q2,:) = (y(q2,:)-min(y(q2,:)));
        y_nn(q2,:) = y_nn(q2,:)./max(y_nn(q2,:));
        y_lpnn(q2,:) = (y_lp(q2,:)-min(y_lp(q2,:)));
        y_lpnn(q2,:) = y_lpnn(q2,:)./max(y_lpnn(q2,:));
    end
    
    y_n = y/max(max(y));
    y_lpn = y_lp/max(max(y_lp));
    
    % find peak intensity (interface)
    r1 = 1:60;     % first interface possible range
    r2 = 61:is-59;
    r3 = is-59:is;
    if(Ip(q1)==0)
        [yInterf1,zInterf1] = max([Attribute(stack(r1)).yM]);
        zInterf1 = zInterf1+r1(1)-1;
    else
        yInterf1 = Attribute(Ip(q1)).yM;
        zInterf1 = Ip(q1);
    end
    
    yInf1_lp(q1,:) = interp1(fexpc,y_lp(zInterf1,:),f,'spline',0);

end
%%

sample = zeros(1,it);
for pri = 1:it
    temps = find([Attribute.test]==pri);
    sample(pri) = temps(1);
end

%%

if (check==1)    
    %%
    % Spectra averaging and normalized    
    yInf1_lpn = yInf1_lp/max((yInf1_lp(spec(1),:)+yInf1_lp(spec(2),:))/2);

    avg = zeros(isp/2,i);
    avg_n = avg;
    avg_nn = avg;

    aR = zeros(isp/2,i);
    aG = zeros(isp/2,i);
    aB = zeros(isp/2,i);

    for q5 = 1:isp/2
        if (spec(2*q5-1)==0)
            Tag(1,2*q5-1)={' '};
            Tag(1,2*q5)={' '};
            Tag2(1,q5)={' '};
            iT = iT-2;
            
        else
            avg(q5,:) = (yInf1_lp(spec(2*q5-1),:)+yInf1_lp(spec(2*q5),:))/2;
            avg_n(q5,:) = (yInf1_lpn(spec(2*q5-1),:)+yInf1_lpn(spec(2*q5),:))/2;
            avg_nn(q5,:) = avg_n(q5,:)./max(avg_n(q5,:));

        end
        
        aR(q5,:) = avg(q5,:).*FR;
        aG(q5,:) = avg(q5,:).*FG;
        aB(q5,:) = avg(q5,:).*FB;
    end
    
    SR = sum(aR,2);
    SG = sum(aG,2);
    SB = sum(aB,2);
    S = [SB,SG,SR];
    Sn = [SB./SG,SG./SG,SR./SG];

    %%
    % For susceptibility calculation
    R = zeros(isp/2,i);
    R_n = R;
    

    R(1,:) = ones(1,i);
    R_n(1,:) = ones(1,i);

    for q5 = 2:isp/2
        if (spec(2*q5-1)==0)
            R(q5,:) = zeros(1,i);
            R_n(q5,:) = zeros(1,i);
        else      
            R(q5,:) = avg(q5,:)./avg(1,:);
            R_n(q5,:) = R(q5,:)/max(R(q5,68130/dde:74960/dde));
        end
        

    end

    refX = (-1.6393E-9*(ln*3).^3+7.0902E-6*(ln*3).^2-10.652E-3*(ln*3)+7.2566)*1E-22; % X(w) @1260
    RX = 3;
    Xref = RX*R.*repmat(refX,[isp/2,1]);
    Xref_n = RX*R_n.*repmat(refX,[isp/2,1]);
    Xref(1,:) = Xref(1,:)/RX;
    Xref_n(1,:) = Xref_n(1,:)/RX;
    
    

    %%
    % Save spetra information
    Tt = zeros(1,iT);
    Tt2 = zeros(1,iT/2);
    q7=1;

    yInf1 = struct;
    yInf1Headings = { 'sample','lp_1','lp_2','Avg','Avg_n','Avg_nn','R','R_n'};
    

    for q5 = 1:isp/2
        if (spec(2*q5-1)==0)
            yInf1(q5).sample = 'no';
            yInf1(q5).lp_1 = zeros(1,i);
            yInf1(q5).lp_2 = zeros(1,i);
  
        else
            Tt(2*q7-1)=2*q5-1;
            Tt(2*q7)=2*q5;
            Tt2(q7) = q5;
            q7 = q7+1;
            
            yInf1(q5).sample = char(Attribute(sample(spec(2*q5-1))).sample);
            yInf1(q5).lp_1 = yInf1_lpn(spec(2*q5-1),:);
            yInf1(q5).lp_2 = yInf1_lpn(spec(2*q5),:); 

        end
        yInf1(q5).Avg = avg(q5,:);
        yInf1(q5).Avg_n = avg_n(q5,:);
        yInf1(q5).Avg_nn = avg_nn(q5,:);        
        yInf1(q5).R = R(q5,:);
        yInf1(q5).R_n = R_n(q5,:);
        yInf1(q5).S = S(q5,:);
        yInf1(q5).S_n = Sn(q5,:);
        
    end

    
    if (PT==1)
        
        Spectra_THG_plot;
        
    end
    if (ST==1)
        
        save([Foldername '\demodata\material_peak.mat'],'yInf1','Tag2a','Tag2','Tt2','spec','-mat');
        
    end
end