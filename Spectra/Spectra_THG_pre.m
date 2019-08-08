% close all;
clearvars -except Foldername;
clc;

allExp=dir(fullfile(Foldername, 'demodata', '*.sif'));

AttCell = cell([1 10]);
Attribute = struct;
AttributeHeadings = { 'sample','test','num','exP','t','slit','zDist','x','y','lpy','xM','yM','title'};

k = 1;
t = 0;

disp('THG spectra preprocessing...');
for j = 1:length(allExp)

% Load spectra files
info = sifreadnk([Foldername '\demodata\' allExp(j).name]);
x = info.axisWavelength; 
i = length(x);

% Filename example
% 20171107 Water 0 180mW 200ms 50 14654 10 u 1.sif
           % 1   2  3     4    5    6    7 8 9
AttCell = textscan(allExp(j).name,'%*n %s %n %n%*s %n%*s %n %n %n %c %n.sif');
                                  %     1  2  3     4     5  6  7  8  9         

% Read spectra files and save parameters for analysis                                  
if (AttCell{2}==0)
    
    t = t+1;
    
else
    if ( AttCell{8}=='u' && zN==(AttCell{6}-AttCell{7}/10) )
        
    elseif ( AttCell{8}=='d' && zN==(AttCell{6}+AttCell{7}/10) )
        
    else
        En = (AttCell{6}-zN)/AttCell{7}*10;
        
        for j3 = k:k+En-1
            Attribute(j3).sample = AttCell{1};
            Attribute(j3).test = t;
            Attribute(j3).num = -1;
            Attribute(j3).exP = AttCell{3};      % mW
            Attribute(j3).t = AttCell{4};        % ms
            Attribute(j3).slit = AttCell{5};     % um
                if (AttCell{8}=='u')
                    Attribute(j3).zDist = zN+AttCell{7}*(j3-k+1)/10;  % um
                elseif (AttCell{8}=='d')
                    Attribute(j3).zDist = zN-AttCell{7}*(j3-k+1)/10;
                else
                    Attribute(j3).zDist = zN;
                end

            Attribute(j3).x = x;               % nm
            Attribute(j3).y = 260+10*randn(size(y));  
                temp = ifft(Attribute(j3).y);
                temp(round(0.15*i):round(0.85*i)) = 0; 
            Attribute(j3).lpy = abs(fft(temp));
            Attribute(j3).xM = 0;
            Attribute(j3).yM = 0;
            Attribute(j3).title = [char(AttCell{1}) '_' num2str(Attribute(j3).test) ',' num2str(Attribute(j3).zDist) 'um'];           
        end
        
        k = k+En;
    end
end


for j2 = k:(k+info.kineticLength-1)
    
    y = double(info.imageData(:,:,(j2-k+1)))';
    
    [MaxValue,Maxposition] = max(y(201:700));
    Maxposition = Maxposition+200;
    
    Attribute(j2).sample = AttCell{1};
    Attribute(j2).test = t;
    Attribute(j2).num = AttCell{2};
    Attribute(j2).exP = AttCell{3};      % mW
    Attribute(j2).t = AttCell{4};        % ms
    Attribute(j2).slit = AttCell{5};     % um
        if (AttCell{8}=='u')
            Attribute(j2).zDist = AttCell{6}+AttCell{7}*(j2-k+1)/10;  % um
        elseif (AttCell{8}=='d')
            Attribute(j2).zDist = AttCell{6}-AttCell{7}*(j2-k+1)/10;
        else
            Attribute(j2).zDist = AttCell{6};
        end

    Attribute(j2).x = x;               % nm
    Attribute(j2).y = y;  
        temp = ifft(y);
        temp(round(0.15*i):round(0.85*i)) = 0;
    Attribute(j2).lpy = abs(fft(temp));
    Attribute(j2).xM = x(Maxposition);
    Attribute(j2).yM = MaxValue;
    Attribute(j2).title = [char(AttCell{1}) '_' num2str(Attribute(j2).test) ',' num2str(Attribute(j2).zDist) 'um'];
end

zN = Attribute(j2).zDist;
k = k+info.kineticLength;
disp(['(' num2str(Attribute(j2).test) ')  ' char(AttCell{1})]); % file order and type of sample
end

% Save spectra data for processing
save([Foldername '\demodata\Attribute_org.mat'],'Attribute','AttributeHeadings','-mat');


