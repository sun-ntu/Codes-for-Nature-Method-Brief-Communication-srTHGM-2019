clearvars -except Foldername;

format long;

datafile = [Foldername '\demodata'];

% Load parameters
load([Foldername '\frequencybase.mat']);  %  c0, de=10, dde=2, s=1e10, f(Hz), ln(nm), ef(Iw datapoint)

tu = 1/(max(f)-min(f));
t = 10^15*(0:tu:(s)^-1);    % frequency sampling
t = t-t(round(numel(t)/2));
w = 2*pi.*f;
n=0;

% Load data
[le,Iwe] = textread([datafile '\20171218_c_5270_36_18_586mW_1268_95']);
[lf,Iwf,pf0] = textread([datafile '\Pulse_Freq_Domain 12_18_time_13_23_894_0.txt'],'%*n %n %n %n');
[tf,Itf,pt0] = textread([datafile '\Pulse_Time_Domain 12_18_time_13_23_894_0.txt'],'%n %n %n');

% Iwec = Iwe';
le = le';
Iwe = Iwe';

temp = fft(Iwe);
temp(300:end-300) =0;
Iwe = abs(ifft(temp));

fe = c0*10^9 ./le/dde;
fs = roundn(fe,de)*dde;

%%
ff = c0*10^9./lf/dde;
ff = roundn(ff,de)*dde;

if (mean(pf0)<0)
    pf0 = -pf0;
end
pf0 = pf0-min(pf0);

pfl = fft(pf0);
pfl(30:226) = 0;
pf0 = abs(ifft(pfl));

pp = find(lf>1090 & lf<1510);
%%

temp = fft(pt0);
temp(75:end-75) =0;
pt0 = abs(ifft(temp));
pt0 = pt0-min(pt0);

%%

sm = 2;

Iw0 = interp1(fs,Iwe,f,'spline',0);

y1 = peak(1,1,ln,Iw0/max(Iw0));
[y2,xc] = min(abs(ln-y1));

pf = interp1(ff(pp(1):sm:pp(end)),pf0(pp(1):sm:pp(end)),f,'pchip',0);
pt = interp1(tf,pt0,t,'pchip',0);

if(n==1)
    pf = pf*-1+10;
end

%%

% Iw0

l = 1150;
u = 1400;
[cm2,cu] = min(abs(ln-u));
[cm3,cl] = min(abs(ln-l));

% plot figure

figure('Name','Measured Spectrum of Laser with Phase')
[ha,h1,h2]=plotyy(ln,abs(Iw0)/max(Iw0),ln(cu:cl),pf(cu:cl));
set(h1,'color',[0.8,0,0],'LineWidth',5)
set(h2,'color',[0.9,0.5,0],'LineWidth',5)
xlim(ha(1),[1200,1332]);
xlim(ha(2),[1200,1332]);
ylim(ha(1),[0,1.1]);
ylim(ha(2),[0,11]);
set(ha,{'ycolor'},{[0.8,0,0];[0.9,0.5,0]})
set(ha(1),'fontsize',15,'linewidth',2,'XTick',1215:15:1320, 'XAxisLocation', 'top','YTick',0.2:0.2:1,'FontWeight','bold')
set(ha(2),'fontsize',15,'linewidth',2,'XTick',1215:15:1320, 'XAxisLocation', 'top','YTick',2:2:10,'FontWeight','bold')
xlabel('Fundamental wavelength (nm)','fontsize',20,'FontWeight','bold');
ylabel(ha(1),'Relative intensity (a.u.)','color','k','fontsize',20,'FontWeight','bold');
ylabel(ha(2),'Phase (rad)','color','k','fontsize',20,'FontWeight','bold');
box on;
set(gcf,'color','w');
set(gcf,'Position', [200 100 805 550]);
pbaspect(ha(1),[1.7 1 1])
pbaspect(ha(2),[1.7 1 1])
