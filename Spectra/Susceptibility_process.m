clearvars -except Foldername;
clc;
warning('off','all');

format long;
%%
% colormap
chsv = hsv(8);
ctmap = chsv([1,2,3,5,6,7,8],:);
ctmap(:,2) = 0.7*ctmap(:,2);

%%
% load data and parameter

bs = [Foldername '\frequencybase.mat'];
load(bs);  %   c0, de=10, dde=2, s=1e10, f(Hz),  ln(nm),   ef(Iw datapoint),
lu = ln/1000; %um
i = length(f);

load([Foldername '\demodata\material_peak.mat']);
load([Foldername '\THG_sim0_real.mat']);
load([Foldername '\lpfilter.mat'])

%%

RaInf = zeros(7,i);

% Average THG Spectra
RaInf(1,:) =yInf1(1).Avg_n;
RaInf(2,:) =yInf1(2).Avg_n;
RaInf(3,:) =yInf1(3).Avg_n;
RaInf(4,:) =yInf1(4).Avg_n;
RaInf(5,:) =yInf1(5).Avg_n;
RaInf(6,:) =yInf1(6).Avg_n;
RaInf(7,:) =yInf1(7).Avg_n;

RaInf = RaInf.*(RaInf>0);

% THG Spectra 1
a1 = zeros(7,i);
a1(1,:) = yInf1(1).lp_1;
a1(2,:) = yInf1(2).lp_1;
a1(3,:) = yInf1(3).lp_1;
a1(4,:) = yInf1(4).lp_1;
a1(5,:) = yInf1(5).lp_1;
a1(6,:) = yInf1(6).lp_1;
a1(7,:) = yInf1(7).lp_1;

a1 = a1.*(a1>0);

% THG Spectra 2
a2 = zeros(7,i);
a2(1,:) = yInf1(1).lp_2;
a2(2,:) = yInf1(2).lp_2;
a2(3,:) = yInf1(3).lp_2;
a2(4,:) = yInf1(4).lp_2;
a2(5,:) = yInf1(5).lp_2;
a2(6,:) = yInf1(6).lp_2;
a2(7,:) = yInf1(7).lp_2;

a2 = a2.*(a2>0);

%%

% interpolation spacing
s = 110;

% calculated based on Fused Silica
base  = RaInf(1,:); 

%%
% Nonlinear Susceptibility calculation
Xref = (-1.6393E-9*(ln*3).^3+7.0902E-6*(ln*3).^2-10.652E-3*(ln*3)+7.2566)*1E-22;

Xexp_FSt = Xref.*sqrt(RaInf(1,:)./base.*Ib3w(1,:)./Ib3w(1,:));
Xexp_FS = interp1(ln(1:s:end),Xexp_FSt(1:s:end),ln,'spline',0);

Xexp_Wt = Xref.*sqrt(RaInf(2,:)./base.*Ib3w(1,:)./Ib3w(2,:));
Xexp_W = interp1(ln(1:s:end),Xexp_Wt(1:s:end),ln,'spline',0);

Xexp_Hbt = Xref.*sqrt(RaInf(3,:)./base.*Ib3w(1,:)./Ib3w(7,:));
Xexp_Hb = interp1(ln(1:s:end),Xexp_Hbt(1:s:end),ln,'spline',0);

Xexp_HbOt = Xref.*sqrt(RaInf(4,:)./base.*Ib3w(1,:)./Ib3w(8,:));
Xexp_HbO = interp1(ln(1:s:end),Xexp_HbOt(1:s:end),ln,'spline',0);

Xexp_Met = Xref.*sqrt(RaInf(5,:)./base.*Ib3w(1,:)./Ib3w(3,:));
Xexp_Me = interp1(ln(1:s:end),Xexp_Met(1:s:end),ln,'spline',0);

Xexp_LAt = Xref.*sqrt(RaInf(6,:)./base.*Ib3w(1,:)./Ib3w(4,:));
Xexp_LA = interp1(ln(1:s:end),Xexp_LAt(1:s:end),ln,'spline',0);

Xexp_OAt = Xref.*sqrt(RaInf(7,:)./base.*Ib3w(1,:)./Ib3w(5,:));
Xexp_OA = interp1(ln(1:s:end),Xexp_OAt(1:s:end),ln,'spline',0);

%% plot

figure('Name','Susceptibility Spectra 1') % Water; Melanin
hold on
plot(3*ln,Xexp_W*1E22,'color',ctmap(2,:),'LineWidth',5)
plot(3*ln,Xexp_Me*1E22,'color',ctmap(5,:),'LineWidth',5)
hold off
xlim([1230 1290]);
legend({'Water','Melanin'},'Location','Northeast','fontsize',12)
xlabel('Wavelength(nm)','fontsize',20,'FontWeight','bold');
ylabel('Nonlinear susceptibility (10^-^2^2 m^2V^-^2)','fontsize',20,'FontWeight','bold');
box on;
set(gca, 'YAxisLocation', 'right')
set(gca,'XTick',1236:6:1284,'FontWeight','bold')
ylim([1 9]);
set(gca,'YTick',2:1.5:8,'FontWeight','bold')
set(gca,'fontsize',15)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])

figure('Name','Susceptibility Spectra 2') % Deoxy-hemoglobin; Oxy-hemoglobin
hold on
plot(3*ln,Xexp_Hb*1E22,'color',ctmap(3,:),'LineWidth',5)
plot(3*ln,Xexp_HbO*1E22,'color',ctmap(4,:),'LineWidth',5)
hold off
xlim([1230 1290]);
legend({'Deoxy-hemoglobin','Oxy-hemoglobin'},'Location','Northeast','fontsize',12)
xlabel('Wavelength(nm)','fontsize',20,'FontWeight','bold');
ylabel('Nonlinear susceptibility (10^-^2^2 m^2V^-^2)','fontsize',20,'FontWeight','bold');
box on;
set(gca, 'YAxisLocation', 'right')
set(gca,'XTick',1236:6:1284,'FontWeight','bold') 
ylim([0.8 2]);
set(gca,'YTick',0.9:0.2:2,'FontWeight','bold')
set(gca,'fontsize',15)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])

figure('Name','Susceptibility Spectra 3') % Linoleic acid; Oleic acid;
hold on
plot(3*ln,Xexp_LA*1E22,'color',ctmap(6,:),'LineWidth',5)
plot(3*ln,Xexp_OA*1E22,'color',ctmap(7,:),'LineWidth',5)
hold off
xlim([1230 1290]); 
legend({'Linoleic acid','Oleic acid'},'Location','Northeast','fontsize',12)
xlabel('Wavelength(nm)','fontsize',20,'FontWeight','bold');
ylabel('Nonlinear susceptibility (10^-^2^2 m^2V^-^2)','fontsize',20,'FontWeight','bold');
box on;
set(gca, 'YAxisLocation', 'right')
set(gca,'XTick',1236:6:1284,'FontWeight','bold')
ylim([1.6 2.8]);
set(gca,'YTick',1.8:0.2:2.6,'FontWeight','bold')
set(gca,'fontsize',15)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])

