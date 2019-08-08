%%  Interface1 with avg material

% Relative THG spectra

figure('Name','Relative THG spectra')
hold on;
for q6=1:isp/2
    if (spec(2*q6)~=0)
        plot(ln,avg_n(q6,:),'Color',ctmap(q6,:),'LineWidth',5);
    end
end
hold off;
xlim([400 444]);
set(gca,'XTick',405:5:440,'FontWeight','bold') 
set(gca,'YTick',0:0.5:1.5,'FontWeight','bold') 
xlabel('Third-Harmonic Wavelength(nm)','fontsize',15,'FontWeight','bold');
ylabel('Relative intensity (a.u.)','fontsize',15,'FontWeight','bold');
legend(Tag2,'Fontsize',12);
box on;
set(gca,'fontsize',17)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])

% Normalized THG spectra

figure('Name','Normalized THG spectra 1') % Fused silica; Water; Melanin
hold on;
for q6=[1 2 5]
    if (spec(2*q6)~=0)
        plot(ln,avg_nn(q6,:),'Color',ctmap(q6,:),'LineWidth',5);
    end
end
hold off;
xlim([410 430]);
ylim([0 1.2]);
xlabel('Wavelength(nm)','fontsize',20,'FontWeight','bold');
ylabel('Normalized Intensity','fontsize',20,'FontWeight','bold');
set(gca,'XTick',412:2:428,'FontWeight','bold') 
set(gca,'YTick',0.2:0.2:1,'FontWeight','bold') 
legend(Tag2(1,[1 2 5]),'Fontsize',12);
box on;
set(gca,'fontsize',15)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])

figure('Name','Normalized THG spectra 2') % De-oxy; Oxy;
hold on;
for q6=[3 4]
    if (spec(2*q6)~=0)
        plot(ln,avg_nn(q6,:),'Color',ctmap(q6,:),'LineWidth',5);
    end
end
hold off;
xlim([410 430]);
ylim([0 1.2]);
xlabel('Wavelength(nm)','fontsize',20,'FontWeight','bold');
ylabel('Normalized Intensity','fontsize',20,'FontWeight','bold');
set(gca,'XTick',412:2:428,'FontWeight','bold') 
set(gca,'YTick',0.2:0.2:1,'FontWeight','bold') 
legend(Tag2(1,3:4),'Fontsize',12);
box on;
set(gca,'fontsize',15)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])

figure('Name','Normalized THG spectra 3') % Linoleic acid; Oleic acid;
hold on;
for q6=[6 7]
    if (spec(2*q6)~=0)
        plot(ln,avg_nn(q6,:),'Color',ctmap(q6,:),'LineWidth',5);
    end
end
hold off;
xlim([410 430]);
ylim([0 1.2]);
xlabel('Wavelength(nm)','fontsize',20,'FontWeight','bold');
ylabel('Normalized Intensity','fontsize',20,'FontWeight','bold');
set(gca,'XTick',412:2:428,'FontWeight','bold') 
set(gca,'YTick',0.2:0.2:1,'FontWeight','bold') 
legend(Tag2(1,6:7),'Fontsize',12);
box on;
set(gca,'fontsize',15)
set(gca,'linewidth',2)
set(gcf,'color','w');
set(gcf, 'Position', [200 100 805 550]);
pbaspect([1.7 1 1])
