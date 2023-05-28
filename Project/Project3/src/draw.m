clear;
p=importdata('FehlbergRK.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/FehlbergRK.eps');
p=importdata('ABF.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/ABF.eps');
clear;
p=importdata('ADM.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/ADM.eps');
clear;
p=importdata('BDF.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/BDF.eps');
clear;
p=importdata('DormandPrinceRK.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/DormandPrinceRK.eps');
clear;
p=importdata('ESDIRK.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/ESDIRK.eps');
clear;
p=importdata('FehlbergRK.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/FehlbergRK.eps');
clear;
p=importdata('GaussLegendreRK.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/GaussLegendreRK.eps');
clear;
p=importdata('RK.data');
plot(p(:,1), p(:,2))
saveas(gcf,'../pic/RK.eps');

clear;
base_path = '../data';
all_file_path = fullfile(base_path);

file = dir(all_file_path);
for i=3:length(file)
    figure(1);
    clf
    p=importdata(['../data/',file(i).name]);
    
    subplot(1,2,1);
    semilogy(p(:,3),p(:,1));
    set(gca,'XDir','reverse');
    xlabel('Time Step');
    ylabel('Max-norm Error');
    subplot(1,2,2);
    plot(p(:,3),p(:,2));
    set(gca,'XDir','reverse');
    xlabel('Time Step');
    ylabel('CPU Time');
    figure(1);
    subplot(1,2,1);
    hold on;
    subplot(1,2,2);
    hold on;
    saveas(gcf,['../pic/',file(i).name,'.eps']);
    file(i).name
    (log2(p(1,1)/p(2,1))+log2(p(2,1)/p(3,1)))/2
end