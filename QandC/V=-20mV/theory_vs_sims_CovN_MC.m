
% coded by Roberto Fernández Galán, rfgalan@case.edu, December 2011

clear
close all
clc

%%

load theoretical_CovN_full_MC.mat

QNatf = QNa;
QKtf = QK;

C0Natf = C0Na;
C0Ktf = C0K;

load simulations_CovN_full_MC.mat

QNasf = QNa;
QKsf = QK;

C0Nasf = C0Na;
C0Ksf = C0K;

load theoretical_CovN_simp_MC.mat

QNats = QNa;
QKts = QK;

C0Nats = C0Na;
C0Kts = C0K;

load simulations_CovN_simp_MC.mat

QNass = QNa;
QKss = QK;

C0Nass = C0Na;
C0Kss = C0K;

%%

climQNa = [min([QNatf(:); QNasf(:); QNats(:); QNass(:)])  max([QNatf(:); QNasf(:); QNats(:); QNass(:)])];

climCNa = [min([C0Natf(:); C0Nasf(:); C0Nats(:); C0Nass(:)])  max([C0Natf(:); C0Nasf(:); C0Nats(:); C0Nass(:)])];

climQK = [min([QKtf(:); QKsf(:); QKts(:); QKss(:)])  max([QKtf(:); QKsf(:); QKts(:); QKss(:)])];

climCK = [min([C0Ktf(:); C0Ksf(:); C0Kts(:); C0Kss(:)])  max([C0Ktf(:); C0Ksf(:); C0Kts(:); C0Kss(:)])];

%%

figure(1)

colormap(gray)

subplot(4,4,1)
imagesc(QNatf,climQNa)
title('theoretical Q_{Na} for full Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
subplot(4,4,2)
imagesc(QNasf,climQNa)
title('simulated Q_{Na} for full Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
subplot(4,4,3)
imagesc(QNats,climQNa)
title('theoretical Q_{Na} for simplified Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
subplot(4,4,4)
imagesc(QNass,climQNa)
title('simulated Q_{Na} for simplified Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(4,4,5)
imagesc(QKtf,climQK)
title('theoretical Q_{K} for full Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
subplot(4,4,6)
imagesc(QKsf,climQK)
title('simulated Q_{K} for full Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
subplot(4,4,7)
imagesc(QKts,climQK)
title('theoretical Q_{K} for simplified Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
subplot(4,4,8)
imagesc(QKss,climQK)
title('simulated Q_{K} for simplified Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
colorbar

subplot(4,4,9)
imagesc(C0Natf,climCNa)
title('theoretical C_{Na} for full Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
subplot(4,4,10)
imagesc(C0Nasf,climCNa)
title('simulated C_{Na} for full Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
subplot(4,4,11)
imagesc(C0Nats,climCNa)
title('theoretical C_{Na} for simplified Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
subplot(4,4,12)
imagesc(C0Nass,climCNa)
title('simulated C_{Na} for simplified Markov chain')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(4,4,13)
imagesc(C0Ktf,climQK)
title('theoretical C_{K} for full Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
subplot(4,4,14)
imagesc(C0Ksf,climQK)
title('simulated C_{K} for full Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
subplot(4,4,15)
imagesc(C0Kts,climQK)
title('theoretical C_{K} for simplified Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
subplot(4,4,16)
imagesc(C0Kss,climQK)
title('simulated C_{K} for simplified Markov chain')
set(gca,'xtick',1:5,'xticklabel',1:5,'ytick',1:5,'yticklabel',1:5)
axis equal
axis tight
colorbar

%%

temp1 = QNatf(7:8,7:8);
temp2 = QNasf(7:8,7:8);
temp3 = QNats(7:8,7:8);
temp4 = QNass(7:8,7:8);
climQNa = [min([temp1(:); temp2(:); temp3(:); temp4(:)])  max([temp1(:); temp2(:); temp3(:); temp4(:)])];

temp1 = C0Natf(7:8,7:8);
temp2 = C0Nasf(7:8,7:8);
temp3 = C0Nats(7:8,7:8);
temp4 = C0Nass(7:8,7:8);
climCNa = [min([temp1(:); temp2(:); temp3(:); temp4(:)])  max([temp1(:); temp2(:); temp3(:); temp4(:)])];

temp1 = QKtf(4:5,4:5);
temp2 = QKsf(4:5,4:5);
temp3 = QKts(4:5,4:5);
temp4 = QKss(4:5,4:5);
climQK = [min([temp1(:); temp2(:); temp3(:); temp4(:)])  max([temp1(:); temp2(:); temp3(:); temp4(:)])];

temp1 = C0Ktf(4:5,4:5);
temp2 = C0Ksf(4:5,4:5);
temp3 = C0Kts(4:5,4:5);
temp4 = C0Kss(4:5,4:5);
climCK = [min([temp1(:); temp2(:); temp3(:); temp4(:)])  max([temp1(:); temp2(:); temp3(:); temp4(:)])];

%%

figure(2)

colormap(gray)

subplot(4,4,1)
imagesc(QNatf(7:8,7:8),climQNa)
title('theoretical Q_{Na} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
subplot(4,4,2)
imagesc(QNasf(7:8,7:8),climQNa)
title('simulated Q_{Na} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
subplot(4,4,3)
imagesc(QNats(7:8,7:8),climQNa)
title('theoretical Q_{Na} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
subplot(4,4,4)
imagesc(QNass(7:8,7:8),climQNa)
title('simulated Q_{Na} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
colorbar

subplot(4,4,5)
imagesc(QKtf(4:5,4:5),climQK)
title('theoretical Q_{K} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
subplot(4,4,6)
imagesc(QKsf(4:5,4:5),climQK)
title('simulated Q_{K} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
subplot(4,4,7)
imagesc(QKts(4:5,4:5),climQK)
title('theoretical Q_{K} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
subplot(4,4,8)
imagesc(QKss(4:5,4:5),climQK)
title('simulated Q_{K} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
colorbar

subplot(4,4,9)
imagesc(C0Natf(7:8,7:8),climCNa)
title('theoretical C_{Na} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
subplot(4,4,10)
imagesc(C0Nasf(7:8,7:8),climCNa)
title('simulated C_{Na} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
subplot(4,4,11)
imagesc(C0Nats(7:8,7:8),climCNa)
title('theoretical C_{Na} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
subplot(4,4,12)
imagesc(C0Nass(7:8,7:8),climCNa)
title('simulated C_{Na} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',7:8,'ytick',1:2,'yticklabel',7:8)
axis equal
axis tight
colorbar

subplot(4,4,13)
imagesc(C0Ktf(4:5,4:5),climQK)
title('theoretical C_{K} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
subplot(4,4,14)
imagesc(C0Ksf(4:5,4:5),climQK)
title('simulated C_{K} for full Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
subplot(4,4,15)
imagesc(C0Kts(4:5,4:5),climQK)
title('theoretical C_{K} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
subplot(4,4,16)
imagesc(C0Kss(4:5,4:5),climQK)
title('simulated C_{K} for simplified Markov chain')
set(gca,'xtick',1:2,'xticklabel',4:5,'ytick',1:2,'yticklabel',4:5)
axis equal
axis tight
colorbar

%%

figure(3)

subplot(2,4,1)
bar(sqrt([QNatf(8:8,8:8) QNats(8:8,8:8)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{Q_{Na}(8,8)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
ylabel('standard deviation')
xlim([0.25 2.75])
ylim([0 3])

subplot(2,4,2)
bar(sqrt([QNasf(8:8,8:8) QNass(8:8,8:8)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{Q_{Na}(8,8)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
xlim([0.25 2.75])
ylim([0 3])

subplot(2,4,3)
bar(sqrt([C0Natf(8:8,8:8) C0Nats(8:8,8:8)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{C_{Na}(8,8)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
xlim([0.25 2.75])
ylim([0 15])

subplot(2,4,4)
bar(sqrt([C0Nasf(8:8,8:8) C0Nass(8:8,8:8)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{C_{Na}(8,8)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
xlim([0.25 2.75])
ylim([0 15])

subplot(2,4,5)
bar(sqrt([QKtf(5:5,5:5) QKts(5:5,5:5)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{Q_{K}(5,5)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
ylabel('standard deviation')
xlim([0.25 2.75])
ylim([0 4])

subplot(2,4,6)
bar(sqrt([QKsf(5:5,5:5) QKss(5:5,5:5)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{Q_{K}(5,5)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
xlim([0.25 2.75])
ylim([0 4])

subplot(2,4,7)
bar(sqrt([C0Ktf(5:5,5:5) C0Kts(5:5,5:5)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{C_{K}(5,5)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
xlim([0.25 2.75])
ylim([0 40])

subplot(2,4,8)
bar(sqrt([C0Ksf(5:5,5:5) C0Kss(5:5,5:5)]),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{C_{K}(5,5)}$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1:2,'xticklabel',{'full'; 'simplified'})
xlim([0.25 2.75])
ylim([0 40])

%%

figure(4)

subplot(2,4,1)
bar(100*(sqrt(QNats(8:8,8:8))/sqrt(QNatf(8:8,8:8))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{\frac{Q_{Na}^S(8,8)}{Q_{Na}^F(8,8)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
ylabel('relative error (%)')
xlim([0.25 1.75])
ylim([-1 1])

subplot(2,4,2)
bar(100*(sqrt(QNass(8:8,8:8))/sqrt(QNasf(8:8,8:8))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{\frac{Q_{Na}^S(8,8)}{Q_{Na}^F(8,8)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
xlim([0.25 1.75])
ylim([-1 1])

subplot(2,4,3)
bar(100*(sqrt(C0Nats(8:8,8:8))/sqrt(C0Natf(8:8,8:8))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{\frac{C_{Na}^S(8,8)}{C_{Na}^F(8,8)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
xlim([0.25 1.75])
ylim([-8 0])

subplot(2,4,4)
bar(100*(sqrt(C0Nass(8:8,8:8))/sqrt(C0Nasf(8:8,8:8))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{\frac{C_{Na}^S(8,8)}{C_{Na}^F(8,8)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
xlim([0.25 1.75])
ylim([-8 0])

subplot(2,4,5)
bar(100*(sqrt(QKts(5:5,5:5))/sqrt(QKtf(5:5,5:5))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{\frac{Q_{K}^S(5,5)}{Q_{K}^F(5,5)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
ylabel('relative error (%)')
xlim([0.25 1.75])
ylim([-1 1])

subplot(2,4,6)
bar(100*(sqrt(QKss(5:5,5:5))/sqrt(QKsf(5:5,5:5))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{\frac{Q_{K}^S(5,5)}{Q_{K}^F(5,5)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
xlim([0.25 1.75])
ylim([-1 1])

subplot(2,4,7)
bar(100*(sqrt(C0Kts(5:5,5:5))/sqrt(C0Ktf(5:5,5:5))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('theoretical $\sqrt{\frac{C_{K}^S(5,5)}{C_{K}^F(5,5)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
xlim([0.25 1.75])
ylim([-6 0])

subplot(2,4,8)
bar(100*(sqrt(C0Kss(5:5,5:5))/sqrt(C0Ksf(5:5,5:5))-1),'facecolor',[1 1 1]*0.5)
set(gca,'fontsize',16)
title('simulated $\sqrt{\frac{C_{K}^S(5,5)}{C_{K}^F(5,5)}}-1$','Interpreter','latex','fontsize',18)
set(gca,'xtick',1,'xticklabel',[])
xlim([0.25 1.75])
ylim([-6 0])

%%

figure(5)

colormap(gray)

subplot(2,4,1)
imagesc(abs((QNatf-QNats)./QNatf),[0 1])
title('|Q_{Na}^{full} - Q_{Na}^{simp}| / |Q_{Na}^{full}|  (theory)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(2,4,2)
imagesc(abs((QNasf-QNass)./QNasf),[0 1])
title('|Q_{Na}^{full} - Q_{Na}^{simp}| / |Q_{Na}^{full}|  (simulations)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(2,4,3)
imagesc(abs((C0Natf-C0Nats)./C0Natf),[0 1])
title('|C_{Na}^{full} - C_{Na}^{simp}| / |C_{Na}^{full}|  (theory)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(2,4,4)
imagesc(abs((C0Nasf-C0Nass)./C0Nasf),[0 1])
title('|C_{Na}^{full} - C_{Na}^{simp}| / |C_{Na}^{full}|  (simulations)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(2,4,5)
imagesc(abs((QKtf-QKts)./QKtf),[0 1])
title('|Q_{K}^{full} - Q_{K}^{simp}| / |Q_{K}^{full}|  (theory)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(2,4,6)
imagesc(abs((QKsf-QKss)./QKsf),[0 1])
title('|Q_{K}^{full} - Q_{K}^{simp}| / |Q_{K}^{full}|  (simulations)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(2,4,7)
imagesc(abs((C0Ktf-C0Kts)./C0Ktf),[0 1])
title('|C_{K}^{full} - C_{K}^{simp}| / |C_{K}^{full}|  (theory)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar

subplot(2,4,8)
imagesc(abs((C0Ksf-C0Kss)./C0Ksf),[0 1])
title('|C_{K}^{full} - C_{K}^{simp}| / |C_{K}^{full}|  (simulations)')
set(gca,'xtick',1:8,'xticklabel',1:8,'ytick',1:8,'yticklabel',1:8)
axis equal
axis tight
colorbar
