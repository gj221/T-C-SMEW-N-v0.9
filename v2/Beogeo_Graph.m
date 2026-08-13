
Lkday = mean(Lk)*24; 
ZBIOG = Zbio/1000;
rsd_mean = mean(rsd);

%%
f1 = figure(1);
subplot(3,1,1)
plot(B(:,1),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gC/m^2')
plot(B(:,2),'g','LineWidth', 1.5);
plot(B(:,3),'r','LineWidth', 1.5);
legend('Ab. Met','Ab Str','Ab Str Lig','Location','best')

subplot(3,1,2)
plot(B(:,6),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gC/m^2')
plot(B(:,7),'g','LineWidth', 1.5);
plot(B(:,8),'r','LineWidth', 1.5);
legend('Be. Met','Be Str','Be Str Lig','Location','best')

subplot(3,1,3)
plot(B(:,4),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gC/m^2')
plot(B(:,5),'g','LineWidth', 1.5);
legend('Ab. Wood','Ab Wood Lign','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

f2 = figure(2);
subplot(2,2,1)
plot(B(:,9),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,10),'k','LineWidth', 1.5);
plot(B(:,11),'y','LineWidth', 1.5);
legend('SOM POC Lign','SOM POC - Cell','SOM MOC','Location','best')
subplot(2,2,2)
plot(B(:,12),'k','LineWidth', 1.5);
hold on; grid on;
plot(B(:,13),'r','LineWidth', 1.5);
legend('DOC-B','DOC-F','Location','best')
subplot(2,2,3)
plot(B(:,18),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,19),'k','LineWidth', 1.5);
plot(B(:,20),'r','LineWidth', 1.5);
plot(B(:,21),'b','LineWidth', 1.5);
plot(B(:,22),'y','LineWidth', 1.5);
title('Carbon Pool')
xlabel('Days'); ylabel('gC/m^2')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Earthworms','Location','best')
subplot(2,2,4)
plot(B(:,14),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,15),'k','LineWidth', 1.5);
plot(B(:,16),'r','LineWidth', 1.5);
plot(B(:,17),'b','LineWidth', 1.5);
title('Carbon Pool')
xlabel('Days'); ylabel('gC/m^2')
legend('EM-POC-B','EM-POC-F','EM-MOC-B','EM-MOC-F','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

f3 = figure(3);
subplot(2,2,1)
plot(B(:,23),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gN/m^2')
plot(B(:,24),'g','LineWidth', 1.5);
plot(B(:,25),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
subplot(2,2,2)
plot(B(:,26),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM','Location','best')
subplot(2,2,3)
plot(B(:,27),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,28),'k','LineWidth', 1.5);
plot(B(:,29),'r','LineWidth', 1.5);
plot(B(:,30),'b','LineWidth', 1.5);
title('Nitrogen Pool')
xlabel('Days'); ylabel('gN/m^2')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Location','best')
subplot(2,2,4)
plot(B(:,31),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,32),'k','LineWidth', 1.5);
plot(B(:,33),'b','LineWidth', 1.5);
title('Nitrogen Pool')
xlabel('Days'); ylabel('gN/m^2')
legend('NH4+ ','NO3-','DON','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

f4 = figure(4);
subplot(3,2,1)
plot(B(:,35),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gP/m^2')
plot(B(:,36),'g','LineWidth', 1.5);
plot(B(:,37),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
subplot(3,2,2)
plot(B(:,38),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM','Location','best')
subplot(3,2,3)
plot(B(:,39),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,40),'k','LineWidth', 1.5);
plot(B(:,41),'r','LineWidth', 1.5);
plot(B(:,42),'b','LineWidth', 1.5);
title('Phosporus Pool')
xlabel('Days'); ylabel('gP/m^2')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Location','best')
subplot(3,2,4)
plot(B(:,43),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,47),'r','LineWidth', 1.5);
title('Phosporus Pool')
xlabel('Days'); ylabel('gP/m^2')
legend('Mineral','DOP','Location','best')
subplot(3,2,5)
plot(B(:,44),'k','LineWidth', 1.5);
hold on; grid on;
title('Phosporus Pool')
xlabel('Days'); ylabel('gP/m^2')
legend('Primary Material','Location','best')
subplot(3,2,6)
plot(B(:,46),'k','LineWidth', 1.5);
hold on; grid on;
plot(B(:,45),'g','LineWidth', 1.5);
xlabel('Days'); ylabel('gP/m^2')
legend('Occluded','Secondary','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);


f5 = figure(5);
subplot(3,2,1)
plot(B(:,48),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gK/m^2')
plot(B(:,49),'g','LineWidth', 1.5);
plot(B(:,50),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
subplot(3,2,2)
plot(B(:,51),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM','Location','best')
subplot(3,2,3)
plot(B(:,52),'g','LineWidth', 1.5);
hold on; grid on;
title('Potassium Pool')
xlabel('Days'); ylabel('gK/m^2')
legend('Mineral Solution ','Location','best')
subplot(3,2,4)
plot(B(:,53),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,54),'k','LineWidth', 1.5);
title('Potassium Pool')
xlabel('Days'); ylabel('gK/m^2')
legend('Excheangeable ','Non-Excheangeable','Location','best')
subplot(3,2,5)
plot(B(:,55),'g','LineWidth', 1.5);
hold on; grid on;
legend('Primary Minerals','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

f15 = figure(15);
subplot(3,2,1)
plot(B(:,56),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gCa/m^2')
plot(B(:,57),'g','LineWidth', 1.5);
plot(B(:,58),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
subplot(3,2,2)
plot(B(:,59),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM','Location','best')
subplot(3,2,3)
plot(B(:,60),'g','LineWidth', 1.5);
hold on; grid on;
title('Calcium Pool')
xlabel('Days'); ylabel('gCa/m^2')
legend('Mineral Solution ','Location','best')
subplot(3,2,4)
plot(B(:,61),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,62),'k','LineWidth', 1.5);
title('Calcium Pool')
xlabel('Days'); ylabel('gCa/m^2')
legend('Excheangeable ','Non-Excheangeable','Location','best')
subplot(3,2,5)
plot(B(:,63),'g','LineWidth', 1.5);
hold on; grid on;
legend('Primary Minerals','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

f16 = figure(16);
subplot(3,2,1)
plot(B(:,64),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gMg/m^2')
plot(B(:,65),'g','LineWidth', 1.5);
plot(B(:,66),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
subplot(3,2,2)
plot(B(:,67),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM','Location','best')
subplot(3,2,3)
plot(B(:,68),'g','LineWidth', 1.5);
hold on; grid on;
title('Magnesium Pool')
xlabel('Days'); ylabel('gMg/m^2')
legend('Mineral Solution ','Location','best')
subplot(3,2,4)
plot(B(:,69),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,70),'k','LineWidth', 1.5);
title('Magnesium Pool')
xlabel('Days'); ylabel('gMg/m^2')
legend('Excheangeable ','Non-Excheangeable','Location','best')
subplot(3,2,5)
plot(B(:,71),'g','LineWidth', 1.5);
hold on; grid on;
legend('Primary Minerals','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

f17 = figure(17);
subplot(3,2,1)
plot(B(:,72),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gSi/m^2')
plot(B(:,73),'g','LineWidth', 1.5);
plot(B(:,74),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
subplot(3,2,2)
plot(B(:,75),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM','Location','best')
subplot(3,2,3)
plot(B(:,76),'g','LineWidth', 1.5);
hold on; grid on;
title('Silicate Pool')
xlabel('Days'); ylabel('gSi/m^2')
legend('Mineral Solution ','Location','best')
subplot(3,2,4)
plot(B(:,77),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,78),'k','LineWidth', 1.5);
title('Silicate Pool')
xlabel('Days'); ylabel('gSi/m^2')
legend('Excheangeable ','Non-Excheangeable','Location','best')
subplot(3,2,5)
plot(B(:,79),'g','LineWidth', 1.5);
hold on; grid on;
legend('Primary Minerals','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);


f6 = figure(6);
plot((B(:,2)+B(:,3)+B(:,1))./B(:,23),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,24),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,25),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,26),'b','LineWidth', 1.5);
plot(B(:,18)./B(:,27),'r','LineWidth', 1.5);
plot(B(:,19)./B(:,28),'y','LineWidth', 1.5);
plot(B(:,20)./B(:,29),'c','LineWidth', 1.5);
plot(B(:,21)./B(:,30),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
plot(B(:,22)./B(:,34),'Color',[0.06 0.7 0.6],'LineWidth', 1.5);
title('C:N Ratio')
xlabel('Days'); ylabel('C:N')
legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Earthworms','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
%%%%%%%%%%%%%%%%%%%%%%%%%


f7 = figure(7);
plot((B(:,2)+B(:,3)+B(:,1))./B(:,35),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,36),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,37),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,38),'b','LineWidth', 1.5);
plot(B(:,18)./B(:,39),'r','LineWidth', 1.5);
plot(B(:,19)./B(:,40),'y','LineWidth', 1.5);
plot(B(:,20)./B(:,41),'c','LineWidth', 1.5);
plot(B(:,21)./B(:,42),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
title('C:P Ratio')
xlabel('Days'); ylabel('C:P')
legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
%%%%%%%%%%%%%%%%%%%%%%%%%


f8 = figure(8);
plot((B(:,2)+B(:,3)+B(:,1))./B(:,48),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,49),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,50),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,51),'b','LineWidth', 1.5);
title('C:K Ratio')
xlabel('Days'); ylabel('C:K')
legend('AG Litter','AG Wood','BG Litter','SOM','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
%%%%%%%%%%%%%%%%%%%%%%%%%

f18 = figure(18);
plot((B(:,2)+B(:,3)+B(:,1))./B(:,56),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,57),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,58),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,59),'b','LineWidth', 1.5);
title('C:Ca Ratio')
xlabel('Days'); ylabel('C:Ca')
legend('AG Litter','AG Wood','BG Litter','SOM','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
%%%%%%%%%%%%%%%%%%%%%%%%%


f19 = figure(19);
plot((B(:,2)+B(:,3)+B(:,1))./B(:,64),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,65),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,66),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,67),'b','LineWidth', 1.5);
title('C:Mg Ratio')
xlabel('Days'); ylabel('C:Mg')
legend('AG Litter','AG Wood','BG Litter','SOM','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
%%%%%%%%%%%%%%%%%%%%%%%%%


f20 = figure(20);
plot((B(:,2)+B(:,3)+B(:,1))./B(:,72),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,73),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,74),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,75),'b','LineWidth', 1.5);
title('C:Si Ratio')
xlabel('Days'); ylabel('C:Si')
legend('AG Litter','AG Wood','BG Litter','SOM','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
%%%%%%%%%%%%%%%%%%%%%%%%%

% f9 = figure(9);
% subplot(2,2,1)
% plot(R_litter,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(R_litter-R_litter_sur,'m','LineWidth', 1.5);
% plot(R_microbe,'b','LineWidth', 1.5);
% title('Respiration Het.')
% plot(R_ew,'g','LineWidth', 1.5);
% xlabel('Days'); ylabel('[gC/m2 day]')
% legend('Litter','Litter below','Microbe','Earthworms','Location','best')
% subplot(2,2,2)
% plot(VOL,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(N2flx,'b','LineWidth', 1.5);
% title('N - Fluxes')
% xlabel('Days'); ylabel('[gN/m2 day]')
% legend('NH_4 Vol.','N_2','Location','best')
% subplot(2,2,3)
% plot(Min_N,'k','LineWidth', 1.5);
% hold on ;  grid on 
% title('N - Fluxes')
% xlabel('Days'); ylabel('[gN/m2 day]')
% legend('Min-N','Location','best')
% subplot(2,2,4)
% plot(Min_P,'k','LineWidth', 1.5);
% hold on ;  grid on 
% title('P - Fluxes')
% xlabel('Days'); ylabel('[gP/m2 day]')
% legend('Min-P','Location','best')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% f10 = figure(10);
% subplot(3,2,1)
% plot(NH4uptake_L,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(NO3uptake_L,'b','LineWidth', 1.5);
% title('N Uptake.')
% xlabel('Days'); ylabel('[gN/m2 day]')
% legend('NH_4','NO_3','Location','best')
% subplot(3,2,2)
% plot(Puptake_L,'g','LineWidth', 1.5);
% hold on; grid on;
% title('P Uptake')
% xlabel('Days'); ylabel('[gP/m2 day]')
% subplot(3,2,3)
% plot(Kuptake_L,'m','LineWidth', 1.5);
% hold on; grid on;
% title('K Uptake')
% xlabel('Days'); ylabel('[gK/m2 day]')
% subplot(3,2,4)
% plot(Cauptake_L,'m','LineWidth', 1.5);
% hold on; grid on;
% title('Ca Uptake')
% xlabel('Days'); ylabel('[gCa/m2 day]')
% subplot(3,2,5)
% plot(Mguptake_L,'m','LineWidth', 1.5);
% hold on; grid on;
% title('Mg Uptake')
% xlabel('Days'); ylabel('[gMg/m2 day]')
% subplot(3,2,6)
% plot(Siuptake_L,'m','LineWidth', 1.5);
% hold on; grid on;
% title('Si Uptake')
% xlabel('Days'); ylabel('[gSi/m2 day]')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% f21 = figure(21);
% subplot(3,2,1)
% plot(Nreserve_L,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(Nreserve_H,'b','LineWidth', 1.5);
% title('N reserve.')
% legend('L-Veg','H-Veg','Location','best')
% xlabel('Days'); ylabel('[gN/m2 day]')
% subplot(3,2,2)
% plot(Preserve_L,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(Preserve_H,'b','LineWidth', 1.5);
% title('P reserve')
% xlabel('Days'); ylabel('[gP/m2 day]')
% subplot(3,2,3)
% plot(Kreserve_L,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(Kreserve_H,'b','LineWidth', 1.5);
% title('K reserve')
% xlabel('Days'); ylabel('[gK/m2 day]')
% subplot(3,2,4)
% plot(Careserve_L,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(Careserve_H,'b','LineWidth', 1.5);
% title('Ca reserve')
% xlabel('Days'); ylabel('[gCa/m2 day]')
% subplot(3,2,5)
% plot(Mgreserve_L,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(Mgreserve_H,'b','LineWidth', 1.5);
% title('Mg reserve')
% xlabel('Days'); ylabel('[gMg/m2 day]')
% subplot(3,2,6)
% plot(Sireserve_L,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(Sireserve_H,'b','LineWidth', 1.5);
% title('Si reserve')
% xlabel('Days'); ylabel('[gSi/m2 day]')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% f11 = figure(11);
% subplot(2,2,1)
% plot(LEAK_DOC,'b','LineWidth', 1.5);
% hold on; grid on;
% title('DOC Leaching')
% xlabel('Days'); ylabel('[gC/m2 day]')
% subplot(2,2,2)
% plot(LEAK_NH4,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(LEAK_NO3,'b','LineWidth', 1.5);
% hold on; grid on;
% plot(LEAK_DON,'m','LineWidth', 1.5);
% title('N Leaching')
% xlabel('Days'); ylabel('[gN/m2 day]')
% legend('NH_4','NO_3','DON','Location','best')
% subplot(2,2,3)
% plot(LEAK_P,'g','LineWidth', 1.5);
% hold on; grid on;
% plot(LEAK_DOP,'m','LineWidth', 1.5);
% title('P Leaching')
% legend('P','PON','Location','best')
% xlabel('Days'); ylabel('[gP/m2 day]')
% subplot(2,2,4)
% plot(LEAK_K,'m','LineWidth', 1.5);
% hold on; grid on;
% title('K Leaching')
% xlabel('Days'); ylabel('[gK/m2 day]')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
% %%%%%%%%%%%%%%%%%%%%%%%%%
% 
% f12 = figure(12);
% subplot(2,2,1)
% plot(LEAK_DOC./(Lkday)*1000,'b','LineWidth', 1.5);
% hold on; grid on;
% title('DOC Conc.')
% xlabel('Days'); ylabel('[mg/l]')
% subplot(2,2,2)
% plot(LEAK_NH4./(Lkday)*1000,'r','LineWidth', 1.5);
% hold on; grid on;
% plot(LEAK_NO3./(Lkday)*1000,'b','LineWidth', 1.5);
% hold on; grid on;
% plot(LEAK_DON./(Lkday)*1000,'m','LineWidth', 1.5);
% title('N Conc.')
% xlabel('Days'); ylabel('[mg/l]')
% legend('NH_4','NO_3','DON','Location','best')
% subplot(2,2,3)
% plot(LEAK_P./(Lkday)*1000*1000,'g','LineWidth', 1.5);
% hold on; grid on;
% plot(LEAK_DOP./(Lkday)*1000*1000,'m','LineWidth', 1.5);
% title('P Conc.')
% legend('P','DOP','Location','best')
% xlabel('Days');  ylabel('[ug/ l]')
% subplot(2,2,4)
% plot(LEAK_K./(Lkday)*1000,'m','LineWidth', 1.5);
% hold on; grid on;
% title('K Conc.')
% xlabel('Days'); ylabel('[mg/ l]')
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);


f13 = figure(13);
subplot(3,2,1)
plot((B(:,18)+B(:,19)+B(:,20) + B(:,21))./(B(:,9)+B(:,10)+B(:,11)+B(:,12)+B(:,13)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Microbial/Substrate','Location','best')
subplot(3,2,2)
plot((B(:,22))./(B(:,18)+B(:,19)+B(:,20) + B(:,21)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Earthworms/Microbial','Location','best')
subplot(3,2,3)
plot(sum(B(:,6:21),2)./(ZBIOG*rsd_mean),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('SOC [gC /kg soil]')
subplot(3,2,4)
plot(sum(B(:,6:21),2),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('SOC [gC /m^2]')
subplot(3,2,5)
plot(((B(:,19) + B(:,20) + B(:,21))./B(:,18)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Fungi/Bacteria','Location','best')
subplot(3,2,6)
plot( B(:,19)./(B(:,20)+B(:,21)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Saprotrophic/Mycorrhiza','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);


Rlitter_sub= R_litter-R_litter_sur; 
rb = R_bacteria./R_microbe; 
rf = 1 -rb; 
TSR = Rlitter_sub + R_microbe + R_ew; 
f14 = figure(14);
subplot(1,1,1)
plot((R_bacteria + rb.*Rlitter_sub)./TSR,'r','LineWidth', 1.5);
hold on; grid on;
plot((rf.*R_microbe + rf.*Rlitter_sub)./TSR,'b','LineWidth', 1.5);
title('Respiration Het.')
plot(R_ew./TSR,'g','LineWidth', 1.5);
xlabel('Days'); ylabel('[%]')
legend('Bacteria','Fungi','Earthworms','Location','best')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);



%%
if exist("B_out",'var')
    
    figure(1);
    subplot(3,1,1)
    plot(length(B),B_out(1),'bo');
    hold on; grid on;  ylabel('gC/m^2')
    plot(length(B),B_out(2),'go');
    plot(length(B),B_out(3),'ro');
    legend('Ab. Met','Ab Str','Ab Str Lig','Location','best')
    
    subplot(3,1,2)
    plot(length(B),B_out(6),'bo');
    hold on; grid on;  ylabel('gC/m^2')
    plot(length(B),B_out(7),'go');
    plot(length(B),B_out(8),'ro');
    legend('Be. Met','Be Str','Be Str Lig','Location','best')
    
    subplot(3,1,3)
    plot(length(B),B_out(4),'bo');
    hold on; grid on;  ylabel('gC/m^2')
    plot(length(B),B_out(5),'go');
    legend('Ab. Wood','Ab Wood Lign','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
    
    figure(2);
    subplot(2,2,1)
    plot(length(B),B_out(9),'go');
    hold on; grid on;
    plot(length(B),B_out(10),'ko');
    plot(length(B),B_out(11),'yo');
    legend('SOM POC Lign','SOM POC - Cell','SOM MOC','Location','best')
    subplot(2,2,2)
    plot(length(B),B_out(12),'ko');
    hold on; grid on;
    plot(length(B),B_out(13),'ro');
    legend('DOC-B','DOC-F','Location','best')
    subplot(2,2,3)
    plot(length(B),B_out(18),'go');
    hold on; grid on;
    plot(length(B),B_out(19),'ko');
    plot(length(B),B_out(20),'ro');
    plot(length(B),B_out(21),'bo');
    plot(length(B),B_out(22),'yo');
    title('Carbon Pool')
    xlabel('Days'); ylabel('gC/m^2')
    legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Earthworms','Location','best')
    subplot(2,2,4)
    plot(length(B),B_out(14),'go');
    hold on; grid on;
    plot(length(B),B_out(15),'ko');
    plot(length(B),B_out(16),'ro');
    plot(length(B),B_out(17),'bo');
    title('Carbon Pool')
    xlabel('Days'); ylabel('gC/m^2')
    legend('EM-POC-B','EM-POC-F','EM-MOC-B','EM-MOC-F','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
    
    figure(3);
    subplot(2,2,1)
    plot(length(B),B_out(23),'bo');
    hold on; grid on;  ylabel('gN/m^2')
    plot(length(B),B_out(24),'go');
    plot(length(B),B_out(25),'ro');
    legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
    subplot(2,2,2)
    plot(length(B),B_out(26),'go');
    hold on; grid on;
    legend('SOM','Location','best')
    subplot(2,2,3)
    plot(length(B),B_out(27),'go');
    hold on; grid on;
    plot(length(B),B_out(28),'ko');
    plot(length(B),B_out(29),'ro');
    plot(length(B),B_out(30),'bo');
    title('Nitrogen Pool')
    xlabel('Days'); ylabel('gN/m^2')
    legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Location','best')
    subplot(2,2,4)
    plot(length(B),B_out(31),'go');
    hold on; grid on;
    plot(length(B),B_out(32),'ko');
    plot(length(B),B_out(33),'bo');
    title('Nitrogen Pool')
    xlabel('Days'); ylabel('gN/m^2')
    legend('NH4+ ','NO3-','DON','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
    
    figure(4);
    subplot(3,2,1)
    plot(length(B),B_out(35),'bo');
    hold on; grid on;  ylabel('gP/m^2')
    plot(length(B),B_out(36),'go');
    plot(length(B),B_out(37),'ro');
    legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
    subplot(3,2,2)
    plot(length(B),B_out(38),'go');
    hold on; grid on;
    legend('SOM','Location','best')
    subplot(3,2,3)
    plot(length(B),B_out(39),'go');
    hold on; grid on;
    plot(length(B),B_out(40),'ko');
    plot(length(B),B_out(41),'ro');
    plot(length(B),B_out(42),'bo');
    title('Phosporus Pool')
    xlabel('Days'); ylabel('gP/m^2')
    legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Location','best')
    subplot(3,2,4)
    plot(length(B),B_out(43),'go');
    hold on; grid on;
    plot(length(B),B_out(47),'ro');
    title('Phosporus Pool')
    xlabel('Days'); ylabel('gP/m^2')
    legend('Mineral','DOP','Location','best')
    subplot(3,2,5)
    plot(length(B),B_out(44),'ko');
    hold on; grid on;
    title('Phosporus Pool')
    xlabel('Days'); ylabel('gP/m^2')
    legend('Primary Material','Location','best')
    subplot(3,2,6)
    plot(length(B),B_out(46),'ko');
    hold on; grid on;
    plot(length(B),B_out(45),'go');
    xlabel('Days'); ylabel('gP/m^2')
    legend('Occluded','Secondary','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
    
    
    figure(5);
    subplot(3,2,1)
    plot(length(B),B_out(48),'bo');
    hold on; grid on;  ylabel('gK/m^2')
    plot(length(B),B_out(49),'go');
    plot(length(B),B_out(50),'ro');
    legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
    subplot(3,2,2)
    plot(length(B),B_out(51),'go');
    hold on; grid on;
    legend('SOM','Location','best')
    subplot(3,2,3)
    plot(length(B),B_out(52),'go');
    hold on; grid on;
    title('Potassium Pool')
    xlabel('Days'); ylabel('gK/m^2')
    legend('Mineral Solution ','Location','best')
    subplot(3,2,4)
    plot(length(B),B_out(53),'go');
    hold on; grid on;
    plot(length(B),B_out(54),'ko');
    title('Potassium Pool')
    xlabel('Days'); ylabel('gK/m^2')
    legend('Excheangeable ','Non-Excheangeable','Location','best')
    subplot(3,2,5)
    plot(length(B),B_out(55),'go');
    hold on; grid on;
    legend('Primary Minerals','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
    
    figure(15);
    subplot(3,2,1)
    plot(length(B),B_out(56),'bo');
    hold on; grid on;  ylabel('gCa/m^2')
    plot(length(B),B_out(57),'go');
    plot(length(B),B_out(58),'ro');
    legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
    subplot(3,2,2)
    plot(length(B),B_out(59),'go');
    hold on; grid on;
    legend('SOM','Location','best')
    subplot(3,2,3)
    plot(length(B),B_out(60),'go');
    hold on; grid on;
    title('Calcium Pool')
    xlabel('Days'); ylabel('gCa/m^2')
    legend('Mineral Solution ','Location','best')
    subplot(3,2,4)
    plot(length(B),B_out(61),'go');
    hold on; grid on;
    plot(length(B),B_out(62),'ko');
    title('Calcium Pool')
    xlabel('Days'); ylabel('gCa/m^2')
    legend('Excheangeable ','Non-Excheangeable','Location','best')
    subplot(3,2,5)
    plot(length(B),B_out(63),'go');
    hold on; grid on;
    legend('Primary Minerals','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
    
    figure(16);
    subplot(3,2,1)
    plot(length(B),B_out(64),'bo');
    hold on; grid on;  ylabel('gMg/m^2')
    plot(length(B),B_out(65),'go');
    plot(length(B),B_out(66),'ro');
    legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
    subplot(3,2,2)
    plot(length(B),B_out(67),'go');
    hold on; grid on;
    legend('SOM','Location','best')
    subplot(3,2,3)
    plot(length(B),B_out(68),'go');
    hold on; grid on;
    title('Magnesium Pool')
    xlabel('Days'); ylabel('gMg/m^2')
    legend('Mineral Solution ','Location','best')
    subplot(3,2,4)
    plot(length(B),B_out(69),'go');
    hold on; grid on;
    plot(length(B),B_out(70),'ko');
    title('Magnesium Pool')
    xlabel('Days'); ylabel('gMg/m^2')
    legend('Excheangeable ','Non-Excheangeable','Location','best')
    subplot(3,2,5)
    plot(length(B),B_out(71),'go');
    hold on; grid on;
    legend('Primary Minerals','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);
    
    figure(17);
    subplot(3,2,1)
    plot(length(B),B_out(72),'bo');
    hold on; grid on;  ylabel('gSi/m^2')
    plot(length(B),B_out(73),'go');
    plot(length(B),B_out(74),'ro');
    legend('Ab. Lit','Ab Wod','Be Lit','Location','best')
    subplot(3,2,2)
    plot(length(B),B_out(75),'go');
    hold on; grid on;
    legend('SOM','Location','best')
    subplot(3,2,3)
    plot(length(B),B_out(76),'go');
    hold on; grid on;
    title('Silicate Pool')
    xlabel('Days'); ylabel('gSi/m^2')
    legend('Mineral Solution ','Location','best')
    subplot(3,2,4)
    plot(length(B),B_out(77),'go');
    hold on; grid on;
    plot(length(B),B_out(78),'ko');
    title('Silicate Pool')
    xlabel('Days'); ylabel('gSi/m^2')
    legend('Excheangeable ','Non-Excheangeable','Location','best')
    subplot(3,2,5)
    plot(length(B),B_out(79),'go');
    hold on; grid on;
    legend('Primary Minerals','Location','best')
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.01 0.01 0.98 0.98]);

    
end