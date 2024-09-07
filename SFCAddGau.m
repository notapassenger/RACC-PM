%% load HB,HD
%% Hua wei data
clear 
addpath(fullfile(pwd,'vast'));
ATF = importdata("ATFRandnAddGau.mat");
% fSpaceChoose = [2:1:40, 41:2:65, 67:16:size(ATF.irMeasured.HB, 3)];
% fSpaceChoose = [2:2:198, 199:16:size(ATF.irMeasured.HB, 3)];
fSpaceChoose = 2:4:201*4;
HBAll = ATF.irMeasured.HB;
HDAll = ATF.irMeasured.HD;
% % HbdesiredPre = ATF.irDesired.HB;
% HbdesiredPre = ATF.irTrue.HB(:, 8, :);
M = size(HBAll, 1);
L = size(HDAll, 2);
HbdesiredPre = ATF.irDesired.HB(:, fSpaceChoose);
% HbdesiredPre = squeeze(ATF.irTrue.HB(:, 8, fSpaceChoose));
% clear ATF;
%%
frePoint = 1:length(fSpaceChoose);
PerformanceChoose = 1:size(HBAll, 4);
nDFT = 3200;
fs = 16000;
f = fs*(fSpaceChoose)/nDFT; 
para.frePoint = frePoint;
para.fs = fs;
para.nDFT = nDFT;
para.freq = f;
para.L = L;
para.M = M;
% load("alpha.mat");
% para.alpha = alpha+5; % Ccon
% para.alpha(69:187) = 30;
load("alphaAddGauFinal.mat");
para.alpha = alpha; % Ccon
% para.alpha(141:150) = alpha(140);
para.mu = 1; % up, AC
% para.rho = 0.05; % rev better
% para.rho = 0.1; % gau better
para.rho = 0.1; % can't equal to 1 % down NSDE
% load('ewNew.mat');
% para.ew = ew+100;
% para.ew(69:187) = 1000;
% load('ewAddGau.mat');
% para.ew = ew;


AC_AP = zeros(length(PerformanceChoose), size(frePoint, 2), length(PerformanceChoose)-1);
NSDE_AP = AC_AP;
AE_AP = AC_AP;
AC_RAP = AC_AP;
NSDE_RAP = AC_AP;
AE_RAP = AC_AP;
AC_VASTNF = AC_AP;
NSDE_VASTNF = AC_AP;
AE_VASTNF = AC_AP;
AC_wcACC = AC_AP;
NSDE_wcACC = AC_AP;
AE_wcACC = AC_AP;
AC_PRACC = AC_AP;
NSDE_PRACC = AC_AP;
AE_PRACC = AC_AP;
AC_ACC = AC_AP;
NSDE_ACC = AC_AP;
AE_ACC = AC_AP;
AC_PM = AC_AP;
NSDE_PM = AC_AP;
AE_PM = AC_AP;

clear ATF;
%%
for k = PerformanceChoose
    rirComputeTh = k;
    rirPerformanceTh = setdiff(PerformanceChoose, rirComputeTh);
    HB = squeeze(HBAll(:, :, fSpaceChoose, rirComputeTh)); 
    HD = squeeze(HDAll(:, :, fSpaceChoose, rirComputeTh));
    HBe = HBAll(:, :, fSpaceChoose, rirPerformanceTh);
    HDe = HDAll(:, :, fSpaceChoose, rirPerformanceTh);
    Hbdesired = HbdesiredPre;
%     for ii = frePoint
%         Hbdesired(:, ii) = exp(-sqrt(-1)*2*pi*f(ii) * 0.05) * HbdesiredPre(:, ii);
%         para.epsilonB(ii) = 0.01*sqrt(trace(HB(:, :, ii)'*HB(:, :, ii))); % better
%         para.epsilonD(ii) = 0.01*sqrt(trace(HD(:, :, ii)'*HD(:, :, ii)));
%     end
%     for j = frePoint
%         tRB = HB(:, :, j)'*HB(:, :, j);
%         tRD = HD(:, :, j)'*HD(:, :, j);
%         for i = 1:length(rirPerformanceTh) 
%             tdeltaRB = HBe(:, :, j, i)'*HBe(:, :, j, i) - tRB;
%             tdeltaRD = HDe(:, :, j, i)'*HDe(:, :, j, i) - tRD;
%             truegammaB(j, i) = norm(tdeltaRB, 'fro');
%             truegammaD(j, i)= norm(tdeltaRD, 'fro');
%             trueepsilonB(j, i) = norm(HB(:, :, j) - squeeze(HBe(:, :, j, i)), 'fro');
%             trueepsilonD(j, i) = norm(HD(:, :, j) - squeeze(HDe(:, :, j, i)), 'fro');
%             computegammaB(j, i) = trueepsilonB(j, i).*trueepsilonB(j, i);
%             computegammaD(j, i) = trueepsilonD(j, i).*trueepsilonD(j, i);
%         end 
%     end
%     deltaHb = zeros(length(rirPerformanceTh), 1);   
%     deltaHd = deltaHb;
%     for j = frePoint
%         for i = 1:length(rirPerformanceTh)
%             deltaHb(i) = norm(HB(:, :, j) - squeeze(HBe(:, :, j, i)), 'fro');
%             deltaHd(i) = norm(HD(:, :, j) - squeeze(HDe(:, :, j, i)), 'fro');
%         end 
% %         para.epsilonB(j) = mean(deltaHb)*(0.02); % better
% %         para.epsilonD(j) = mean(deltaHd)*(0.02);
%         para.epsilonB(j) = mean(deltaHb)*(0.05); % better
%         para.epsilonD(j) = mean(deltaHd)*(0.05); 
%     end
    
%     % ACC-PM method
%     kappa = 0.7;
%     wAP = ACC_PM(HB, HD, Hbdesired, kappa, para);
%     % RPM method
%     for i = frePoint
%         para.ew(i) = norm(wAP(:, i))^2; 
%     end
%     wRAP = NewRAP(HB, HD, Hbdesired, para);
%     % VAST
%     para.V = floor(L/2);
% %     para.mu = 1;
%     q = VAST(HB, HD, Hbdesired, para);
%     wVASTNF = q{1};
%     % worst-case ACC
%     wwcACC = wcACC(HB, HD, Hbdesired, para); % using eig function, wwcACC is normalized
% %     POTDC ACC
%     wPRACC = POTDC_RACC(HB, HD, para);
    % ACC
    wACC = ACC(HB, HD, para);
    % PM
    wPM = PM(HB, HD, Hbdesired, para);
% compute wRAP, real operation
    performanceTimes = 1:length(rirPerformanceTh);
    for j = frePoint
        for i = performanceTimes
            Hbe = HBe(:, :, j, i);
            Hde = HDe(:, :, j, i);
%             PerformanceIndex = CaculateAC_NSDE;
%             [AC_AP(k, j, i), NSDE_AP(k, j, i), AE_AP(k, j, i)] = PerformanceIndex.AC_NSDE(wAP(:, j), Hbdesired(:, j), Hbe, Hde);
%             PerformanceIndex = CaculateAC_NSDE;
%             [AC_RAP(k, j, i), NSDE_RAP(k, j, i), AE_RAP(k, j, i)] = PerformanceIndex.AC_NSDE(wRAP(:, j), Hbdesired(:, j), Hbe, Hde);
%             PerformanceIndex = CaculateAC_NSDE;
%             [AC_VASTNF(k, j, i), NSDE_VASTNF(k, j, i), AE_VASTNF(k, j, i)] = PerformanceIndex.AC_NSDE(wVASTNF(:, j), Hbdesired(:, j), Hbe, Hde);
%             PerformanceIndex = CaculateAC_NSDE;
%             [AC_wcACC(k, j, i), NSDE_wcACC(k, j, i), AE_wcACC(k, j, i)] = PerformanceIndex.AC_NSDE(wwcACC(:, j), Hbdesired(:, j), Hbe, Hde);
% %             PerformanceIndex = CaculateAC_NSDE;
% %             [AC_PRACC(k, j, i), NSDE_PRACC(k, j, i), AE_PRACC(k, j, i)] = PerformanceIndex.AC_NSDE(wPRACC(:, j), Hbdesired(:, j), Hbe, Hde);
%             % AC sup
%             PerformanceIndex = CaculateAC_NSDE;
%             [AC_ACC(k, j, i), NSDE_ACC(k, j, i), AE_ACC(k, j, i)] = PerformanceIndex.AC_NSDE(wACC(:, j), Hbdesired(:, j), HB(:, :, j), HD(:, :, j));
%             % NSDE inf
%             PerformanceIndex = CaculateAC_NSDE;
%             [AC_PM(k, j, i), NSDE_PM(k, j, i), AE_PM(k, j, i)] = PerformanceIndex.AC_NSDE(wPM(:, j), Hbdesired(:, j), HB(:, :, j), HD(:, :, j));
            % AC sup
            PerformanceIndex = CaculateAC_NSDE;
            [AC_ACC1(k, j, i), NSDE_ACC1(k, j, i), AE_ACC1(k, j, i)] = PerformanceIndex.AC_NSDE(wACC(:, j), Hbdesired(:, j), Hbe, Hde);
            % NSDE inf
            PerformanceIndex = CaculateAC_NSDE;
            [AC_PM1(k, j, i), NSDE_PM1(k, j, i), AE_PM1(k, j, i)] = PerformanceIndex.AC_NSDE(wPM(:, j), Hbdesired(:, j), Hbe, Hde);
        end    
    end
end

%% plot results
load('resultsGauNew.mat');
AC_AP = results.AC_AP;
NSDE_AP = results.NSDE_AP;
AE_AP = results.AE_AP;
AC_RAP = results.AC_RAP;
NSDE_RAP = results.NSDE_RAP;
AE_RAP = results.AE_RAP;
AC_VASTNF = results.AC_VASTNF;
NSDE_VASTNF = results.NSDE_VASTNF;
AE_VASTNF = results.AE_VASTNF;
AC_PRACC = results.AC_PRACC;
NSDE_PRACC = results.NSDE_PRACC;
AE_PRACC = results.AE_PRACC;
AC_wcACC = results.AC_wcACC;
NSDE_wcACC = results.NSDE_wcACC;
AE_wcACC = results.AE_wcACC;
AC_ACC = results.AC_ACC;
NSDE_ACC = results.NSDE_ACC;
AE_ACC = results.AE_ACC;
AC_PM = results.AC_PM;
NSDE_PM = results.NSDE_PM;
AE_PM = results.AE_PM;
f = results.parameters.freq;
frePoint = results.parameters.frePoint;
%% GAU
% NSDE
figure(5)
subplot(3,1,1);
% figure('Name', 'Bright zone error with respect to frequency','Position', [656 212 767 396])  
% set(gcf,'DefaultTextInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'Fontname','Times New Roman','FontSize',15)

for jj = 1:5
    switch jj
        case 1
            method = NSDE_PM;
            style = '-';
            color = '#464747';
            linewidth = 1.1;
            Makerline = 'o';
            MarkerIndices = [2:4:51, 52:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5,'Color','k'); hold on;
        case 2
            method = NSDE_PM1;
            style = '-.';
            color = '#FFB6C1';
            linewidth = 1.5;
            Makerline = 'hexagram';
            MarkerIndices = [5:4:51, 58:20:201];
            MarkerSize = 6;
        case 3
            method = NSDE_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.6;
            Makerline = 'x';
            MarkerIndices = [3:4:51, 54:20:201];
            MarkerSize = 7;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5); hold on;
        case 4
            method = NSDE_VASTNF;
            style = '-.';
%             color = [0.3010 0.7450 0.9330];
%             color = '#4169E1';
            color = '#384bfc';
            linewidth = 1.6;
            Makerline = '*';
            MarkerIndices = [4:4:51, 56:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5); hold on;
        case 5
            method = NSDE_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.1;
            Makerline = 'square';
            MarkerIndices = [5:4:51, 58:20:201];
            MarkerSize = 5;
        
    end
        semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
%                 semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), ...
%             'LineStyle', style, 'linewidth', linewidth); hold on;

end
hold off
% xlabel({'$f (Hz)$'},'Interpreter', 'latex', 'FontSize',15); 
% ylabel({'$\varepsilon _{\rm{NSDE}}$ (dB)'},'Interpreter','latex'); 
ylabel({'$\rm{NRE}$ (dB)'},'Interpreter','latex', 'FontSize',15);
grid on;
% legend('ACC-PM ', 'RPM', 'VASTNF', 'wcACC', 'PRACC');
% legend('ACC-PM ', 'RACC-PM', 'VAST-NF', 'PM(True)');
% legend('PM(True)', 'ACC-PM ', 'VAST-NF', 'RACC-PM');
legend('PM(True)', 'PM(Measured)', 'ACC-PM ', 'VAST-NF', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
% set (gcf,'Position', [200,200,750,280])
xlim([100 f(max(frePoint))]);
ylim([-11 1]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
% end

subplot(3,1,2);
% AC
% figure('Name', 'Acoustic contrast with respect to frequency','Position', [656 212 767 396])
% set(gcf,'DefaultTextInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'Fontname','Times New Roman','FontSize',15)

for jj = 1:8
    switch jj
        case 1
            method = AC_ACC;
            style = '-';
            color = '#464747';
            linewidth = 1.1;
            Makerline = 'o';
            MarkerIndices = [2:6:51, 52:20:201];
            MarkerSize = 5;
        case 2
            method = AC_ACC1;
            style = '-.';
            color = '#17c9e8';
            linewidth = 1.5;
            Makerline = '+';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 7;
        case 3
            method = AC_PM1;
            style = '-.';
            color = '#FFB6C1';
            linewidth = 1.5;
            Makerline = 'hexagram';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 6;
        case 4
            method = AC_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.5;
            Makerline = 'x';
            MarkerIndices = [3:6:51, 54:20:201];
            MarkerSize = 7;
        case 5
            method = AC_VASTNF;
            style = '-.';
%             color = [0.3010 0.7450 0.9330];
%             color = '#4169E1';
            color = '#384bfc';
            linewidth = 1.5;
            Makerline = '*';
            MarkerIndices = [4:6:51, 56:20:201];
            MarkerSize = 5;
        case 6
            method = AC_wcACC;
            style = '-.';
            color = '#edb120';
            linewidth = 1.5;
            Makerline = '^';
            MarkerIndices = [5:6:51, 58:20:201];
            MarkerSize = 5;
        case 7
            method = AC_PRACC;
            style = '-.';
            color = '#c547de';
            linewidth = 1.4;
            Makerline = 'diamond';
            MarkerIndices = [6:6:51, 60:20:201];
            MarkerSize = 5;
        case 8
            method = AC_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.1;
            Makerline = 'square';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 5;   
    end
    semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
end
hold off

% semilogx(f(frePoint), AC_AP,'Color', [0 0 0], 'LineStyle',':','linewidth',1.5); hold on;
% semilogx(f(frePoint), AC_RAP,'Color', [0.6350 0.0780 0.1840], 'LineStyle','--','linewidth',1.5); hold on;
% semilogx(freq, mean_AC_robust,'Color', [0 0.4470 0.7410], 'LineStyle','-.','linewidth',1.5); hold on;
% xlabel({'$f (Hz)$'},'Interpreter','latex', 'FontSize',15); 
ylabel({'${\rm{AC}}$ (dB)'},'Interpreter','latex', 'FontSize',15); 
grid on;
% legend('ACC(True)', 'ACC-PM ', 'VAST-NF', 'wc-RACC', 'POTDC-RACC', 'RACC-PM');
legend('ACC(True)', 'ACC(Measured)', 'PM(Measured)', 'ACC-PM ', 'VAST-NF', 'wc-RACC', 'POTDC-RACC', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
% set (gcf,'Position', [200,200,750,280])
xlim([100 f(max(frePoint))]);
ylim([0 25]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
% ylim([5 25]);

subplot(3,1,3);
% AE
% figure('Name', 'Array effort with respect to frequency','Position', [656 212 767 396])  
% set(gcf,'DefaultTextInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'Fontname','Times New Roman','FontSize',15)
for jj = 1:7
    switch jj
        case 1
            method = AE_ACC1;
            style = '-.';
            color = '#17c9e8';
            linewidth = 1.5;
            Makerline = '+';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 7;
        case 2
            method = AE_PM1;
            style = '-.';
            color = '#FFB6C1';
            linewidth = 1.5;
            Makerline = 'hexagram';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 6;
        case 3
            method = AE_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.5;
            Makerline = 'x';
            MarkerIndices = [3:6:51, 54:20:201];
            MarkerSize = 7;
        case 4
            method = AE_VASTNF;
            style = '-.';
            color = '#384bfc';
            linewidth = 1.5;
            Makerline = '*';
            MarkerIndices = [4:6:51, 56:20:201];
            MarkerSize = 5;
        case 5
            method = AE_wcACC;
            style = '-.';
            color = '#edb120';
            linewidth = 1.5;
            Makerline = '^';
            MarkerIndices = [5:6:51, 58:20:201];
            MarkerSize = 5;
        case 6
            method = AE_PRACC;
            style = '-.';
            color = '#c547de';
            linewidth = 1.4;
            Makerline = 'diamond';
            MarkerIndices = [6:6:51, 60:20:201];
            MarkerSize = 5;
        case 7
            method = AE_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.1;
            Makerline = 'square';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 5;
    end
    semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
end
hold off
xlabel({'$f (Hz)$'},'Interpreter','latex', 'FontSize',15); 
ylabel({'$\rm{AE}(dB)$'},'Interpreter','latex', 'FontSize',15); 
grid on;
legend('ACC(Measured)', 'PM(Measured)', 'ACC-PM', 'VAST-NF', 'wc-RACC', 'POTDC-RACC', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
set (gcf,'Position', [300,50,780,840])
xlim([100 f(max(frePoint))]);
ylim([-15 18]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
%% Rev
% NSDE
figure(4)
subplot(3,1,1);
% figure('Name', 'Bright zone error with respect to frequency','Position', [656 212 767 396])  
% set(gcf,'DefaultTextInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'Fontname','Times New Roman','FontSize',15)

for jj = 1:4
    switch jj
        case 1
            method = NSDE_PM;
            style = '-';
            color = '#464747';
            linewidth = 1.2;
            Makerline = 'o';
            MarkerIndices = [2:4:51, 52:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5,'Color','k'); hold on;
        case 2
            method = NSDE_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.6;
            Makerline = 'x';
            MarkerIndices = [3:4:51, 54:20:201];
            MarkerSize = 7;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5); hold on;
        case 3
            method = NSDE_VASTNF;
            style = '-.';
%             color = [0.3010 0.7450 0.9330];
%             color = '#4169E1';
            color = '#384bfc';
            linewidth = 1.6;
            Makerline = '*';
            MarkerIndices = [4:4:51, 56:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5); hold on;
        case 4
            method = NSDE_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.2;
            Makerline = 'square';
            MarkerIndices = [5:4:51, 58:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5,'Color','r'); hold on;
    end
        semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
%                 semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), ...
%             'LineStyle', style, 'linewidth', linewidth); hold on;

end
hold off
% xlabel({'$f (Hz)$'},'Interpreter', 'latex', 'FontSize',15); 
% ylabel({'$\varepsilon _{\rm{NSDE}}$ (dB)'},'Interpreter','latex'); 
ylabel({'$\rm{NSDE}$ (dB)'},'Interpreter','latex', 'FontSize',15);
grid on;
% legend('ACC-PM ', 'RPM', 'VASTNF', 'wcACC', 'PRACC');
% legend('ACC-PM ', 'RACC-PM', 'VAST-NF', 'PM(True)');
legend('PM(True)', 'ACC-PM ', 'VAST-NF', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
% set (gcf,'Position', [200,200,750,280])
xlim([100 f(max(frePoint))]);
ylim([-13 0]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
% end

subplot(3,1,2);
% AC
% figure('Name', 'Acoustic contrast with respect to frequency','Position', [656 212 767 396])
% set(gcf,'DefaultTextInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'Fontname','Times New Roman','FontSize',15)

for jj = 1:6
    switch jj
        case 1
            method = AC_ACC;
            style = '-';
            color = '#464747';
            linewidth = 1.2;
            Makerline = 'o';
            MarkerIndices = [2:6:51, 52:20:201];
            MarkerSize = 5;
        case 2
            method = AC_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.5;
            Makerline = 'x';
            MarkerIndices = [3:6:51, 54:20:201];
            MarkerSize = 7;
        case 3
            method = AC_VASTNF;
            style = '-.';
%             color = [0.3010 0.7450 0.9330];
%             color = '#4169E1';
            color = '#384bfc';
            linewidth = 1.5;
            Makerline = '*';
            MarkerIndices = [4:6:51, 56:20:201];
            MarkerSize = 5;
        case 4
            method = AC_wcACC;
            style = '-.';
            color = '#edb120';
            linewidth = 1.5;
            Makerline = '^';
            MarkerIndices = [5:6:51, 58:20:201];
            MarkerSize = 5;
        case 5
            method = AC_PRACC;
            style = '-.';
            color = '#c547de';
            linewidth = 1.4;
            Makerline = 'diamond';
            MarkerIndices = [6:6:51, 60:20:201];
            MarkerSize = 5;
        case 6
            method = AC_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.2;
            Makerline = 'square';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 5;

    end
    semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
end
hold off

% semilogx(f(frePoint), AC_AP,'Color', [0 0 0], 'LineStyle',':','linewidth',1.5); hold on;
% semilogx(f(frePoint), AC_RAP,'Color', [0.6350 0.0780 0.1840], 'LineStyle','--','linewidth',1.5); hold on;
% semilogx(freq, mean_AC_robust,'Color', [0 0.4470 0.7410], 'LineStyle','-.','linewidth',1.5); hold on;
% xlabel({'$f (Hz)$'},'Interpreter','latex', 'FontSize',15); 
ylabel({'${\rm{AC}}$ (dB)'},'Interpreter','latex', 'FontSize',15); 
grid on;
legend('ACC(True)', 'ACC-PM ', 'VAST-NF', 'wc-RACC', 'POTDC-RACC', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
% set (gcf,'Position', [200,200,750,280])
xlim([100 f(max(frePoint))]);
% ylim([0 25]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
% ylim([5 25]);

subplot(3,1,3);
% AE
% figure('Name', 'Array effort with respect to frequency','Position', [656 212 767 396])  
% set(gcf,'DefaultTextInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex');
% set(gca,'Fontname','Times New Roman','FontSize',15)
for jj = 1:5
    switch jj
        case 1
            method = AE_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.5;
            Makerline = 'x';
            MarkerIndices = [3:6:51, 54:20:201];
            MarkerSize = 7;
        case 2
            method = AE_VASTNF;
            style = '-.';
            color = '#384bfc';
            linewidth = 1.5;
            Makerline = '*';
            MarkerIndices = [4:6:51, 56:20:201];
            MarkerSize = 5;
        case 3
            method = AE_wcACC;
            style = '-.';
            color = '#edb120';
            linewidth = 1.5;
            Makerline = '^';
            MarkerIndices = [5:6:51, 58:20:201];
            MarkerSize = 5;
        case 4
            method = AE_PRACC;
            style = '-.';
            color = '#c547de';
            linewidth = 1.4;
            Makerline = 'diamond';
            MarkerIndices = [6:6:51, 60:20:201];
            MarkerSize = 5;
        case 5
            method = AE_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.2;
            Makerline = 'square';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 5;
%         case 6
%             method = AE_ACC;
%             style = '--';
%         case 7
%             method = AE_PM;
%             style = ':';
    end
    semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
end
hold off
xlabel({'$f (Hz)$'},'Interpreter','latex', 'FontSize',15); 
ylabel({'$\rm{AE}(dB)$'},'Interpreter','latex', 'FontSize',15); 
grid on;
legend('ACC-PM ', 'VAST-NF', 'wc-RACC', 'POTDC-RACC', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
set (gcf,'Position', [300,50,750,840])
xlim([100 f(max(frePoint))]);
ylim([-8 55]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
%%
% NSDE
subplot(3,1,1);
figure('Name', 'Bright zone error with respect to frequency','Position', [656 212 767 396])  
set(gcf,'DefaultTextInterpreter','latex')
set(gca,'TickLabelInterpreter','latex');
set(gca,'Fontname','Times New Roman','FontSize',15)

for jj = 1:4
    switch jj
        case 1
            method = NSDE_PM;
            style = '-';
            color = '#464747';
            linewidth = 1.2;
            Makerline = 'o';
            MarkerIndices = [2:4:51, 52:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5,'Color','k'); hold on;
        case 2
            method = NSDE_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.6;
            Makerline = 'x';
            MarkerIndices = [3:4:51, 54:20:201];
            MarkerSize = 7;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5); hold on;
        case 3
            method = NSDE_VASTNF;
            style = '-.';
%             color = [0.3010 0.7450 0.9330];
%             color = '#4169E1';
            color = '#384bfc';
            linewidth = 1.6;
            Makerline = '*';
            MarkerIndices = [4:4:51, 56:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5); hold on;
        case 4
            method = NSDE_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.2;
            Makerline = 'square';
            MarkerIndices = [5:4:51, 58:20:201];
            MarkerSize = 5;
%             semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), 'LineStyle',style,'linewidth',1.5,'Color','r'); hold on;
    end
        semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
%                 semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), ...
%             'LineStyle', style, 'linewidth', linewidth); hold on;

end
hold off
xlabel({'$f (Hz)$'},'Interpreter', 'latex', 'FontSize',15); 
% ylabel({'$\varepsilon _{\rm{NSDE}}$ (dB)'},'Interpreter','latex'); 
ylabel({'$\rm{NSDE}$ (dB)'},'Interpreter','latex', 'FontSize',15);
grid on;
% legend('ACC-PM ', 'RPM', 'VASTNF', 'wcACC', 'PRACC');
% legend('ACC-PM ', 'RACC-PM', 'VAST-NF', 'PM(True)');
legend('PM(True)', 'ACC-PM ', 'VAST-NF', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
set (gcf,'Position', [200,200,750,280])
xlim([100 f(max(frePoint))]);
ylim([-11 1]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
% end

subplot(3,1,2);
% AC
figure('Name', 'Acoustic contrast with respect to frequency','Position', [656 212 767 396])
set(gcf,'DefaultTextInterpreter','latex')
set(gca,'TickLabelInterpreter','latex');
set(gca,'Fontname','Times New Roman','FontSize',15)

for jj = 1:6
    switch jj
        case 1
            method = AC_ACC;
            style = '-';
            color = '#464747';
            linewidth = 1.2;
            Makerline = 'o';
            MarkerIndices = [2:6:51, 52:20:201];
            MarkerSize = 5;
        case 2
            method = AC_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.5;
            Makerline = 'x';
            MarkerIndices = [3:6:51, 54:20:201];
            MarkerSize = 7;
        case 3
            method = AC_VASTNF;
            style = '-.';
%             color = [0.3010 0.7450 0.9330];
%             color = '#4169E1';
            color = '#384bfc';
            linewidth = 1.5;
            Makerline = '*';
            MarkerIndices = [4:6:51, 56:20:201];
            MarkerSize = 5;
        case 4
            method = AC_wcACC;
            style = '-.';
            color = '#edb120';
            linewidth = 1.5;
            Makerline = '^';
            MarkerIndices = [5:6:51, 58:20:201];
            MarkerSize = 5;
        case 5
            method = AC_PRACC;
            style = '-.';
            color = '#c547de';
            linewidth = 1.4;
            Makerline = 'diamond';
            MarkerIndices = [6:6:51, 60:20:201];
            MarkerSize = 5;
        case 6
            method = AC_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.2;
            Makerline = 'square';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 5;

    end
    semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
end
hold off

% semilogx(f(frePoint), AC_AP,'Color', [0 0 0], 'LineStyle',':','linewidth',1.5); hold on;
% semilogx(f(frePoint), AC_RAP,'Color', [0.6350 0.0780 0.1840], 'LineStyle','--','linewidth',1.5); hold on;
% semilogx(freq, mean_AC_robust,'Color', [0 0.4470 0.7410], 'LineStyle','-.','linewidth',1.5); hold on;
xlabel({'$f (Hz)$'},'Interpreter','latex', 'FontSize',15); 
ylabel({'${\rm{AC}}$ (dB)'},'Interpreter','latex', 'FontSize',15); 
grid on;
legend('ACC(True)', 'ACC-PM ', 'VAST-NF', 'wc-RACC', 'POTDC-RACC', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
set (gcf,'Position', [200,200,750,280])
xlim([100 f(max(frePoint))]);
ylim([0 25]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);
% ylim([5 25]);

subplot(3,1,3);
% AE
figure('Name', 'Array effort with respect to frequency','Position', [656 212 767 396])  
set(gcf,'DefaultTextInterpreter','latex')
set(gca,'TickLabelInterpreter','latex');
set(gca,'Fontname','Times New Roman','FontSize',15)
for jj = 1:5
    switch jj
        case 1
            method = AE_AP;
            style = '-.';
            color = "#60e01b";
            linewidth = 1.5;
            Makerline = 'x';
            MarkerIndices = [3:6:51, 54:20:201];
            MarkerSize = 7;
        case 2
            method = AE_VASTNF;
            style = '-.';
            color = '#384bfc';
            linewidth = 1.5;
            Makerline = '*';
            MarkerIndices = [4:6:51, 56:20:201];
            MarkerSize = 5;
        case 3
            method = AE_wcACC;
            style = '-.';
            color = '#edb120';
            linewidth = 1.5;
            Makerline = '^';
            MarkerIndices = [5:6:51, 58:20:201];
            MarkerSize = 5;
        case 4
            method = AE_PRACC;
            style = '-.';
            color = '#c547de';
            linewidth = 1.4;
            Makerline = 'diamond';
            MarkerIndices = [6:6:51, 60:20:201];
            MarkerSize = 5;
        case 5
            method = AE_RAP;
            style = '-';
            color = 'r';
            linewidth = 1.2;
            Makerline = 'square';
            MarkerIndices = [7:6:51, 58:20:201];
            MarkerSize = 5;
%         case 6
%             method = AE_ACC;
%             style = '--';
%         case 7
%             method = AE_PM;
%             style = ':';
    end
    semilogx(f(frePoint), squeeze(mean(mean(method, 3), 1)), Makerline, "MarkerIndices", MarkerIndices, ...
            'MarkerSize', MarkerSize, 'LineStyle',style,'linewidth',linewidth,'Color', color); hold on;
end
hold off
xlabel({'$f (Hz)$'},'Interpreter','latex', 'FontSize',15); 
ylabel({'$\rm{AE}(dB)$'},'Interpreter','latex', 'FontSize',15); 
grid on;
legend('ACC-PM ', 'VAST-NF', 'wc-RACC', 'POTDC-RACC', 'RACC-PM');
set(gca,'Fontname','Times New Roman','FontSize',12)
set(gca,'linewidth',1.5)
set (gcf,'Position', [200,200,750,280])
xlim([100 f(max(frePoint))]);
ylim([-15 18]);
xticks([100,1000,min(f(max(frePoint)), 4000)]);

%%
for i = frePoint
    ew(i) = norm(wAP(:, i))^2; 
    AEtACC(i) = norm(wACC(:, i))^2; 
    AEtPM(i) = norm(wPM(:, i))^2;
    AEtwcACC(i) = norm(wwcACC(:, i))^2;
end
alpha = squeeze(mean(mean(AC_ACC, 3), 1));
save alphaAddGauFinal.mat alpha;
save ewAddGau.mat ew;
%%
results.AC_AP = AC_AP;
results.NSDE_AP = NSDE_AP;
results.AE_AP = AE_AP;
results.AC_RAP = AC_RAP;
results.NSDE_RAP = NSDE_RAP;
results.AE_RAP = AE_RAP;
results.AC_VASTNF = AC_VASTNF;
results.NSDE_VASTNF = NSDE_VASTNF;
results.AE_VASTNF = AE_VASTNF;
results.AC_wcACC = AC_wcACC;
results.NSDE_wcACC = NSDE_wcACC;
results.AE_wcACC = AE_wcACC;
results.AC_PRACC = AC_PRACC;
results.NSDE_PRACC = NSDE_PRACC;
results.AE_PRACC = AE_PRACC;
results.AC_ACC = AC_ACC;
results.NSDE_ACC = NSDE_ACC;
results.AE_ACC = AE_ACC;
results.AC_PM = AC_PM;
results.NSDE_PM = NSDE_PM;
results.AE_PM = AE_PM;
results.parameters = para;
results.parameters.fSpaceChoose = fSpaceChoose;
results.parameters.PerformanceChoose = PerformanceChoose;
results.parameters.kappa = kappa;
results.remark = {'HbdesiredPre = squeeze(ATF.irTrue.HB(:, 8, fSpaceChoose))', 'better results', 'added 50 Gau'};
save resultsGau.mat results;
%%
AC_PRACCPre = AC_PRACC(1:11, :, :);