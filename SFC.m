%% This is a MATLAB code repository for the manuscript below. 
% (1) 'A Robust Hybrid Algorithm Using Worst-CaseOptimization for Pesonal Sound
% Zones'.
% (2) sound field control added gaussian noise.
% (3) If it is not convenient to run the calculation procedure code, ...
% you can start directly from step 5.
%% 1.load RIR data
clear 
addpath(genpath(pwd));
ATF = importdata("ATFAddGau.mat");
fSpaceChoose = 2:4:201*4;
HBAll = ATF.irMeasured.HB;
HDAll = ATF.irMeasured.HD;
M = size(HBAll, 1);
L = size(HDAll, 2);
HbdesiredPre = ATF.irDesired.HB(:, fSpaceChoose);
clear ATF;
%% 2.Experiment and parameter setting
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
load("alphaAddGauFinal.mat");
% load("alphaAddRev.mat");
para.alpha = alpha; 
para.mu = 1; 
para.rho = 0.1; % can't equal to 1 

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
AC_ACC1 = AC_AP;
NSDE_ACC1 = AC_AP;
AE_ACC1 = AC_AP;
AC_PM1 = AC_AP;
NSDE_PM1 = AC_AP;
AE_PM1 = AC_AP;

%% 3.compute control filter and performance
for k = PerformanceChoose
    rirComputeTh = k;
    rirPerformanceTh = setdiff(PerformanceChoose, rirComputeTh);
    HB = squeeze(HBAll(:, :, fSpaceChoose, rirComputeTh)); 
    HD = squeeze(HDAll(:, :, fSpaceChoose, rirComputeTh));
    HBe = HBAll(:, :, fSpaceChoose, rirPerformanceTh);
    HDe = HDAll(:, :, fSpaceChoose, rirPerformanceTh);
    Hbdesired = HbdesiredPre;
    for ii = frePoint
        Hbdesired(:, ii) = exp(-sqrt(-1)*2*pi*f(ii) * 0.05) * HbdesiredPre(:, ii);
        para.epsilonB(ii) = 0.01*sqrt(trace(HB(:, :, ii)'*HB(:, :, ii))); 
        para.epsilonD(ii) = 0.01*sqrt(trace(HD(:, :, ii)'*HD(:, :, ii)));
        % rev, posion
%         para.epsilonB(ii) = 0.0001*sqrt(trace(HB(:, :, ii)'*HB(:, :, ii))); 
%         para.epsilonD(ii) = 0.0001*sqrt(trace(HD(:, :, ii)'*HD(:, :, ii)));
    end   
    % ACC-PM method
    kappa = 0.7;
    wAP = ACC_PM(HB, HD, Hbdesired, kappa, para);
    % RACC-PM method
    for i = frePoint
        para.ew(i) = norm(wAP(:, i))^2; 
    end
    wRAP = RACC_PM(HB, HD, Hbdesired, para);
    % VAST-NF
    para.V = floor(L/2);
    q = VAST(HB, HD, Hbdesired, para);
    wVASTNF = q{1};
    % worst-case RACC
    wwcACC = wcACC(HB, HD, Hbdesired, para); % using eig function, wwcACC is normalized
    % POTDC-RACC
    wPRACC = POTDC_RACC(HB, HD, para);
    % ACC
    wACC = ACC(HB, HD, para);
    % PM
    wPM = PM(HB, HD, Hbdesired, para);

    performanceTimes = 1:length(rirPerformanceTh);
    for j = frePoint
        for i = performanceTimes
            Hbe = HBe(:, :, j, i);
            Hde = HDe(:, :, j, i);
            PerformanceIndex = CaculateAC_NSDE;
            [AC_AP(k, j, i), NSDE_AP(k, j, i), AE_AP(k, j, i)] = PerformanceIndex.AC_NSDE(wAP(:, j), Hbdesired(:, j), Hbe, Hde);
            PerformanceIndex = CaculateAC_NSDE;
            [AC_RAP(k, j, i), NSDE_RAP(k, j, i), AE_RAP(k, j, i)] = PerformanceIndex.AC_NSDE(wRAP(:, j), Hbdesired(:, j), Hbe, Hde);
            PerformanceIndex = CaculateAC_NSDE;
            [AC_VASTNF(k, j, i), NSDE_VASTNF(k, j, i), AE_VASTNF(k, j, i)] = PerformanceIndex.AC_NSDE(wVASTNF(:, j), Hbdesired(:, j), Hbe, Hde);
            PerformanceIndex = CaculateAC_NSDE;
            [AC_wcACC(k, j, i), NSDE_wcACC(k, j, i), AE_wcACC(k, j, i)] = PerformanceIndex.AC_NSDE(wwcACC(:, j), Hbdesired(:, j), Hbe, Hde);
            PerformanceIndex = CaculateAC_NSDE;
            [AC_PRACC(k, j, i), NSDE_PRACC(k, j, i), AE_PRACC(k, j, i)] = PerformanceIndex.AC_NSDE(wPRACC(:, j), Hbdesired(:, j), Hbe, Hde);
            % ACC, True
            PerformanceIndex = CaculateAC_NSDE;
            [AC_ACC(k, j, i), NSDE_ACC(k, j, i), AE_ACC(k, j, i)] = PerformanceIndex.AC_NSDE(wACC(:, j), Hbdesired(:, j), HB(:, :, j), HD(:, :, j));
            % PM, True
            PerformanceIndex = CaculateAC_NSDE;
            [AC_PM(k, j, i), NSDE_PM(k, j, i), AE_PM(k, j, i)] = PerformanceIndex.AC_NSDE(wPM(:, j), Hbdesired(:, j), HB(:, :, j), HD(:, :, j));
            % ACC, Measured 
            PerformanceIndex = CaculateAC_NSDE;
            [AC_ACC1(k, j, i), NSDE_ACC1(k, j, i), AE_ACC1(k, j, i)] = PerformanceIndex.AC_NSDE(wACC(:, j), Hbdesired(:, j), Hbe, Hde);
            % PM, Measured
            PerformanceIndex = CaculateAC_NSDE;
            [AC_PM1(k, j, i), NSDE_PM1(k, j, i), AE_PM1(k, j, i)] = PerformanceIndex.AC_NSDE(wPM(:, j), Hbdesired(:, j), Hbe, Hde);
        end    
    end
end
%% 4.save results
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
results.AC_ACC1 = AC_ACC1;
results.NSDE_ACC1 = NSDE_ACC1;
results.AE_ACC1 = AE_ACC1;
results.AC_PM1 = AC_PM1;
results.NSDE_PM1 = NSDE_PM1;
results.AE_PM1 = AE_PM1;
results.parameters = para;
results.parameters.fSpaceChoose = fSpaceChoose;
results.parameters.PerformanceChoose = PerformanceChoose;
results.parameters.kappa = kappa;
currentTime = datetime('now', 'TimeZone', 'local', 'Format', 'yyyy_MM_dd');
results.remark = {'HbdesiredPre = squeeze(ATF.irTrue.HB(:, 8, fSpaceChoose))', currentTime, 'added 50 Gau'};
save results/resultsGaunew.mat results;
%% 5.plot results
% If it is not convenient to run the above code, you can directly load the
% results of our run, then run the next section.

load('resultsGau.mat');
% load('resultsRev(0.3-0.6s).mat');
% load('resultsPosion.mat');
%% 6.
plotResults(results);




