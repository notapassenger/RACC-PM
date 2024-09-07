function [wwcACC] = wcACC(HB, HD, Hbdesired, para)
    frePoint = para.frePoint;
    epsilonB = para.epsilonB;
    epsilonD = para.epsilonD;
    for i = frePoint
        HB1 = squeeze(HB(:, :, i));
        HD1 = squeeze(HD(:, :, i));
        Rb = HB1'*HB1;
        Rd = HD1'*HD1;
%         gammaB = 2*epsilonB(i)*sqrt(norm(HB1, 'fro')) + epsilonB(i)^2;
%         gammaD = 2*epsilonD(i)*sqrt(norm(HD1, 'fro')) + epsilonD(i)^2;
        gammaB = epsilonB(i)^2;
        gammaD = epsilonD(i)^2;
%         gammaB = 0.5*max(eig(Rb));
%         gammaD = 0.05*(norm(Rd, 'fro'))/100;
        A = inv(Rd + gammaD*eye(size(Rd))) * (Rb - gammaB*eye(size(Rb)));
        [~, wwcACC(:, i)] = MaxEigenvector(A);
    end
end

% function[AC] = ACC(Hb1, Hd1, HbE1, HdE1)
% % DL diagonal loading, 0 or 1AC = zeros(1, size(Hb1, 3));
%     for i = 1:size(Hb1, 3)
%     %     scale = 1e-4;
%         Hbe = squeeze(HbE1(:, :, i));
%         Hde = squeeze(HdE1(:, :, i));
%         % compute Filter coefficient with Inevitable solution method
%         Rb = (squeeze(Hb1(:, :, i))') * squeeze(Hb1(:, :, i));
%         Rd = (squeeze(Hd1(:, :, i))') * squeeze(Hd1(:, :, i));
%         belta = 0; %max(eig(Rd))*(1e-5);
%         R = inv(Rd + belta*eye(size(Hd1, 2))) * (Rb);
%         [~, w] = MaxEigenvector(R);
%     %     wAll(:, i) = w;
%     %     AC(i) = 10*log10(MaxMu);
%         PerformanceIndex = CaculateAC_NSDE;
%         [AC(i)] = PerformanceIndex.AC(w, Hbe, Hde);
%     end
%     
% end