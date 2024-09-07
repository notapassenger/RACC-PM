function PerformanceIndex = CaculateAC_NSDE
    PerformanceIndex.AC_NSDE = @AC_NSDE;
    PerformanceIndex.AC = @AC;
%     PerformanceIndex.AC = @AE;
end
function[AC, NSDE, AE] =  AC_NSDE(w, pdB, Hbe, Hde)
    % input: w, pd, HbE, HdE
    Mb = size(Hbe, 1);
    Md = size(Hde, 1);
    % compute AC
    mu = ((Md*(w')*(Hbe')*Hbe*w) / (Mb*(w')*(Hde')*Hde*w));
    AC = real(10*log10(mu));
%     % compute NSDE, normalization
%     PM = (power(norm(Hbe*w-pdB), 2) / power(norm(pdB), 2));
%     NSDE = 10*log10(PM);
    % compute NSDE, normalization, contain dark zone
    pd = [pdB; zeros(Md, 1)];
    He = [Hbe; Hde];
    PM = (power(norm(He*w-pd), 2) / power(norm(pd), 2));
    NSDE = 10*log10(PM);
    % compute AE
    w0 = zeros(size(w));
    w0(8) = 1;
    AE = 10*log10(((w')*w)*((w0')*(Hbe')*Hbe*w0)/((w')*(Hbe')*Hbe*w));
end

function[AC, AE] = AC(w, Hbe, Hde)
    Mb = size(Hbe, 1);
    Md = size(Hde, 1);
    % compute AC
    mu = ((Md*(w')*(Hbe')*Hbe*w) / (Mb*(w')*(Hde')*Hde*w));
    AC = real(10*log10(mu));
    % compute AE
    w0 = zeros(size(w));
    w0(14) = 1;
    AE = 10*log10(((w')*w)*((w0')*(Hbe')*Hbe*w0)/((w')*(Hbe')*Hbe*w));
end
