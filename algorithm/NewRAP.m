function[wRAP] = NewRAP(Hb1, Hd1, pBd, para)
    %
    mu = para.mu; 
    rho = para.rho; % can't equal to 1
    alpha = para.alpha;
    ew = para.ew;
    L = para.L;
    M = para.M;
    frePoint = para.frePoint;
    epsilonB = para.epsilonB;
    epsilonD = para.epsilonD;
%     freq = para.freq;
    %
    for i = frePoint
%         scale = 1e2; % to meet cvx data requirements 
        Hb = squeeze(Hb1(:, :, i));
        Hd = squeeze(Hd1(:, :, i));
        scale = max(max(max(abs(Hb))), max(max(abs(Hd))));
        pBdPre = pBd(:, i);
        yPre = [squeeze(pBdPre); zeros(M, 1)]/scale;
        y = [real(yPre); imag(yPre)];
        Hb = squeeze(Hb1(:, :, i))/scale;
        Hd = squeeze(Hd1(:, :, i))/scale;
        lambda1 = sqrt(1-rho);
        lambda2 = sqrt(mu + rho*alpha(i));
        H = [lambda1*Hb; lambda2*Hd];
        epsilon = sqrt((lambda1*epsilonB(i)/scale)^2 + (lambda2*epsilonD(i)/scale)^2);
%         H = [Hb; sqrt(kappa1+kappa2*Ccon)*Hd];
%         epsilon = sqrt((epsilonB(i)/scale)^2 + (kappa1+kappa2*Ccon)*(epsilonD(i)/scale)^2);
        % cvx
        c = [zeros(2*L, 1); 1; 0];  
        d1 = -y;
        d2 = zeros(2*L, 1);
        d3 = d2;
        e1 = [zeros(2*L, 1); 1; -1];
        e2 = [zeros(2*L, 1); 0; 1/epsilon];
        e3 = zeros((2*L+2), 1);
        f1 = 0;
        f2 = f1;
        f3 = sqrt(ew(i));
        Heb = [real(H), -imag(H); imag(H), real(H)];
        C1 = [Heb, zeros(4*M, 2)];
        C2 = [eye(2*L, 2*L), zeros(2*L, 1), zeros(2*L, 1)];
        C3 = C2;
        cvx_begin 
            variable z(2*L+2)
            minimize (c'*z) 
            subject to
                norm(C1*z + d1) <= e1'*z + f1;
                norm(C2*z + d2) <= e2'*z + f2;
                norm(C3*z + d3) <= e3'*z + f3;
        cvx_end
        wRAP(:, i) = z(1:L) + z(L+1:2*L)*1i; 
    end
end