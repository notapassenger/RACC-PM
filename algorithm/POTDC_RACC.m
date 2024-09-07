function[wPRACC] = POTDC_RACC(Hb1, Hd1, para)
    L = size(Hb1, 2);
    scale = 1e3;
    frePoint = para.frePoint;
    epsilonB = para.epsilonB;
    epsilonD = para.epsilonD;
    for i = frePoint
        Rb = squeeze(Hb1(:, :, i))'*scale * squeeze(Hb1(:, :, i))*scale;
        Rd = squeeze(Hd1(:, :, i))'*scale * squeeze(Hd1(:, :, i))*scale;
        e = 1e-12;
%         eta = 0.025 * sqrt(real(trace(Rb))); % pre
        eta = epsilonB(i); % pre
        EigvalueB = max(eig(Rb));
%         gammaD = 0.05 * norm(Rd, 'fro')/1e+2;% pre
%         gammaD = 2*epsilonD(i)*sqrt(norm(Hd1(:, :, i), 'fro')) + epsilonD(i)^2;
        gammaD = epsilonD(i)^2;
        alphaL = real(1 / power((1 - eta / sqrt(EigvalueB)), 2));
        I = eye(size(Rd));
        EigvalueBD = max(eig(inv(Rd + gammaD * I) * Rb));
%         w0 = ones(N, 1);
%         w0 = rand(N, 1, "like", 1i);
        w0 = rand(L, 1);
        alphaU = real(EigvalueBD * (w0') * (Rd + gammaD * I) * w0);
        
        % opt compute
        NoIterations = 20;
        alphaOpt = zeros(1, NoIterations);
        WOpt = zeros(L, L, NoIterations);
        for j = 1:NoIterations  
            % Calculate variables Dopt and αopt by solving (15).
            if j == 1
                alphaC = (alphaU - alphaL) * rand(1,1) + alphaL;
            end
            cvx_begin 
%             cvx_solver mosek
                variable W(L, L) complex 
                variable alphaopt 
                minimize(real(trace((Rd + gammaD * I) * W)))
                subject to
                    real(trace(Rb * W)) == alphaopt;
                    W == semidefinite(L, L);
                    alphaL <= alphaopt;
                    alphaopt <= alphaU;
                    power(eta, 2) * real(trace(W)) + (sqrt(alphaC) - 1) + (alphaopt) * (1 / sqrt(alphaC) - 1) <= 0;
            cvx_end
            WOpt(:, :, j) = W; 
            alphaOpt(j) = alphaopt;
            alphaC = alphaopt;
            if (j >= 2) && (real(trace((Rd + gammaD * I) * WOpt(:, :, j - 1))) - real(trace((Rd + gammaD * I) * WOpt(:, :, j))) <= e)
                break;  
            end
        end 
        wPRACC(:, i) = W(:, 1) / W(1, 1);
    end
end
