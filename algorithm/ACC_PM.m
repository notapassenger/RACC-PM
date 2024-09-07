function[w] = ACC_PM(Hb1, Hd1, pdB1, kappa, para) 
    frePoint = para.frePoint;
    freq = para.freq;
    for i = frePoint
        Hb = squeeze(Hb1(:, :, i));
        Hd = squeeze(Hd1(:, :, i));
        pdB = pdB1(:, i);
%         pdB = squeeze(exp(-sqrt(-1)*2*pi*freq(i) * 0.01) * pdB1(:, i));
        % compute Filter coefficient with Inevitable solution method
        R = kappa*(Hd')*Hd + (1-kappa)*(Hb')*Hb;
        belta = 0;%max(eig(R))*(1e-5);
    %     [a1, a2] = eig(R);
    %     mma = min(min(a2));
    %     if mma < 0
    %         R = R + 10*abs(mma)*eye(size(R));
    %     end
        w(:, i) = inv(R + belta*eye(size(Hb, 2))) * (1-kappa) * (Hb') * pdB;  
    %     w = (R + belta*eye(size(Hb2, 2))) \ ((1-kappa) * (Hb2') * pdB); 
    %     w = R \ ((1-kappa) * (Hb2') * pdB);
    end
end
