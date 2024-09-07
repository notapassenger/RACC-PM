function[w] = PM(Hb1, Hd1, HbDesired1, para)
%% 
    freq = para.freq;
    frePoint = para.frePoint;
    for i = frePoint
%         w0 = 1; % initialize virtual source filter gain(coefficient).
        pdB = HbDesired1(:, i);
%         pdB = exp(-sqrt(-1)*2*pi*freq(i) * 0.01) * HbDesired1(:, i);
%         pdB = HbDesired1(:, i) * w0; % desired sound pressure in bright zone.
        % closed solution
        pdD = zeros(size(Hd1, 1), 1);
        pd = [pdB; pdD];
        Hbd = [squeeze(Hb1(:, :, i)); squeeze(Hd1(:, :, i))];
        % compute Filter coefficient with Inevitable solution method
        R = (Hbd') * Hbd;
        belta = 0;%max(eig(R))*(1e-5);
        w(:, i) = inv(R+belta*eye(size(Hb1, 2))) * (Hbd') * pd;
    end
end

