function[w] = ACC(Hb1, Hd1, para)
% function [q] = VAST(HB, HD, Hbdesired, para)
% DL diagonal loading, 0 or 1AC = zeros(1, size(Hb1, 3));
    frePoint = para.frePoint;
    for i = frePoint
    %     scale = 1e-4; 
        % compute Filter coefficient with Inevitable solution method
        Rb = (squeeze(Hb1(:, :, i))') * squeeze(Hb1(:, :, i));
        Rd = (squeeze(Hd1(:, :, i))') * squeeze(Hd1(:, :, i));
        belta = 0;%max(eig(Rd))*(1e-5);
        R = inv(Rd + belta*eye(size(Hd1, 2))) * (Rb);
        [~, w(:, i)] = MaxEigenvector(R);
    end
    
end


