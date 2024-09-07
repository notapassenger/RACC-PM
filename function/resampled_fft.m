function[Hb, f, rirResampled] = resampled_fft(rir, nDFT, fsResampled, fs, method, directPath)
%% RACC
% %     rirResampled = 3*resample(rir, fsResampled, fs);
% %     rirResampled = resample(rir, fsResampled, fs);
%     rirResampled = rir;
% %     if nargin > 4
% %         rirResampled((directPath + 12):end) = 0;
% %     end
%     Y = fft(rirResampled, nDFT)/nDFT;
%     P1 = 2*Y(1:floor(nDFT/2) + 1);
% %     P1 = Y(1:floor(nDFT/2) + 1);
%     f = fsResampled*(0:floor(nDFT/2))/nDFT; 
%     Hb= P1';
%% RPM
    if method == 0
        rirResampled = 3*resample(rir, fsResampled, fs);
    %     rirResampled = resample(rir, fsResampled, fs);
        if nargin > 5
            rirResampled((directPath + 12):end) = 0;
        end
        Y = fft(rirResampled, nDFT);
    %     P1 = 2*Y(1:floor(nDFT/2) + 1);
        P1 = Y(1:floor(nDFT/2) + 1);
        f = fsResampled*(0:floor(nDFT/2))/nDFT; 
        Hb= P1';
    elseif method == 1
        rirResampled = resample(rir, fsResampled, fs);
        if nargin > 5
            rirResampled((directPath + 12):end) = 0;
        end
        Y = fft(rirResampled, nDFT);
        P1 = 2*Y(1:floor(nDFT/2) + 1);
        f = fsResampled*(0:floor(nDFT/2))/nDFT; 
        Hb= P1';
    else
        disp('input error');
    end
end