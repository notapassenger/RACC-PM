%% add Gau data
close all;clear;clc

%% Array Setup
load("array.mat");
% room
L = [4.5, 4.5, 2.2]; 
% LoundspeakerPosition
s = array.s;
% s = array.s;
% MicrophonePosition
bCtrPtsPositions = array.bCtrPtsPositions; % Bright
dCtrPtsPositions = array.dCtrPtsPositions; % Dark

% Position grapha
figure(1)
scatter(s(:, 1), s(:, 2), 10, 'filled');
hold on;
scatter(bCtrPtsPositions(:, 1), bCtrPtsPositions(:, 2), 5, 'filled');
hold on;
scatter(dCtrPtsPositions(:, 1), dCtrPtsPositions(:, 2), 5, 'filled');
hold off;
title('ArrayPosition');
grid on;
axis equal;

%% FFT set 
c = 343;
fs = 48000;
fsResampled = fs/3;
beta = 0.3; % RT
nDFT = 3200; % n, DFT
n = 2967*3; % rir length in Time domain
LVirtualSource = 8;

%% Room and ArrayPosition(Place the array in the room)

figure(2)
scatter(s(:, 1), s(:, 2), 10, 'filled');
hold on;
scatter(bCtrPtsPositions(:, 1), bCtrPtsPositions(:, 2), 5, 'filled');
hold on;
scatter(dCtrPtsPositions(:, 1), dCtrPtsPositions(:, 2), 5, 'filled');
hold on;
rectangle('Position',[0 0 L(1) L(2)]);
hold off;
title('Room and ArrayPosition');
grid on;
axis equal;
%% generate desired rir matrix, only Reserved direct path, Resampled, and FFT
directPath = zeros(1, size(bCtrPtsPositions, 1));
HB = [];
RIR = [];
method = 1;
for T = 1:size(bCtrPtsPositions, 1)
    distance = norm(bCtrPtsPositions(T, 1:2) - s(LVirtualSource, 1:2));
    directPath(T) = ceil(distance/c * fsResampled);
end
for T = 1:size(bCtrPtsPositions, 1)
    rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(LVirtualSource, :), L, beta, n);
    [HB(T, :), f,  ~] = resampled_fft(rir, nDFT, fsResampled, fs, method, directPath(T)); 
end
ATF.irDesired.HB = HB;
HB = zeros(size(bCtrPtsPositions, 1), size(s, 1), nDFT/2 +1);
HD = HB;


%%  add gau
% Dark
method = 1;
DisturbanceTimesTh = 50;
snrPre = 15:(25-15)/(DisturbanceTimesTh-1):25;
rowrank = randperm(size(snrPre, 2)); 
snr = snrPre(1, rowrank);

% snr = (25-15).*rand(1, 50) + 15;

HD = zeros(size(bCtrPtsPositions, 1), size(s, 1), 1601, DisturbanceTimesTh);
for i = 1:DisturbanceTimesTh
    for T = 1:size(dCtrPtsPositions, 1)
        for TT = 1:size(s, 1)
            rir = rir_generator(c, fs, dCtrPtsPositions(T, :), s(TT, :), L, beta, n);
            rir = awgn(rir, snr(i),'measured');
            [a, f, ~] = resampled_fft(rir, nDFT, fsResampled, fs, method);
            HD(T, TT, :, i) = a;
        end                           
    end
end
ATF.irMeasured.HD = HD;
% Bright
HB = zeros(size(bCtrPtsPositions, 1), size(s, 1), 1601, DisturbanceTimesTh);
for i = 1:DisturbanceTimesTh
    for T = 1:size(bCtrPtsPositions, 1)
        for TT = 1:size(s, 1)
            rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(TT, :), L, beta, n);
            rir = awgn(rir, snr(i),'measured');
            [a, f, ~] = resampled_fft(rir, nDFT, fsResampled, fs, method);
            HB(T, TT, :, i) = a;
        end                           
    end
end
ATF.irMeasured.HB = HB;
save ATFRandnAddGau.mat ATF;

