%% 
close all;clear;clc
addpath('functions');
addpath('toolbox');
%% Array Setup
load("array.mat");
% room
L = [4.5, 4.5, 2.2]; 
% LoundspeakerPosition
s = array.s;
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
% plot_new(a, s, bCtrPtsPositions, dCtrPtsPositions);
% strFigureSet = ["1"; "scatter"; ];
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
% generate true rir matrix, resampled(48k->16k) and FFT
HB = zeros(size(bCtrPtsPositions, 1), size(s, 1), nDFT/2 +1);
HD = HB;
% control points
for T = 1:size(bCtrPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(TT, :), L, beta, n);
        [HB(T, TT, :), f, ~] = resampled_fft(rir, nDFT, fsResampled, fs, method);
    end                           
end
for T = 1:size(dCtrPtsPositions, 1)
    for TT = 1:size(s, 1)
        rir = rir_generator(c, fs, dCtrPtsPositions(T, :), s(TT, :), L, beta, n);
        [HD(T, TT, :), f, ~] = resampled_fft(rir, nDFT, fsResampled, fs, method);
    end                           
end
ATF.irTrue.HB = HB;
ATF.irTrue.HD = HD;
ATF.room.c = c;
ATF.room.L = L;
ATF.room.beta = beta;
ATF.irTrue.f = f;
ATF.irTrue.fs = fs;
ATF.irTrue.fsResampled = fsResampled;
ATF.irTrue.nDFT = nDFT;

% Dark
method = 1;
DisturbanceTimesTh = 50;
HD = zeros(size(bCtrPtsPositions, 1), size(s, 1), 1601, DisturbanceTimesTh);
betaNew = (0.4-0.3)*rand(1, DisturbanceTimesTh) + 0.3;
for i = 1:DisturbanceTimesTh
    dCtrPtsPositions = array.dCtrPtsPositions;
    % control points, dark zone
    for T = 1:size(dCtrPtsPositions, 1)
        for TT = 1:size(s, 1)
            rir = rir_generator(c, fs, dCtrPtsPositions(T, :), s(TT, :), L, betaNew(i), n);
            [a, f, ~] = resampled_fft(rir, nDFT, fsResampled, fs, method);
            HD(T, TT, :, i) = a;
        end                           
    end
end
ATF.irMeasured.HD = HD;

% Bright
HB = zeros(size(bCtrPtsPositions, 1), size(s, 1), 1601, DisturbanceTimesTh);
for i = 1:DisturbanceTimesTh
    % control points, dark zone
    for T = 1:size(bCtrPtsPositions, 1)
        for TT = 1:size(s, 1)
            rir = rir_generator(c, fs, bCtrPtsPositions(T, :), s(TT, :), L, betaNew(i), n);
            [a, f, ~] = resampled_fft(rir, nDFT, fsResampled, fs, method);
            HB(T, TT, :, i) = a;
        end                           
    end
end
ATF.irMeasured.HB = HB;
ATF.irMeasured.RevRange = betaNew;
ATF.remark = {'betaNew = (0.4-0.3)*rand(1, DisturbanceTimesTh) + 0.3;', 'method = 1;DisturbanceTimesTh = 50;'};
save ATFAddRev.mat ATF;
