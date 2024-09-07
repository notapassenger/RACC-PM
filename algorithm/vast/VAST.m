function [q] = VAST(HB, HD, Hbdesired, para)
    %% Initialization
%     close all;clear;clc
 
    % % Add paths for relevant files
    % addpath(fullfile(pwd,'subcodes'))
    
    % Initialization
%     datafilename = 'srir_fs16kHz_vast_dft__rir_generator.mat';
%     
%     srir = load(fullfile(pwd,'rirdata',datafilename));
    
    
    % Processing details
%     general = srir.general;
    
    % Array geometry
%     array = srir.array;
    
    % Room
%     room = srir.room;
    
    % Zone
%     zone = srir.zone;
    
    % Simulated room impulse responses (RIRs) by the rir_generator
    % Note that rir_generator can be downloaded from the following link:
    % https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator
    % https://github.com/ehabets/RIR-Generator
%     irMeasured = srir.irMeasured;
    
    % The impulse response of the virtual source
    % This is the RIR of the 8th loudspeaker in irMeasured
%     irVirsrc = srir.irVirsrc;
    
    % Some variables
%     varout = srir.varout;
%     
%     clear srir
    fs = para.fs;
%     fSpace = para.fSpace;
    fSpace = 3;
    nDFT = para.nDFT;
    L = para.L;
    freq = para.freq;
    frePoint = para.frePoint;
    pdB = Hbdesired;
%     for i = frePoint
%         pdB(:, i) = exp(-sqrt(-1)*2*pi*freq(i) * 0.01) .* Hbdesired(:, i);
%     end

    ATFs.Hml1 = permute(HB(:, :, frePoint), [3 1 2]);
    ATFs.Hml2 = permute(HD(:, :, frePoint), [3 1 2]);
%     ATFs.Dm1 = permute(pdB, [2 1]);
    ATFs.Dm1 = permute(pdB(:, frePoint), [2 1]);
    ATFs.Dm2 = ATFs.Dm1;

    array.numLoudspk = L;
    zone.number = 2;
    zone.numCtrPtsBr = size(ATFs.Dm1, 2);
    general.lenConFilter = 2*(size(ATFs.Hml1, 1)-1);
    general.fs = fs;
    general.idx.vast_nf = 1;
    general.idx.vast_t = 3;

    varout.Kbins = size(ATFs.Hml1, 1);
    varout.LJ = array.numLoudspk*general.lenConFilter;
    varout.dF = fs*fSpace/nDFT; % e*
    varout.freq = freq;
    varout.fs = fs;
    varout.jn_no = 2;
    varout.nctrpts = size(ATFs.Hml1, 1);
    varout.nfft = nDFT;
    varout.nloudspks = L;

    % The following parameters may not be used
    varout.lenir = 2000;
    varout.multiplierfactor = 20;
    varout.nmethods = 7;
    varout.nzones = zone.number;
    varout.pickidx = 20;
    varout.pidx = 1:varout.pickidx:varout.pickidx*(size(ATFs.Hml1, 1)-1);
    varout.sfactor = 3;
    %%
    % System geometry illustrated in Fig. 4
    % showSystemGeometry(array,zone,room)
    
    % Initialize control filters
    ctrfilt = getCtrfilt(general, varout);
    targetdB = -10;
    for jj = general.idx.vast_nf:general.idx.vast_t
        ctrfilt{jj}.cvxopt_properties.findopt = false;
        ctrfilt{jj}.cvxopt_properties.opttype = 'min_sd';
        ctrfilt{jj}.cvxopt_properties.const = 'nsb';
        ctrfilt{jj}.cvxopt_properties.tarval = 10^(targetdB/10);
        ctrfilt{jj}.mu = 1;
        ctrfilt{jj}.V = ctrfilt{jj}.Vmax;
    end
    
%     % Initialization process for Sec. V-E
%     cutlength = general.lenConFilter;
%     [irMeasured_cut, irVirsrc_cut] = cutRIRs(general, varout, cutlength, ...
%         irMeasured, irVirsrc);
%     
%     iscuttrue = false;
%     if iscuttrue
%         mrirs = irMeasured_cut;
%         drirs = irVirsrc_cut;
%     else
%         mrirs = irMeasured;
%         drirs = irVirsrc;
%     end
    
    
    % Initialization for VAST-NF
    % 1 kHz
    taroption.taridx = 30; % 32-2
    taroption.tarfreq = (taroption.taridx-1)*varout.dF;
    taroption.journal_exp_1 = false;
    
    
    % Initialization for strings in figures
%     str_ref = {'*Note*',...
%         'In this implementation, the simulated RIRs by rir_generator are used.'};
%     str_ref1 = 'Therefore, the results are similar but different from those shown in Fig. ';
%     bgc = [1,1,1,0.8];
    
    % Note ------------------------------------------------------------------ %
    % The initialization process is done.
    % The code below provides the figures in the manuscript.
    % ------------------------------------------------------------------ Note %
    
    
    %%
    % Note ------------------------------------------------------------------ %
    % In the implementation below,
    % the simulated RIRs by rir_generator are used.
    % Therefore, the results are similar but different from those shown in
    % the manuscript.
    % ------------------------------------------------------------------ Note %
    
    % Note ------------------------------------------------------------------ %
    % In this section, the results of the experiments in 
    %   Sec. V-B Trade-off oAC and SDE
    %   Sec. V-D Performance with respect to the two user parameters
    % are shown.
    % Therefore, Figs. 5, 6, and 8 will be plotted.
    % ------------------------------------------------------------------ Note %
%     close all
    exp1_ctrfilt = ctrfilt{general.idx.vast_nf};
    exp1_taroption = taroption;
    % exp1_taroption.journal_exp_1 = True;
    exp1_taroption.journal_exp_1 = false; % yq
    exp1_ctrfilt.incl_dcnyq = true;
    exp1_ctrfilt.V = para.V;
    exp1_ctrfilt.mu = para.mu;
   
    
    % VAST-NF
    % [q] = calculatefVAST(general, array, zone, exp1_ctrfilt, mrirs, drirs, [], 'narrow', exp1_taroption);
    [q] = calculatefVAST(general, array, zone, exp1_ctrfilt, ATFs, [], 'narrow', exp1_taroption);
end

