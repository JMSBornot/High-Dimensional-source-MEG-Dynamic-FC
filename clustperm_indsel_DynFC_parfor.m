% clustperm_indsel_DynFC_parfor

networkdir = '\\sceis_cl1fs\shared\Functional Brain Mapping\JoseSanchez\RikWakDynFC';
cd(networkdir);

subj = 3;

% Np = 1e4; % number of Monte-Carlo replications
% Np = 1e3 + 1;
Np = 1e4 + 1;

eicdir = fullfile(networkdir, sprintf('Subj%02d_BlockwiseFC_SignRank',subj));
% statdir = fullfile(networkdir, 'BrainFC_SignRank', sprintf('Subj%02d',subj));
statdir = fullfile(networkdir, sprintf('BrainFC_SignRank_MC%dk',round(Np/1000)), sprintf('Subj%02d',subj));
if ~exist(eicdir,'dir')
    error('Directory %s must have been created!\n', eicdir);
end
if exist(statdir,'dir')
    error('Directory %s should not exist yet!\n', statdir);
else
    mkdir(statdir);
end

% pval = [1e-7 1e-8 1e-9];
% pval = [1e-5 1e-6 1e-7];
pval = [1e-5 1e-6 1e-7 1e-8 1e-9];

Nth = 2*length(pval);

Fs = 250; % sampling frequency (Hz)
freqband = [0 100];
Ns = 426; % number of samples
tms = 1000*(0:Ns-1)/Fs - 500; % epoch time in ms
winsize = 200; % window size (ms)
winsec = [-300 -100;
             0  200;
           100  300;
           200  400;
           300  500;
           400  600;
           500  700;
           600  800;
           700  900;
           800 1000;
           900 1100];
if any(diff(winsec,1,2) ~= winsize)
    error('The length of subwindows is %d ms.', winsize);
end
Nw = size(winsec,1) - 1; % number of windows (minus baseline)
NFFT = Fs*winsize/1000; % number of samples in one window;
winhann = hanning(NFFT);
freq = Fs*(0:(NFFT/2)-1)/NFFT;
indsel = find(freq >= freqband(1) & freq <= freqband(2));
freq = freq(indsel);
Nf = length(freq);
freqsel = [5:5:45 55:5:80]; % selected frequency (Hz);
[flag,iloc_freqsel] = ismember(freqsel, freq);
if ~all(flag)
    error('Some selected frequenciencies do not appear in the computed ones.');
end

Ndip = 16403; % number of dipoles
nFC = nchoosek(Ndip,2);
% isize = [250*ones(1,65) 153];
isize = [100*ones(1,163) 103];
if (Ndip ~= sum(isize))
    error('Incompatible settings!');
end
nblock = length(isize);
Niter = nchoosek(nblock+1,2);

%%
parfor it = 1:length(freqsel)*Nw
% for it = 1:length(freqsel)*Nw
    % setting matrix block subindices
    itrow = zeros(Niter+1,1); %#ok<*UNRCH>
    itcol = [nblock; zeros(Niter,1)];
    for cont = 1:Niter
        if (itcol(cont) == nblock)
            itrow(cont+1) = itrow(cont) + 1;
            itcol(cont+1) = itrow(cont+1);
        else
            itrow(cont+1) = itrow(cont);
            itcol(cont+1) = itcol(cont) + 1;
        end
    end
    itrow(1) = [];
    itcol(1) = [];

    t = mod(it-1,Nw) + 2; % add 1 more because the first is the baseline window, the other 2,3,..., are the window of interest
    k = ceil(it/Nw);
    cellind = cell(Np, Nth);
    nbase = 0;
    for iter = 1:Niter
        if (mod(iter,1000) == 1)
            fprintf("TOI=%02d, FOI=%02d, iter=%d\n", t-1, k, iter); drawnow;
        end
        itr = itrow(iter);
        itc = itcol(iter);
        fname = fullfile(eicdir, sprintf('EIC_rn%03dcn%03d_%dto%dms_%02dHz.mat',itr,itc,winsec(t,1),winsec(t,2),freqsel(k))); %#ok<PFBNS>
        [cind, ncol, lowTh, highTh] = load_cellind_ncol(fname);
        
%         keyboard
        
        cind = cellfun(@(X) X+nbase, cind, 'UniformOutput', false);
        %cind = cellfun(@(X) uint32(X)+nbase, cind, 'UniformOutput', false);
        %cind = cellfun(@(X) uint64(X)+nbase, cind, 'UniformOutput', false);
        
        cellind = cellfun(@(X,Y) cat(2,X,Y), cellind, cind, 'UniformOutput', false);
        %cellind = cellfun(@(X,Y) cat(2,X,uint32(Y)), cellind, cind, 'UniformOutput', false);
        %cellind = cellfun(@(X,Y) cat(2,X,uint64(Y)), cellind, cind, 'UniformOutput', false);
        
        nbase = nbase + ncol;
    end
    if (nbase ~= nFC)
        error('This is unexpected!');
    end
    fname = fullfile(statdir, sprintf('BrainFC_%dto%dms_%02dHz.mat',winsec(t,1),winsec(t,2),freqsel(k)));
    save_cellind_BrainFC(fname, cellind, lowTh, highTh);
    fprintf("Done!\n"); drawnow;
end