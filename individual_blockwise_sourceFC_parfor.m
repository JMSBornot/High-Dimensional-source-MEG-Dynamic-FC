networkdir = '\\sceis_cl1fs\shared\Functional Brain Mapping\JoseSanchez\RikWakDynFC';
cd(networkdir);

subj = 3;

% thread_pause = 10; % pause crashed threads for 10 secs

% Np = 1e4; % number of Monte-Carlo replications
% Np = 1e3 + 1;
Np = 1e4 + 1;

% pval = 1e-7; % p-value for suprathreshold ranks
% pval = [1e-7 1e-8 1e-9];
% pval = [1e-5 1e-6 1e-7];
pval = [1e-5 1e-6 1e-7 1e-8 1e-9];

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
Nw = size(winsec,1); % number of windows
NFFT = Fs*winsize/1000; % number of samples in one window;
winhann = hanning(NFFT);
freq = Fs*(0:(NFFT/2)-1)/NFFT;
indsel = find(freq >= freqband(1) & freq <= freqband(2));
freq = freq(indsel);
Nf = length(freq);
freqsel = [5:5:45 55:5:80]; % selected frequency (Hz);
[flag,iloc_freqsel] = ismember(freqsel, freq);
if ~all(flag)
    error('Some selected frequencies do not appear in the computed ones.');
end

Ndip = 16403; % number of dipoles

% isize = [250*ones(1,65) 153];
isize = [100*ones(1,163) 103];

if (Ndip ~= sum(isize))
    error('Incompatible settings!');
end
nblock = length(isize);
Niter = nchoosek(nblock+1,2);

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

%% Run the FC analysis subject-wise
tmpdir = fullfile(networkdir, sprintf('TEMP_Subj%02d',subj));
if exist(tmpdir,'dir')
    Ntrials = 141; % insert the number of known trials here
else
    mkdir(tmpdir);
    fprintf('Subject #%d, partition data into blocks and save blocks to the network...\n', subj);
    Ntrials = zeros(1,6);
    for sess = 1:6
        Xcell = load(fullfile(networkdir,'InvSol',sprintf('inverse_solution_sub%02d_run0%d.mat',subj,sess)));
        Ntrials(sess) = size(Xcell.X,3);
        Xcell = mat2cell(Xcell.X, isize, size(Xcell.X,2), size(Xcell.X,3));
        for it = 1:nblock
            item = Xcell{it};
            save(fullfile(tmpdir,sprintf('session0%d_block%02d.mat',sess,it)), 'item');
        end
    end
    disp("Done!");
    Ntrials = sum(Ntrials);
end
fprintf('Subject #%02d has %d trials.\n', subj, Ntrials);

eicdir = fullfile(networkdir, sprintf('Subj%02d_BlockwiseFC_SignRank',subj));
if ~exist(eicdir,'dir')
    mkdir(eicdir);
end

% Permutation of signs for signrank (paired) stat
Sperm = int8(reshape(randsample([-1 1], Ntrials*Np, true), [Ntrials Np]));
Sperm(:,1) = 1;

% Select supra-threshold values
fname = sprintf('lookupTable_signrank_N%03d.mat',Ntrials);
[rank, pval_exact] = read_lookup_table(fname);
n = floor(length(rank)/2);
lowTh = zeros(1,length(pval));
highTh = zeros(1,length(pval));
for ith = 1:length(pval)
    lowTh(ith) = interp1(log(pval_exact(n+1:end)), rank(n+1:end), log(pval(ith)), 'linear');
    highTh(ith) = interp1(log(pval_exact(1:n)), rank(1:n), log(pval(ith)), 'linear');
end

%% parallel blockwise FC analysis
% for iter = 1:Niter
parfor iter = 1:Niter
    format compact;
    %         tic
    indprev = find(triu(ones(isize(1)),1)); %#ok<PFBNS>
    indlast = find(triu(ones(isize(end)),1));
    % Reading and preparing the data
    itr = itrow(iter);
    itc = itcol(iter);
    
    fname = fullfile(eicdir, sprintf('EIC_rn%03dcn%03d_%dto%dms_%02dHz.mat',itr,itc,winsec(end,1),winsec(end,2),freqsel(end)));
    if ~exist(fname, 'file')
        X = cell(1,1,6);
        Y = cell(1,1,6);
        for sess = 1:6
            tmpX = load(fullfile(tmpdir,sprintf('session0%d_block%02d.mat',sess,itr)));
            X{sess} = tmpX.item;
            if (itr ~= itc)
                tmpY = load(fullfile(tmpdir,sprintf('session0%d_block%02d.mat',sess,itc)));
                Y{sess} = tmpY.item;
            end
        end
        
        X = single(permute(cell2mat(X), [2 3 1]));
        X = detrend(X);
        [~,Ntr,NdX] = size(X);
        if (itr == itc)
            NdY = NdX;
            ncol = 0.5*NdX*(NdX-1);
        else
            Y = single(permute(cell2mat(Y), [2 3 1]));
            Y = detrend(Y);
            NdY = size(Y,3);
            ncol = NdX*NdY;
        end
        if (Ntr ~= Ntrials)
            error('Number of trials mismatch!');
        end
        
        % Sign of measurements
        Sdiff = zeros(length(freqsel), Ntr, ncol, 'int8');
        % Mask of differences(true or false) between succesive measurements
        Mask = false(length(freqsel), Ntr, ncol);
        % Sorting indices of magnitude of measurements
        IndSort = zeros(length(freqsel), Ntr, ncol, 'int32');
        
        % --- Compute EIC for the baseline
        flag = (tms >= winsec(1,1)) & (tms < winsec(1,2)); %#ok<PFBNS>
        Xw = fft(bsxfun(@times, X(flag,:,:), winhann));
        Xw = Xw(indsel,:,:);
        if (itr == itc)
            Yw = conj(Xw);
        else
            Yw = fft(bsxfun(@times, Y(flag,:,:), winhann));
            Yw = conj(Yw(indsel,:,:));
        end
        EIC_baseline = imag(bsxfun(@times, Xw, reshape(Yw,[Nf Ntr 1 NdY])));
        EIC_baseline = reshape(EIC_baseline, [], NdX*NdY);
        % select only block upper-triangular indices
        if (itr == itc)
            if (itr < nblock)
                EIC_baseline = EIC_baseline(:,indprev);
            else
                EIC_baseline = EIC_baseline(:,indlast);
            end
        end
        % complete Hilbert transform calculation
        EIC_baseline = reshape(EIC_baseline, Nf, []);
        EIC_baseline = abs(hilbert(EIC_baseline)).^2;
        EIC_baseline(isnan(EIC_baseline)) = 0;
        EIC_baseline = reshape(EIC_baseline, [Nf Ntr ncol]);
        EIC_baseline = EIC_baseline(iloc_freqsel,:,:);
        
        % --- Computing EIC for post-stimuli windows
        for t = 2:Nw
            % Hann-FFT for selected time window
            flag = (tms >= winsec(t,1)) & (tms < winsec(t,2));
            Xw = fft(bsxfun(@times, X(flag,:,:), winhann));
            Xw = Xw(indsel,:,:);
            if (itr == itc)
                Yw = conj(Xw);
            else
                Yw = fft(bsxfun(@times, Y(flag,:,:), winhann));
                Yw = conj(Yw(indsel,:,:));
            end
            nbase = 0;
            for id = 1:NdY
                % imaginary xspectrum
                if (itr == itc)
                    if (id == 1), continue; end
                    EIC = imag(bsxfun(@times, Xw(:,:,1:id-1), Yw(:,:,id)));
                    ind = nbase + (1:id-1);
                    nbase = nbase + id - 1; % for next iteration
                else
                    EIC = imag(bsxfun(@times, Xw, Yw(:,:,id)));
                    ind = nbase + (1:NdX);
                    nbase = nbase + NdX; % for next iteration
                end
                % complete Hilbert transform calculation
                EIC = reshape(EIC, Nf, []);
                EIC = abs(hilbert(EIC)).^2;
                EIC(isnan(EIC)) = 0;
                EIC = reshape(EIC, Nf, Ntr, []);
                EIC = EIC(iloc_freqsel,:,:);
                % preparing input to compute signrank stat
                EIC = EIC - EIC_baseline(:,:,ind);
                Sdiff(:,:,ind) = int8(sign(EIC));
                [EIC, IndSort(:,:,ind)] = sort(abs(EIC),2);
                Mask(:,2:end,ind) = (diff(EIC,1,2)==0);
            end
            if (nbase ~= ncol)
                error('Unexpected dimension mismatch!');
            end
            % compute signrank stat
            for k = 1:length(freqsel)
                i1 = squeeze(Mask(k,:,:));
                i2 = squeeze(IndSort(k,:,:)) - 1;
                i3 = squeeze(Sdiff(k,:,:));
                cellind = compute_indices_signrank_cppopt(i1, i2, i3, Sperm, [lowTh;highTh]);
                
                fname = fullfile(eicdir, sprintf('EIC_rn%03dcn%03d_%dto%dms_%02dHz.mat',itr,itc,winsec(t,1),winsec(t,2),freqsel(k)));
                save_cellind(fname, cellind, ncol, lowTh, highTh);
                % save_cellrank_cellind(fname, cellrank, cellind, ncol, lowTh, highTh);
                % fprintf("itr=%02d, itc=%02d, TOI=%02d, FOI=%02d, NNZ=%d\n", itr, itc, t, k, sum(sum(cellfun(@length, cellind)))); drawnow;
            end
        end
    end
end