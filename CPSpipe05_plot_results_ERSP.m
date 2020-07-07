% CPSpipe05_plot_results_ERSP
addpath('C:\WORK\MATLAB\Utiles');

% flagDraw = true;
flagDraw = false;

flagTB = true; % ERSP trial basis correction
% flagTB = false;

Fs = 250; % sampling frequency (Hz)
freqband = [0 80];
Ns = 426; % number of samples
tms = 1000*(0:Ns-1)/Fs - 500; % epoch time in ms
winsize = 200; % window size (ms)
winsec = [-400:100:-200 0:100:900]';
winsec = [winsec winsec+winsize];
iwinBase = find(winsec(:,1) < 0);
Nw = size(winsec,1); % number of windows
NFFT = Fs*winsize/1000; % number of samples in one window;
winhann = hanning(NFFT);
freq = Fs*(0:(NFFT/2)-1)/NFFT;
indsel = find(freq >= freqband(1) & freq <= freqband(2));
freq = freq(indsel);
Nf = length(freq);
dirname = 'C:\WORK\MATLAB\Papers\RikWakDynFC';
load(fullfile(dirname,'fuzzy_cmeans_partition_myParcelAAL4mm.mat'), ...
    'labels', 'maskROI_refined');

%%
ERSP_perc_ROIs_allsubj = cell(1,16);
for subj = 1:16
    disp(subj);
    
    if flagDraw
        close all; %#ok<UNRCH>
    end
    
    X = cell(1,1,6);
    for sess = 1:6
        tmp = load(fullfile(dirname,'InvSol',...
            sprintf('inverse_solution_sub%02d_run0%d.mat',subj,sess)));
        X{sess} = tmp.X;
    end
    X = single(permute(cell2mat(X), [2 1 3]));
    [Ns, Nd, Nt] = size(X);
    Y = zeros(Nw, Nf, Nd, Nt, 'single');
    for t = 1:Nw
        % Hann-FFT for selected time window
        flag = (tms >= winsec(t,1)) & (tms < winsec(t,2));
        Xw = fft(bsxfun(@times, X(flag,:,:), winhann));
        Xw = Xw(indsel,:,:);
        Y(t,:,:,:) = Xw;
    end
    if ~flagTB
        ERS = mean(abs(Y).^2, 4); % Eq.(1) from GrandChamp and Delorme (2011)
        muBase = mean(ERS(iwinBase,:,:)); % Eq.(3)
        ERSP_perc = 100*bsxfun(@rdivide, ERS, muBase); % Eq.(5)
    else
        muBase_fullTB = mean(abs(Y).^2); % Eq.(9) modified for Full-TB calculation
        Pk = bsxfun(@rdivide, abs(Y).^2, muBase_fullTB); % Eq.(11) modified for Full-TB calculation
        ERS = mean(Pk, 4); % Eq.(1) modified for Full-TB calculation
        muBase = mean(ERS(iwinBase,:,:)); % Eq.(3) modified for Full-TB calculation
        ERSP_perc_fullTB = 100*bsxfun(@rdivide, ERS, muBase); % Eq.(5) modified for Full-TB calculation
    end
    % ERSP_perc(:,ismember(freq,[0 50]),:) = 0;
%     [i,j,k] = ind2sub([Nw Nf Nd], find(ERSP_perc==max(max(max(ERSP_perc)))));
%     disp([i j k ERSP_perc(i,j,k)]);
%     [i,j,k] = ind2sub([Nw Nf Nd], find(ERSP_perc==min(min(min(nonzeros(ERSP_perc)))),1));
%     disp([i j k ERSP_perc(i,j,k)]);
    
    %% Summarise ERSP per ROIs
    nROIs = length(labels);
    ERSP_perc_ROIs = zeros(Nw, Nf, nROIs);
    for it = 1:nROIs
        flag = (maskROI_refined == it);
        if ~flagTB
            ERSP_perc_ROIs(:,:,it) = mean(ERSP_perc(:,:,flag),3);
        else
            ERSP_perc_ROIs(:,:,it) = mean(ERSP_perc_fullTB(:,:,flag),3);
        end
    end
    ERSP_perc_ROIs_allsubj{subj} = ERSP_perc_ROIs;
    
    %% Plot ERSP for all ROIs
    if flagDraw
        figure; %#ok<UNRCH>
        hax = partition_axes(10, 15, [0.02 0.02 0.002 0.002], [0.002 0.002])';
        indsel_axes = [1:133 136:148];
        xtick = winsec(1:Nw,1) + 100;
        for it = 1:nROIs
            axes(hax(indsel_axes(it))); %#ok<*LAXES>
            surf(xtick, freq(2:Nf), ERSP_perc_ROIs(:,2:end,it)');
            view(2); shading interp; axis tight;
            %     imagesc(ERSP_perc_ROIs(:,2:end,it)');
            %     set(gca, 'Ydir', 'normal');
            caxis([50 200]);
            axis off;
            drawnow;
        end
        axis(hax([134 135 149 150]), 'off');
        %     axis(hax([119 148]), 'on');
        %     axis(hax(119), 'XTick', [0 400 800]);
        %     set(hax(148), 'YAxisLocation', 'right');
        axes('Position',[0.88 0.03 0.05 0.15]);
        colormap jet;
        caxis([50 200]);
        colorbar('FontSize',8);
        axis off;
        if ~flagTB
            print('-dpng', '-r300', sprintf('timefreq_ROIs_subj%02d',subj));
        else
            print('-dpng', '-r300', sprintf('timefreqTB_ROIs_subj%02d',subj));
        end
    end
end