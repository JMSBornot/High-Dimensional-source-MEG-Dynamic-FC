%% load data
load('ConnAll.mat')
load('ERSP_perc_ROIs_allsubj.mat')
clearvars -except ConnAll labels ERSP_perc_ROIs_allsubj
size(ConnAll{1})
size(ERSP_perc_ROIs_allsubj{1})

%% make ERSP data compatible with FC data
nROIs = 146;
rcorr = zeros(1,nROIs);
pval = zeros(1,nROIs);
rlo = zeros(1,nROIs);
rup = zeros(1,nROIs);
Rcorr = zeros(nROIs);
Pval = zeros(nROIs);
Rlo = zeros(nROIs);
Rup = zeros(nROIs);
Rall = zeros(100,1);
Pall = zeros(100,1);
RLOall = zeros(100,1);
RUPall = zeros(100,1);
labelall = cell(100,1);
Nt = 10;
Nf = 15;
freq_ERSP = 0:5:80;
freq_FC = [5:5:45, 55:5:80];
ffreq = ismember(freq_ERSP, freq_FC);
twin_ERSP = [-400:100:-200, 0:100:900];
twin_FC = 0:100:900;
ftwin = ismember(twin_ERSP, twin_FC);
fdiag = (eye(nROIs)>0);

% % indfreqsel = 2:4;
% % indfreqsel = 1:15;
indfreqsel = 2:5;
% % indfreqsel = 13:15;
% % indfreqsel = 15;
% indfreqsel = 1:6;

% % subjects = [12 13 5 14];
% % subjects = [1 3 4 6 7];
% % subjects = [8 9 10 11 16];
subjects = [5 12 15];

% subjects = [2 7 9 15 16];
% subjects = 14:16;

indThrSel = [1 1 2];
% indThrSel = [1 1 1];
% indThrSel = [2 2 2];

indocc = 65:86;
indpar = 93:114;    
    
cmap = [75 0 130; 0 0 255; 0 255 0; 255 255 0; 255 127 0; 255 0 0]/255;

fprintf('\n--- Correlation between ERSP and Short-Range FC desynch:\n');
figure; hax = partition_axes(length(subjects), 1, [0.03 0.01 0.05 0.05], [0.02 0.02])';
cont = 0;
for subj = subjects
    cont = cont + 1;
    % prepare ERSP data in same format as FC data
    X = permute(ERSP_perc_ROIs_allsubj{subj}(ftwin,ffreq,:), [3 1 2]);
    X = X(:,:,indfreqsel);
    % read short-range FC info
    Y = reshape(ConnAll{subj}(:,:,:,:,indThrSel(cont)), nROIs^2, []);
    Y = reshape(Y(fdiag,:), [nROIs Nt Nf]);
    Y = Y(:,:,indfreqsel);
    % corr
    for it = 1:nROIs
        x = X(it,:,:);
        y = Y(it,:,:);
        [R,P,RLO,RUP] = corrcoef(x(:), y(:));
        rcorr(it) = R(1,2);
        pval(it) = P(1,2);
        rlo(it) = RLO(1,2);
        rup(it) = RUP(1,2);
    end
    axes(hax(cont)); %#ok<*LAXES>
    errorbar(1:nROIs, rcorr, rcorr-rlo, rup-rcorr, '.', 'LineWidth', 2, 'Color', [.2 .2 .2]);
    if (cont == length(subjects))
        set(gca, 'XTick', 2:2:146, 'XTickLabelRotation', 90);
    else
        set(gca, 'XTick', 2:2:146, 'XTickLabel', []);
    end
    set(gca, 'YTick', -0.8:0.4:0.8, 'FontSize', 12);
    if (mod(indThrSel(cont),2) == 1)
        axis([0 147 -1.2 0.5]);
    else
        axis([0 147 -0.5 1.2]);
    end
    sz = 5*ones(1,146);
    col = zeros(1,146);
    flag = (pval <= 1e-2 & pval > 1e-4);
    sz(flag) = 10;
    col(flag) = 1;
    flag = (pval <= 1e-4 & pval > 1e-6);
    sz(flag) = 20;
    col(flag) = 2;
    flag = (pval <= 1e-6 & pval > 1e-8);
    sz(flag) = 30;
    col(flag) = 3;
    flag = (pval <= 1e-8 & pval > 1e-10);
    sz(flag) = 40;
    col(flag) = 4;
    flag = (pval <= 1e-10);
    sz(flag) = 50;
    col(flag) = 5;
    hold on; scatter(1:nROIs, rcorr, sz, col, 's', 'filled'); hold off;
    colormap(cmap); caxis([0 5]);
%     print('-dpng', '-r300', 'ERSP_shortFC_correlation_146ROIs');
%     errorbar(1:nROIs, [1 2 1.5], [1 1 0.5], 'o', 'CapSize', 18, 'LineWidth', 2, 'MarkerSize', 12, 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
end

fprintf('\n--- Correlation between ERSP and Avg Long-Range FC desynch:\n');
figure; hax = partition_axes(length(subjects), 1, [0.03 0.01 0.05 0.05], [0.02 0.02])';
cont = 0;
for subj = subjects
    cont = cont + 1;
    % prepare ERSP data in same format as FC data
    X = permute(ERSP_perc_ROIs_allsubj{subj}(ftwin,ffreq,:), [3 1 2]);
    X = X(:,:,indfreqsel);
    % read avg long-range FC info
    Y = reshape(ConnAll{subj}(:,:,:,:,indThrSel(cont)), nROIs^2, []);
    Y(fdiag,:) = 0; % removing the diag (zeroing short-range FC)
    Y = reshape(Y, [nROIs nROIs Nt Nf]);
    Y = squeeze(sum(Y)/(nROIs-1)); % avg long-range ROI's FC to the rest of the brain
    Y = Y(:,:,indfreqsel);
    % corr
    for it = 1:nROIs
        x = X(it,:,:);
        y = Y(it,:,:);
        [R,P,RLO,RUP] = corrcoef(x(:), y(:));
        rcorr(it) = R(1,2);
        pval(it) = P(1,2);
        rlo(it) = RLO(1,2);
        rup(it) = RUP(1,2);
    end
    fprintf('corr (occipital): mean = %.2f, std = %.2f, min = %.2f, max = %.2f\n', ...
        mean(rcorr(indocc)), std(rcorr(indocc)), min(rcorr(indocc)), max(rcorr(indocc)));
    fprintf('corr (parietal): mean = %.2f, std = %.2f, min = %.2f, max = %.2f\n', ...
        mean(rcorr(indpar)), std(rcorr(indpar)), min(rcorr(indpar)), max(rcorr(indpar)));
    fprintf('pval: max_occipital = %.2E, max_parietal = %.2E\n', max(pval(indocc)), max(pval(indpar)));
    axes(hax(cont)); %#ok<*LAXES>
    errorbar(1:nROIs, rcorr, rcorr-rlo, rup-rcorr, '.', 'LineWidth', 2, 'Color', [.2 .2 .2]);
    if (cont == length(subjects))
        set(gca, 'XTick', 2:2:146, 'XTickLabelRotation', 90);
    else
        set(gca, 'XTick', 2:2:146, 'XTickLabel', []);
    end
    set(gca, 'YTick', -0.8:0.4:0.8, 'FontSize', 12);
    if (mod(indThrSel(cont),2) == 1)
        axis([0 147 -1.2 0.5]);
    else
        axis([0 147 -0.5 1.2]);
    end
    sz = 5*ones(1,146);
    col = zeros(1,146);
    flag = (pval <= 1e-2 & pval > 1e-4);
    sz(flag) = 10;
    col(flag) = 1;
    flag = (pval <= 1e-4 & pval > 1e-6);
    sz(flag) = 20;
    col(flag) = 2;
    flag = (pval <= 1e-6 & pval > 1e-8);
    sz(flag) = 30;
    col(flag) = 3;
    flag = (pval <= 1e-8 & pval > 1e-10);
    sz(flag) = 40;
    col(flag) = 4;
    flag = (pval <= 1e-10);
    sz(flag) = 50;
    col(flag) = 5;
    hold on; scatter(1:nROIs, rcorr, sz, col, 's', 'filled'); hold off;
    colormap(cmap); caxis([0 5]);
%     print('-dpng', '-r300', 'ERSP_largeFC_correlation_146ROIs');
end