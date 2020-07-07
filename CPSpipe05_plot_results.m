function [ConnAll,labels] = CPSpipe05_plot_results(subj, Ndip, Thr, Opt)
% %### Examples (Common code):
%
% rootname = 'C:\WORK\MATLAB\Papers\RikWakDynFC\RESULTS';
% fname_cluster = 'C:\\\\WORK\\\\MATLAB\\\\Papers\\\\RikWakDynFC\\\\BrainFC_CPS_MC%dk\\\\Subj%02d\\\\cps_%%dto%%dms_%%02dHz.mat';
% Ndip = 16403;
% pval = [1e-5 1e-6 1e-7 1e-8 1e-9];
% fname_indtriu = 'C:\WORK\MATLAB\Papers\RikWakDynFC\indtriu_13530_16403.mat';
% load(fname_indtriu, 'indtriu');
% Opt.indtriu = indtriu;
% clear indtriu;
% Opt.winsec = [0  200;
%     100  300;
%     200  400;
%     300  500;
%     400  600;
%     500  700;
%     600  800;
%     700  900;
%     800 1000;
%     900 1100];
% Opt.freqsel = [5:5:45 55:5:80];
% Np = 1000; % Monte-Carlo simulations
% % Np = 10000;
% subjects = 1:16;
% Thr = zeros(Np, 2*length(pval), size(Opt.winsec,1), length(Opt.freqsel), length(subjects));
% for subj = subjects
%     Opt.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
%     for k = 1:length(Opt.freqsel)
%         for t = 1:size(Opt.winsec,1)
%             fname = sprintf(Opt.fname_cluster, Opt.winsec(t,1), Opt.winsec(t,2), Opt.freqsel(k));
%             load(fname, 'maxClusterExt');
%             Thr(:,:,t,k,subj) = maxClusterExt;
%         end
%     end
% end
% ThrClusterSize = reshape(Thr, [Np 2 length(pval) size(Opt.winsec,1) length(Opt.freqsel) length(subjects)]);
% ThrClusterSize = permute(ThrClusterSize, [1 3 2 4 5 6]); 
% ThrClusterSize = reshape(ThrClusterSize, Np, length(pval), [], length(subjects));
% ThrClusterSize = max(ThrClusterSize, [], 3);
% ThrClusterSize = squeeze(prctile(ThrClusterSize, 95));
% cutoff = 0.01; % value used to represent only the graph connections that exceed this value.
% Opt.fname_coords = 'C:\WORK\MATLAB\Papers\RikWakDynFC\myParcelAAL4mm.mat';
% Opt.fname_parcel = 'C:\WORK\MATLAB\Papers\RikWakDynFC\fuzzy_cmeans_partition_myParcelAAL4mm.mat';
% Opt.dirout = 'C:\\WORK\\MATLAB\\Papers\\RikWakDynFC\\OUTPUT\\Subj%02d\\%s';
%
% %### Example #1: Plot dipoles-FC for a single subject, suprathreshold value and time-freq bin.
% out = {'time-freq', [1 2]}; % Plot time-freq bin (1,2).
% Opt01 = Opt;
% Opt01.typeFC = 'dipoles'; % Opt01.typeFC = 'ROIs';
% Opt01.typeRun = {'plot' cutoff out};
% % Opt01.typeRun = {'print' cutoff out};
% % Opt01.fileout = 'timefreq_%s_%dto%dms_%dHz.png'; % for printing option provide output file pattern.
% Opt01.plotOpt.alpha = 0.2;
% Opt01.plotOpt.nodesize = 2;
% Opt01.plotOpt.edgesize = 1;
% Opt01.plotOpt.axis_status = 'off';
% Opt01.plotOpt.caxis = [0 0.8];
% Opt01.plotOpt.SortMethod = 'depth'; % Opt01.plotOpt.SortMethod = 'child';
% Opt01.indthr = 3; % select the 3rd suprathreshold.
% Opt01.statEffect = 'NEG'; % plot only for negative effect.
% subj = 1; % output for subject #1.
% Opt01.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
% CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt01);
%
% %### Example #2: Plot time-slice (freq-slice) for a single subject and suprathreshold value.
% out = {'time-slice', 3}; % Select time-slice #3.
% % out = {'freq-slice', 2}; % Select freq-slice #2.
% Opt02 = Opt;
% Opt02.typeFC = 'ROIs'; % Opt02.typeFC = 'dipoles';
% Opt02.typeRun = {'plot' cutoff out};
% % Opt02.typeRun = {'print' cutoff out};
% % Opt02.fileout = 'timeslice_%s_%dto%dms.png';
% % % Opt02.fileout = 'freqslice_%s_%dHz.png';
% Opt02.plotOpt.alpha = 0.2;
% Opt02.plotOpt.nodesize = 2;
% Opt02.plotOpt.edgesize = 1;
% Opt02.plotOpt.axis_status = 'off';
% Opt02.plotOpt.caxis = [0 0.8];
% Opt02.plotOpt.dimontage = [3 5]; % dimension of subplot montage for time-slice.
% % Opt02.plotOpt.dimontage = [3 4]; % dimension of subplot montage for freq-slice.
% Opt02.plotOpt.SortMethod = 'depth'; % Opt02.plotOpt.SortMethod = 'child';
% Opt02.indthr = 3; % select the 3rd suprathreshold.
% Opt02.statEffect = 'NEG'; % plot only for negative effect.
% subj = 1; % output for subject #1.
% Opt02.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
% CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt02);
%
% %### Example #3: Plot time-transition (freq-transition) for a single subject and suprathreshold value.
% out = {'time-trans', 1:10, 1.5}; % transition of time slices 1:10 with pause for one and half seconds.
% % out = {'freq-trans', 1:15, 1.5}; % transition of freq slices 1:15 with pause for one and half seconds.
% % out = {'time-trans', 1:10, 'pause'}; % transition of time slices 1:10 with pause until user press any key.
% % out = {'time-trans', 1:10, 'keyboard'}; % transition of time slices 1:10 with pause until user type dbcont and press ENTER in MATLAB command input.
% Opt03 = Opt;
% Opt03.typeFC = 'ROIs';
% Opt03.typeRun = {'plot' cutoff out};
% % Opt03.typeRun = {'print' cutoff out};
% % Opt03.fileout = 'timetrans_%s_%dto%dms.png';
% % % Opt03.fileout = 'freqtrans_%s_%dHz.png';
% Opt03.plotOpt.alpha = 0.2;
% Opt03.plotOpt.nodesize = 2;
% Opt03.plotOpt.edgesize = 1;
% Opt03.plotOpt.axis_status = 'off';
% Opt03.plotOpt.caxis = [0 0.8];
% Opt03.plotOpt.dimontage = [3 5]; % dimension of subplot montage for time-slice.
% % Opt03.plotOpt.dimontage = [3 4]; % dimension of subplot montage for freq-slice.
% Opt03.plotOpt.SortMethod = 'depth';
% Opt03.indthr = 3; % select the 3rd suprathreshold.
% Opt03.statEffect = 'NEG'; % plot only for negative effect.
% subj = 1; % output for subject #1.
% Opt03.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
% CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt03);
%
% %### Example #4: Display networks of relevant ROIs for a single subject and suprathreshold value.
% out = {'time-trans', 1:10, 12}; % select upto 12 or 13 ROIs in the displayed ROI-ROI connectivity matrices.
% % out = {'freq-trans', 1:15, 12};
% Opt04 = Opt;
% Opt04.typeFC = 'ROIs'; % input 'dipoles' trigger an error
% Opt04.typeRun = {'disp' cutoff out};
% Opt04.indthr = 3; % select the 3rd suprathreshold.
% Opt04.statEffect = 'NEG'; % plot only for negative effect.
% subj = 1; % output for subject #1.
% Opt04.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
% CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt04);
%
% %### Example #5: Print node-degree for all ROIs, subjects and suprathreshold values.
% out = {'ROI-Brain', 'weighted'}; % ROI weighted degree (with rest of Brain's ROIs).
% % out = {'ROI-Brain', 'binary'}; % ROI number of connections (with rest of Brain's ROIs).
% Opt05 = Opt;
% Opt05.typeFC = 'ROIs'; % input 'dipoles' trigger an error
% Opt05.typeRun = {'print' 0 out}; % no cutoff
% Opt05.fileout = 'WeightDegreeAll_%s.png';
% Opt05.indthr = 1:5; % select data for all the suprathresholds, separately.
% Opt05.statEffect = 'BOTH'; % print plots for both positive ('POS') and negative 'NEG' effect.
% for subj = subjects
%     Opt05.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
%     CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt05);
% end
%
% %### Example #6: Only return all ROI-ROI connectivity matrices for a single subject.
% out = {};
% Opt06 = Opt;
% Opt06.typeFC = 'ROIs'; % input 'dipoles' trigger an error
% Opt06.typeRun = {'return' 0 out}; % no cutoff
% Opt06.indthr = 1:5; % select data for all the suprathresholds, separately.
% Opt06.statEffect = 'BOTH';
% subj = 1;
% Opt06.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
% Conn = CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt06);
%
% %### Example #7: Return ROI-ROI connectivity matrices for all subjects and
% %###             compute inter-individual connectiviy distance and imagesc
% %###             of the node degree strength (rows) for all subject (columns).
% out = {};
% Opt07 = Opt;
% Opt07.typeFC = 'ROIs'; % input 'dipoles' trigger an error
% Opt07.typeRun = {'return' 0 out}; % no cutoff
% Opt07.indthr = 1:5; % select data for all the suprathresholds, separately.
% Opt07.statEffect = 'BOTH'; % retrieve connectivity information for both negative and possitive effect.
%
% %# First, compute all the subjects connectivity matrices
% subj = 1;
% Opt07.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
% [Conn, labels] = CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt07);
% ND = zeros(size(Conn,1), length(subjects), size(Conn,5));
% ND(:,subj,:) = sum(sum(sum(Conn>0,4),3),2); % binary
% % ND(:,subj,:) = sum(sum(sum(Conn,4),3),2); % weighted
% ConnAll = cell(1,length(subjects));
% ConnAll{subj} = Conn;
% for subj = subjects(2:end)
%     disp(subj);
%     Opt07.fname_cluster = sprintf(fname_cluster, round(Np/1000), subj);
%     Conn = CPSpipe05_plot_results(subj, Ndip, ThrClusterSize(:,subj), Opt07);
%     ND(:,subj,:) = sum(sum(sum(Conn>0,4),3),2); % binary
%     % ND(:,subj,:) = sum(sum(sum(Conn,4),3),2); % weighted
%     ConnAll{subj} = Conn;
% end
% % save ConnAll ConnAll ND labels -v7.3
%
% %# Now, compute inter-individual connectiviy distance:
% % load('ConnAll.mat', 'ConnAll', 'ND', 'labels');
% ns = length(ConnAll);
% SM = zeros(ns); % Similarity matrix
%
% %# Use either chi-squared,
% % for it1 = 1:ns-1
% %     disp(it1)
% %     x = (ConnAll{it1}(:,:,:,:,1) > 0); % using conn for the less conservative threshold (10^-5)
% %     x = x(:);
% %     for it2 = it1+1:ns
% %         y = (ConnAll{it2}(:,:,:,:,1) > 0);
% %         y = y(:);
% %         [table,chi2,p] = crosstab(x, y);
% %         SM(it1,it2) = chi2;
% %     end
% % end
%
% %# or the correlation criteria
% nROI = size(ConnAll{1}(:,:,:,:,1), 1);
% flag = (triu(ones(nROI)) > 0);
% for it1 = 1:ns-1
%     disp(it1)
%     x = reshape(ConnAll{it1}(:,:,:,:,1), nROI^2, []);
%     x = x(flag,:);
%     for it2 = it1+1:ns
%         y = reshape(ConnAll{it2}(:,:,:,:,1), nROI^2, []);
%         y = y(flag,:);
%         tmp = corrcoef(x(:), y(:));
%         SM(it1,it2) = tmp(1,2);
%     end
% end
%
% % Dist = 10*15*146*146 - SM'; % for chi-squared
% Dist = 1 - SM'; % Dist values in [0; 2] range for corr criteria
%
% figure; imagesc(SM + SM'); axis image; colorbar
% set(gca, 'XTick', 1:16, 'YTick', 1:16);
% set(gca, 'FontSize', 14);
% set(gca, 'XAxisLocation', 'top');
% % print('-dpng', '-r300', 'hist_clustering');
% Dist = Dist(tril(ones(ns),-1)>0);
% Z = linkage(Dist);
% figure; h = dendrogram(Z);
% set(h, 'LineWidth', 2);
% set(gca, 'FontSize', 14);
% % print('-dpng', '-r300', 'dendogram_clustering');
%
% %# Finally, plot the node-strength information using imagesc:
% it = 1; % select data for a particular suprathreshold value
% figure('units', 'normalized', 'position', [0.01 0.3 0.98 0.6]);
% hax1 = axes('Position',[0.03 0.20 0.94 0.35]);
% hax2 = axes('Position',[0.03 0.60 0.94 0.35]);
% axes(hax1);
% imagesc(ND(:,:,it)');
% % imagesc(bsxfun(@rdivide, ND(:,:,it), max(ND(:,:,it)))');
% colorbar('FontSize', 10);
% set(gca, 'XTick', 1:146, 'YTick', 1:16, 'XTickLabel', strrep(labels,'_','\_'), 'FontSize', 6, ...
%     'FontWeight', 'bold', 'XTickLabelRotation', 90, 'YTickLabelRotation', 90, 'Ydir', 'normal');
% ylabel('Desynchronization', 'FontSize', 12);
% axes(hax2);
% imagesc(ND(:,:,it+1)');
% % imagesc(bsxfun(@rdivide, ND(:,:,it+1), max(ND(:,:,it+1)))');
% colorbar('FontSize', 10);
% set(gca, 'XTick', 1:146, 'YTick', 1:16, 'FontSize', 6, 'FontWeight', 'bold', ...
%     'XTickLabelRotation', 90, 'YTickLabelRotation', 90, 'Ydir', 'normal');
% ylabel('Synchronization', 'FontSize', 12);
% % print('-dpng', '-r300', 'NodeDegree_Desynch_Synch');
global defaults

if isempty(defaults)
    addpath('C:\WORK\MATLAB\Utiles');
    addpath('C:\WORK\MATLAB\spm12');
    spm('defaults','EEG');
end

if ~exist('Opt','var')
    error('Please, define Opt struct with all the required parameter settings.');
end
if ~isfield(Opt,'winsec')
    error('Define the option ''Opt.winsec''.');
end
if ~isfield(Opt,'freqsel')
    error('Define the option ''Opt.freqsel''.');
end
if ~isfield(Opt,'anatcol')
    Opt.anatcol = [238 203 193]/255;
end
if ~isfield(Opt,'typeFC')
    error('Define the option ''Opt.typeFC''.');
end
if ~isfield(Opt,'fname_coords')
    error('Define the option ''Opt.fname_coords'' with source coords and labels.');
end
if ~isfield(Opt,'fname_cluster')
    error('Define the option ''Opt.fname_cluster'' with the file pattern of input data.');
end
if strcmp(Opt.typeFC, 'ROIs') && ~isfield(Opt,'fname_parcel')
    error('Define the option ''Opt.fname_parcel'' with refined parcellation.');
end
if ~isfield(Opt,'indthr')
    error('Define ''Opt.indthr'' with the indices of suprathreshold values for analysing corresponding data.');
end
if ~isfield(Opt,'statEffect') || ~ischar(Opt.statEffect)
    error('Define ''Opt.statEffect'' with the type of statistical effects to analyse: ie, positive (''POS''), negative (''NEG'') or both (''BOTH'').');
end
if ~isfield(Opt,'typeRun')
    error('Define the option ''Opt.typeRun''.');
end
if ~iscell(Opt.typeRun) || (length(Opt.typeRun)~=3) || ~ischar(Opt.typeRun{1}) || ...
        ~isscalar(Opt.typeRun{2}) || ~iscell(Opt.typeRun{3})
    error(['Opt.typeRun must be a cell with three items: a string with the run type, ' ...
        'a scalar with the graph connections cutoff, and a cell with further options, ' ...
        'in this order.']);
end
if ~strcmp(Opt.typeRun{1},'plot') && ~strcmp(Opt.typeRun{1},'disp') && ...
        ~strcmp(Opt.typeRun{1},'print') && ~strcmp(Opt.typeRun{1},'return')
    error('Opt.typeRun{1} must be either ''plot'', ''disp'', ''print'' or ''return''.');
end
if strcmp(Opt.typeRun{1},'plot') || strcmp(Opt.typeRun{1},'print')
    if ~isfield(Opt,'plotOpt')
        Opt.plotOpt = [];
    end
    if ismember(Opt.typeRun{3}{1},{'time-slice','time-trans','freq-slice','freq-trans'}) && ...
        ~isfield(Opt.plotOpt,'dimontage')
        error('For this call, define ''Opt.plotOpt.dimontage''. Check the function help.');
    end
end

mesh = gifti(fullfile(defaults.tbx.dir{1}, '..', 'canonical', 'cortex_5124.surf.gii'));
Nw = size(Opt.winsec,1);
% fname_coords = 'C:\WORK\MATLAB\Papers\RikWakDynFC\myParcelAAL4mm.mat';
load(Opt.fname_coords, 'mni_coords', 'pal_labels');

if strcmp(Opt.typeFC, 'ROIs')
    load(Opt.fname_parcel, 'idx', 'centers', 'mask', 'nROIs');
    % Creating the refined mask
    masktmp = zeros(Ndip,1);
    cnt = 0;
    for it = 1:2:nROIs % left and right ROIs are interleaved
        indL = find(mask == it);
        indR = find(mask == it+1);
        if isempty(idx{it})
            masktmp(indL) = cnt+1;
            masktmp(indR) = cnt+2;
            cnt = cnt + 2;
        else
            for k = 1:max(idx{it})
                flag = (idx{it} == k);
                masktmp(indL(flag)) = cnt+1;
                flag = (idx{it+1} == k);
                masktmp(indR(flag)) = cnt+2;
                cnt = cnt + 2;
            end
        end
    end
    labels = cell(cnt,1);
    mniROIs_coords = zeros(cnt,3);
    ndipROIs = zeros(cnt,1);
    cnt = 0;
    for it = 1:2:nROIs % left and right ROIs are interleaved
        indL = find(mask == it);
        indR = find(mask == it+1);
        if isempty(idx{it})
            labels{cnt+1} = pal_labels{it};
            labels{cnt+2} = pal_labels{it+1};
            mniROIs_coords(cnt+1,:) = centers{it};
            mniROIs_coords(cnt+2,:) = centers{it+1};
            ndipROIs(cnt+1) = length(indL);
            ndipROIs(cnt+2) = length(indR);
            cnt = cnt + 2;
        else
            for k = 1:max(idx{it})
                labels{cnt+1} = sprintf('%s_%d', pal_labels{it}, k);
                labels{cnt+2} = sprintf('%s_%d', pal_labels{it+1}, k);
                mniROIs_coords(cnt+1,:) = centers{it}(k,:);
                mniROIs_coords(cnt+2,:) = centers{it+1}(k,:);
                flag = (idx{it} == k);
                ndipROIs(cnt+1) = length(indL(flag));
                flag = (idx{it+1} == k);
                ndipROIs(cnt+2) = length(indR(flag));
                cnt = cnt + 2;
            end
        end
    end
    nROIs = cnt;
    maskROI = masktmp;
    MAXCONN = ndipROIs*ndipROIs';
    MAXCONN = MAXCONN - diag(diag(MAXCONN)) + diag(0.5.*ndipROIs.*(ndipROIs-1));
    xyz_coords = mniROIs_coords;
elseif strcmp(Opt.typeFC, 'dipoles')
    xyz_coords = mni_coords;
else
    error('This option has not be defined.');
end

if strcmp(Opt.statEffect,'NEG')
    indcol = 2*Opt.indthr - 1;
elseif strcmp(Opt.statEffect,'POS')
    indcol = 2*Opt.indthr - 1;
elseif strcmp(Opt.statEffect,'BOTH')
    indcol = [2*Opt.indthr(:)'-1; 2*Opt.indthr(:)'];
    indcol = indcol(:);
else
    error('Define ''Opt.statEffect'' with the type of statistical effects to analyse: ie, positive (''POS''), negative (''NEG'') or both (''BOTH'').');
end
ncol = length(indcol);

switch Opt.typeFC
    case 'dipoles'
        
        return;
    case 'ROIs'
        ConnAll = zeros(nROIs, nROIs, Nw, length(Opt.freqsel), ncol);
        for t = 1:Nw
            for k = 1:length(Opt.freqsel)
                fname = sprintf(Opt.fname_cluster, Opt.winsec(t,1), Opt.winsec(t,2), Opt.freqsel(k));
                load(fname, 'networkInd', 'clust_cell', 'clustext_cell', 'pval');
                for it = 1:ncol
                    indsel = [];
                    ithr = ceil(indcol(it)/2);
                    ind = find(clustext_cell{indcol(it)} > Thr(ithr));
                    for ic = ind
                        indsel = union(indsel, networkInd{indcol(it)}(clust_cell{indcol(it)}==ic));
                    end
                    [ii,jj] = ind2sub([Ndip Ndip], Opt.indtriu(indsel));
                    Conn = sparse(maskROI(ii), maskROI(jj), 1, nROIs, nROIs);
                    Conn = Conn + Conn' - diag(diag(Conn));
                    Conn = Conn./MAXCONN;
                    ConnAll(:,:,t,k,it) = Conn;
                end
            end
        end
end

%% Computations from HERE to the END only applies to ROIs
str = {'negative', 'positive'};
if strcmp(Opt.typeRun{1}, 'print')
    for it = 1:ncol
        ithr = ceil(indcol(it)/2);
        dirtmp = sprintf(Opt.dirout, subj, num2str(pval(ithr)));
        if ~exist(dirtmp, 'dir')
            mkdir(dirtmp);
        end
    end
end
switch Opt.typeRun{1}
    case {'plot' 'print'}
        switch Opt.typeRun{3}{1}
            case 'time-freq'
                t = Opt.typeRun{3}{2}(1);
                k = Opt.typeRun{3}{2}(2);
                for it = 1:ncol
                    figure(100+it);
                    plot_bvolconn(mesh, xyz_coords, Opt.anatcol, ConnAll(:,:,t,k,it), Opt.plotOpt);
                    ithr = ceil(indcol(it)/2);
                    title(sprintf('Suprathreshold = %s: TWOI = %d to %d msec, FOI = %d Hz', ...
                        num2str(pval(ithr)), Opt.winsec(t,1), Opt.winsec(t,2), Opt.freqsel(k)));
                    drawnow;
                    if strcmp(Opt.typeRun{1}, 'print')
                        dirtmp = sprintf(Opt.dirout, subj, num2str(pval(ithr)));
                        print('-dpng', '-r300', fullfile(dirtmp,sprintf(Opt.fileout,...
                            str{mod(indcol(it)-1,2)+1},Opt.winsec(t,1),Opt.winsec(t,2),Opt.freqsel(k))));
                    end
                end
            case 'time-slice'
                t = Opt.typeRun{3}{2};
                for it = 1:ncol
                    figure(100+it);
                    hax = partition_axes(Opt.plotOpt.dimontage(1), Opt.plotOpt.dimontage(2), ...
                        [0.03 0.1 0.05 0.05], [0.02 0.02])';
                    for k = 1:length(Opt.freqsel)
                        axes(hax(k)); %#ok<*LAXES>
                        plot_bvolconn(mesh, xyz_coords, Opt.anatcol, ConnAll(:,:,t,k,it), Opt.plotOpt);
                        title(sprintf('%d Hz', Opt.freqsel(k)), 'FontSize', 6);
                        drawnow;
                    end
                    caxismax = max(max(max(squeeze(ConnAll(:,:,t,:,it)))));
                    if (caxismax > 0)
                        for k = 1:length(Opt.freqsel)
                            axes(hax(k)); %#ok<*LAXES>
                            caxis([0 caxismax]);
                        end
                        axes('Position',[0.7 0.375 0.285 0.25]);
                        colormap jet;
                        caxis([0 caxismax]);
                        colorbar('FontSize',6);
                        axis off;
                        drawnow;
                    end
                    ithr = ceil(indcol(it)/2);
                    set(gcf, 'Name', sprintf('TWOI = %d to %d msec, suprathreshold = %s', ...
                        Opt.winsec(t,1), Opt.winsec(t,2), num2str(pval(ithr))), 'NumberTitle', 'off');
                    if strcmp(Opt.typeRun{1}, 'print')
                        dirtmp = sprintf(Opt.dirout, subj, num2str(pval(ithr)));
                        print('-dpng', '-r300', fullfile(dirtmp,sprintf(Opt.fileout,...
                            str{mod(indcol(it)-1,2)+1},Opt.winsec(t,1),Opt.winsec(t,2))));
                    end
                end
            case 'freq-slice'
                k = Opt.typeRun{3}{2};
                for it = 1:ncol
                    figure(100+it); clf;
                    hax = partition_axes(Opt.plotOpt.dimontage(1), Opt.plotOpt.dimontage(2), ...
                        [0.03 0.1 0.05 0.05], [0.02 0.02])';
                    for t = 1:Nw
                        axes(hax(t));
                        plot_bvolconn(mesh, xyz_coords, Opt.anatcol, ConnAll(:,:,t,k,it), Opt.plotOpt);
                        title(sprintf('%d-%d ms', Opt.winsec(t,1), Opt.winsec(t,2)), 'FontSize', 6);
                        drawnow;
                    end
                    caxismax = max(max(max(ConnAll(:,:,:,k,it))));
                    if (caxismax > 0)
                        for t = 1:Nw
                            axes(hax(t)); %#ok<*LAXES>
                            caxis([0 caxismax]);
                        end
                        axes('Position',[0.7 0.375 0.285 0.25]);
                        colormap jet;
                        caxis([0 caxismax]);
                        colorbar('FontSize',6);
                        axis off;
                        drawnow;
                    end
                    ithr = ceil(indcol(it)/2);
                    set(gcf, 'Name', sprintf('FOI = %d Hz, suprathreshold = %s', ...
                        Opt.freqsel(k), num2str(pval(ithr))), 'NumberTitle', 'off');
                    if strcmp(Opt.typeRun{1}, 'print')
                        dirtmp = sprintf(Opt.dirout, subj, num2str(pval(ithr)));
                        print('-dpng', '-r300', fullfile(dirtmp,sprintf(Opt.fileout,...
                            str{mod(indcol(it)-1,2)+1},Opt.freqsel(k))));
                    end
                end
            case 'time-trans'
                for it = 1:ncol
                    ithr = ceil(indcol(it)/2);
                    for t = Opt.typeRun{3}{2}
                        figure(101); clf;
                        hax = partition_axes(Opt.plotOpt.dimontage(1), Opt.plotOpt.dimontage(2), ...
                            [0.03 0.1 0.05 0.05], [0.02 0.02])';
                        for k = 1:length(Opt.freqsel)
                            axes(hax(k));
                            plot_bvolconn(mesh, xyz_coords, Opt.anatcol, ConnAll(:,:,t,k,it), Opt.plotOpt);
                            title(sprintf('%d Hz', Opt.freqsel(k)), 'FontSize', 6);
                            drawnow;
                        end
                        caxismax = max(max(max(squeeze(ConnAll(:,:,t,:,it)))));
                        if (caxismax > 0)
                            for k = 1:length(Opt.freqsel)
                                axes(hax(k)); %#ok<*LAXES>
                                caxis([0 caxismax]);
                            end
                            axes('Position',[0.7 0.375 0.285 0.25]);
                            colormap jet;
                            caxis([0 caxismax]);
                            colorbar('FontSize',6);
                            axis off;
                            drawnow;
                        end
                        set(gcf, 'Name', sprintf('TWOI = %d to %d msec, suprathreshold = %s', ...
                            Opt.winsec(t,1), Opt.winsec(t,2), num2str(pval(ithr))), 'NumberTitle', 'off');
                        if strcmp(Opt.typeRun{1}, 'print')
                            dirtmp = sprintf(Opt.dirout, subj, num2str(pval(ithr)));
                            print('-dpng', '-r300', fullfile(dirtmp,sprintf(Opt.fileout,...
                                str{mod(indcol(it)-1,2)+1},Opt.winsec(t,1),Opt.winsec(t,2))));
                        end
                    end
                end
            case 'freq-trans'
                for it = 1:ncol
                    ithr = ceil(indcol(it)/2);
                    for k = Opt.typeRun{3}{2}
                        figure(101); clf;
                        hax = partition_axes(Opt.plotOpt.dimontage(1), Opt.plotOpt.dimontage(2), ...
                            [0.03 0.1 0.05 0.05], [0.02 0.02])';
                        for t = 1:Nw
                            axes(hax(t));
                            plot_bvolconn(mesh, xyz_coords, Opt.anatcol, ConnAll(:,:,t,k,it), Opt.plotOpt);
                            title(sprintf('%d to %d msec', Opt.winsec(t,1), Opt.winsec(t,2)), 'FontSize', 6);
                            drawnow;
                        end
                        caxismax = max(max(max(ConnAll(:,:,:,k,it))));
                        if (caxismax > 0)
                            for t = 1:Nw
                                axes(hax(t)); %#ok<*LAXES>
                                caxis([0 caxismax]);
                            end
                            axes('Position',[0.7 0.375 0.285 0.25]);
                            colormap jet;
                            caxis([0 caxismax]);
                            colorbar('FontSize',6);
                            axis off;
                            drawnow;
                        end
                        set(gcf, 'Name', sprintf('FOI = %d Hz, suprathreshold = %s', ...
                            Opt.freqsel(k), num2str(pval(ithr))), 'NumberTitle', 'off');
                        if strcmp(Opt.typeRun{1}, 'print')
                            dirtmp = sprintf(Opt.dirout, subj, num2str(pval(ithr)));
                            print('-dpng', '-r300', fullfile(dirtmp,sprintf(Opt.fileout,...
                                str{mod(indcol(it)-1,2)+1},Opt.freqsel(k))));
                        end
                    end
                end
            case 'ROI-Brain'
                for it = 1:ncol
                    figure('units', 'normalized', 'position', [0.01 0.3 0.98 0.6]);
                    axes('Position',[0.015 0.65 0.98 0.3]);
                    if strcmp(Opt.typeRun{3}{2}, 'weighted')
                        Conn = sum(sum(sum(ConnAll(:,:,:,:,it),4),3),2);
                    elseif strcmp(Opt.typeRun{3}{2}, 'binary')
                        Conn = sum(sum(sum(ConnAll(:,:,:,:,it)>0,4),3),2);
                    else
                        error('Unexpected!');
                    end
                    plot(Conn);
                    ind = find(Conn > 0.05*max(Conn));
                    if ~isempty(ind)
                        set(gca, 'XTick', ind, 'XTickLabel', strrep(labels(ind),'_','\_'));
                        set(gca, 'FontSize', 6, 'FontWeight', 'bold', 'XTickLabelRotation', 90);
                    end
                    drawnow;
                    % Opt05.fileout = 'WeightDegreeAll_%s.png';
                    if strcmp(Opt.typeRun{1}, 'print')
                        ithr = ceil(indcol(it)/2);
                        dirtmp = sprintf(Opt.dirout, subj, num2str(pval(ithr)));
                        print('-dpng', '-r300', fullfile(dirtmp,...
                                sprintf(Opt.fileout,str{mod(indcol(it)-1,2)+1})));
                    end
                end
        end
    case 'disp'
        switch Opt.typeRun{3}{1}
            case 'time-trans'
                for it = 1:ncol
                    ithr = ceil(indcol(it)/2);
                    fprintf(2,'<strong>Suprathreshold = %s: %s effects:</strong>\n', ...
                        num2str(pval(ithr)), str{mod(indcol(it)-1,2)+1});
                    for t = Opt.typeRun{3}{2}
                        fprintf(2,'<strong>TWOI = %d to %d msec</strong>\n', Opt.winsec(t,1), Opt.winsec(t,2));
                        for k = 1:length(Opt.freqsel)
                            fprintf('<strong>FOI = %d Hz</strong>\n', Opt.freqsel(k));
                            Conn = ConnAll(:,:,t,k,it);
                            Conn(Conn < Opt.typeRun{2}) = 0; % threshold by this cutoff value
                            Conn = roundn(Conn, floor(log10(Opt.typeRun{3})));
                            if (nnz(Conn) == 0)
                                disp('This connectivity matrix is empty');
                                continue;
                            end
                            [ii,jj,vv] = find(Conn);
                            tmpjj = sum(Conn)';
                            tmpii = sum(Conn,2);
                            [~,indsel] = sortrows([vv tmpjj(jj)+tmpii(ii)],[-1 -2]);
                            ind = 0;
                            while (ind < length(indsel)) && (length(unique([ii(indsel(1:ind));jj(indsel(1:ind))])) < Opt.typeRun{3}{3})
                                ind = ind + 1;
                            end
                            indsel = unique([ii(indsel(1:ind)); jj(indsel(1:ind))]);
                            rn = strtrim(cellstr(num2str(indsel)));
                            cn = cellfun(@(x) ['kk' x], rn, 'UniformOutput', false);
                            disp(array2table(Conn(indsel,indsel), 'VariableNames', cn, 'RowNames', rn));
                        end
                    end
                end
            case 'freq-trans'
                for it = 1:ncol
                    ithr = ceil(indcol(it)/2);
                    fprintf(2,'<strong>Suprathreshold = %s: %s effects:</strong>\n', ...
                        num2str(pval(ithr)), str{mod(indcol(it)-1,2)+1});
                    for k = Opt.typeRun{3}{2}
                        fprintf(2,'<strong>FOI = %d Hz</strong>\n', Opt.freqsel(k));
                        for t = 1:Nw
                            fprintf('<strong>TWOI = %d to %d msec</strong>\n', Opt.winsec(t,1), Opt.winsec(t,2));
                            Conn = ConnAll(:,:,t,k,it);
                            Conn(Conn < Opt.typeRun{2}) = 0; % threshold by this cutoff value
                            Conn = roundn(Conn, floor(log10(Opt.typeRun{3})));
                            if (nnz(Conn) == 0)
                                disp('This connectivity matrix is empty');
                                continue;
                            end
                            [ii,jj,vv] = find(Conn);
                            tmpjj = sum(Conn)';
                            tmpii = sum(Conn,2);
                            [~,indsel] = sortrows([vv tmpjj(jj)+tmpii(ii)],[-1 -2]);
                            ind = 0;
                            while (ind < length(indsel)) && (length(unique([ii(indsel(1:ind));jj(indsel(1:ind))])) < Opt.typeRun{3}{3})
                                ind = ind + 1;
                            end
                            indsel = unique([ii(indsel(1:ind)); jj(indsel(1:ind))]);
                            rn = strtrim(cellstr(num2str(indsel)));
                            cn = cellfun(@(x) ['kk' x], rn, 'UniformOutput', false);
                            disp(array2table(Conn(indsel,indsel), 'VariableNames', cn, 'RowNames', rn));
                        end
                    end
                end
        end
    case 'return'
        return;
end

% out = {'time-trans', 1:10, 1.5}; % transition of time slices 1:10 with pause for one and half seconds.
% out = {'freq-trans', 1:15, 1.5}; % transition of freq slices 1:15 with pause for one and half seconds.
% out = {'time-trans', 1:10, 'pause'}; % transition of time slices 1:10 with pause until user press any key.
% out = {'time-trans', 1:10, 'keyboard'}; % transition of time slices 1:10 with pause until user type dbcont and press ENTER in MATLAB command input.

% out = {'ROI-Brain', 'weighted'}; % ROI weighted degree (with rest of Brain's ROIs).
% out = {'ROI-Brain', 'binary'}; % ROI number of connections (with rest of Brain's ROIs).

