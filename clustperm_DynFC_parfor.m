% clustperm_DynFC
% addpath('C:\WORK\MATLAB\Tutorial\MexCompiling');

myCluster = parcluster('HPCServerProfile1');
optCluster = parforOptions(myCluster,'RangePartitionMethod','fixed','SubrangeSize',40);

% subjects = 2:16;
% subjects = 13;
% subjects = 16;
% subjects = 13:14;
% subjects = 15:16;
subjects = 3;

% networkdir = 'C:\WORK\MATLAB\Papers\RikWakDynFC';
networkdir = '\\sceis_cl1fs\shared\Functional Brain Mapping\JoseSanchez\RikWakDynFC';
cd(networkdir);

% Np = 1e3 + 1;
Np = 1e4 + 1;

% pval = [1e-5 1e-6 1e-7];
pval = [1e-5 1e-6 1e-7 1e-8 1e-9];
Nth = 2*length(pval);

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
Nw = size(winsec,1);
freqsel = [5:5:45 55:5:80]; % selected frequency (Hz);

%% creating other auxiliar variables
Ndip = 16403; % number of dipoles
isize = [100*ones(1,163) 103];
if (Ndip ~= sum(isize))
    error('Incompatible settings!');
end
nblock = length(isize);
Niter = nchoosek(nblock+1,2);

% % indices for mapping from FC indices to Ndip x Ndip connectivity matrix
% IND = mat2cell(triu(ones(Ndip),1), isize, isize);
% cont = 0;
% for itr = 1:nblock
%     for itc = itr:nblock
%         cont = cont+1;
%         IND{itr,itc} = cont*IND{itr,itc};
%     end
% end
% IND = cell2mat(IND);
% indtriu = cell(Niter, 1);
% for it = 1:Niter
%     if mod(it,100) == 0
%         disp([it Niter]);
%     end
%     indtriu{it} = find(IND == it);
% end
% indtriu = cell2mat(indtriu);
% % save indtriu_13530_16403 indtriu
load(sprintf('indtriu_%d_%d.mat',Niter,Ndip), 'indtriu');

% mapping from vertex index to its correponding row/column index
ivert2irow = repmat(uint16(1:Ndip)',1,Ndip);
ivert2irow = ivert2irow(indtriu);
ivert2icol = repmat(uint16(1:Ndip),Ndip,1);
ivert2icol = ivert2icol(indtriu);
% Creating connection mapping indices
load('myParcelAAL4mm.mat', 'mni_coords')
spatial_res = 4;

% [~, ~, E, ivertNeighBeg] = compute_localization_indices(mni_coords, spatial_res);
% % [E, ivertNeighBeg] = compute_localization_indices(mni_coords, spatial_res);
[E2, ivertNeighBeg2, E, ivertNeighBeg] = compute_localization_indices(mni_coords, spatial_res);

E = int32(E)-1;
E2 = int32(E2)-1;

ivert2irow = ivert2irow - 1; % indices have to start at 0 like in C/C++
ivert2icol = ivert2icol - 1;

ivertNeighBeg = int32(ivertNeighBeg) - 1;
ivertNeighBeg2 = int32(ivertNeighBeg2) - 1;

clear indtriu mni_coords;

%% Calculating cluster-permutation statistics
disp('Computing clusters for surrogate data ...');
[tlist, klist, slist] = meshgrid(2:Nw, 1:length(freqsel), subjects);
tlist = tlist(:);
klist = klist(:);
slist = slist(:);
Nbins = length(tlist);
% for it = 1:Nbins
for it = 1:Nbins
    if mod(it,20) == 1, fprintf('%d of %d\n',it,Nbins); end
    % main variables
    t = tlist(it);
    k = klist(it);
    subj = slist(it);
    % check directories
    statdir = fullfile(networkdir, sprintf('BrainFC_SignRank_MC%dk',round(Np/1000)), sprintf('Subj%02d',subj));
    fname_indsel = fullfile(statdir, sprintf('BrainFC_%dto%dms_%02dHz.mat',winsec(t,1),winsec(t,2),freqsel(k)));
    if ~exist(statdir,'dir')
        error('File must have been created: %s\n', fname_indsel);
    end
    clusterdir = fullfile(networkdir, sprintf('BrainFC_CPS_MC%dk',round(Np/1000)), sprintf('Subj%02d',subj));
    if ~exist(clusterdir,'dir')
        mkdir(clusterdir);
    end
    % main code
    fname_cluster = fullfile(clusterdir, sprintf('cps_%dto%dms_%02dHz.mat',winsec(t,1),winsec(t,2),freqsel(k)));
    if ~exist(fname_cluster, 'file')
        [cellind, highTh, lowTh] = load_cellind(fname_indsel);
        maxClusterExt = zeros(Np-1, Nth);
        %parfor it1 = 2:Np
        parfor (it1 = 2:Np, optCluster)
            for it2 = 1:Nth
                indsel = cellind{it1,it2};
                if ~isempty(indsel)
                    % indsel = uint32(indsel)-1;
                    indsel = indsel - 1;
                    [clust, clust_ext] = compute_networkcluster_cppopt(E, ivert2irow, ivert2icol, ivertNeighBeg, indsel);
                    %[clust, clust_ext] = compute_networkcluster(E2, ivert2irow, ivert2icol, ivertNeighBeg2, indsel);
                    clust = clust+1;
                    nc = max(clust); % number of clusters
                    clust_ext = clust_ext(1:nc);
                    maxClusterExt(it1-1,it2) = max(clust_ext);
                end
            end
        end
        save_cps_bin(fname_cluster, maxClusterExt);
        %save(fname_cluster, 'maxClusterExt');
    end
end
disp('Done!');

%% Adding cluster info for the original data
disp('Computing clusters for original data ...');
[tlist, klist, slist] = meshgrid(2:Nw, 1:length(freqsel), subjects);
tlist = tlist(:);
klist = klist(:);
slist = slist(:);
Nbins = length(tlist);
for it = 1:Nbins
%     warning off;
    if mod(it,20) == 1, fprintf('%d of %d\n',it,Nbins); end
    % main variables
    t = tlist(it);
    k = klist(it);
    subj = slist(it);
    % check directories
    statdir = fullfile(networkdir, sprintf('BrainFC_SignRank_MC%dk',round(Np/1000)), sprintf('Subj%02d',subj));
    fname_indsel = fullfile(statdir, sprintf('BrainFC_%dto%dms_%02dHz.mat',winsec(t,1),winsec(t,2),freqsel(k)));
    if ~exist(fname_indsel,'file')
        error('File must have been created: %s\n', fname_indsel);
    end
    clusterdir = fullfile(networkdir, sprintf('BrainFC_CPS_MC%dk',round(Np/1000)), sprintf('Subj%02d',subj));
    fname_cluster = fullfile(clusterdir, sprintf('cps_%dto%dms_%02dHz.mat',winsec(t,1),winsec(t,2),freqsel(k)));
    if ~exist(fname_cluster, 'file')
        error('File must have been created: %s\n', fname_cluster);
    end
    % [cellind, highTh, lowTh] = load_cellind(fname_indsel);
    load(fname_indsel, 'cellind', 'highTh', 'lowTh');
    load(fname_cluster, 'maxClusterExt');
    clust_cell = cell(1,Nth);
    clustext_cell = cell(1,Nth);
    networkInd = cell(1,Nth);
    for it2 = 1:Nth
        networkInd{it2} = cellind{1,it2};
        if ~isempty(networkInd{it2})
            % indsel = uint32(networkInd{it2})-1;
            indsel = networkInd{it2} - 1;
            [clust, clust_ext] = compute_networkcluster(E2, ivert2irow, ivert2icol, ivertNeighBeg2, indsel);
            clust = clust+1;
            nc = max(clust); % number of clusters
            clust_ext = clust_ext(1:nc);
            clust_cell{it2} = clust;
            clustext_cell{it2} = clust_ext;
        end
    end
    save(fname_cluster, 'networkInd', 'clust_cell', 'clustext_cell', 'maxClusterExt', 'pval', 'highTh', 'lowTh');
end
disp('Done!');