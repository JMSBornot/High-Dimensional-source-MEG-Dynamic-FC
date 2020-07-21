function run_dataset(subjects)

myCluster = parcluster('HPCServerProfile1');

% subj = 7;
for subj = subjects
    % run blockwise FC measures and stats
    individual_blockwise_sourceFC_parfor;
    % merge together all the blockwise indices
    clustperm_indsel_DynFC_parfor
    % delete temporary directories
    job = createJob(myCluster);
    createTask(job, @rmdir, 0, {{tmpdir,'s'} {eicdir,'s'}});
    submit(job);
end