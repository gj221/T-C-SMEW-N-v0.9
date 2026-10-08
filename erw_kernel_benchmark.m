function erw_kernel_benchmark(n_rounds)
%%% Time one ERW_function step under TC_ERW_KERNEL=matlab, fast and mex.
%%%
%%% The step is the step-22053 fixture of test_erw_appelo_exchange_floor
%%% with its residual anions scaled to put pore pH between about 11 and 5.8,
%%% called as both frames call it (suppress_no_converge_warning and
%%% skip_erw_data set).  The modes alternate within every round and each
%%% time is the median over n_rounds (default 15), so other load on the
%%% node affects all three alike.  'lsqnonlin' marks the steps where the
%%% reduced Appelo seed did not pass the gate and the solver ran.
%%% Run from the package root:
%%%   matlab -singleCompThread -batch "erw_kernel_benchmark"

if nargin < 1 || isempty(n_rounds), n_rounds = 15; end
root = fileparts(mfilename('fullpath'));
code_dir = fullfile(root,'v2','T&C_Code');
addpath(code_dir,fullfile(code_dir,'Inputs'),fullfile(code_dir,'Parameters'), ...
    fullfile(code_dir,'Tests'));
old_mode = getenv('TC_ERW_KERNEL');
restore = onCleanup(@() setenv('TC_ERW_KERNEL',old_mode));
warning('off','ERW:NoConverge');

erw_source = fileread(fullfile(code_dir,'ERW_function.m'));
F = load(fullfile(code_dir,'Tests','Fixtures','erw_appelo_al_floor_step22053.mat'),'fixture');
F = F.fixture;
header_start = regexp(erw_source,'ERW_function\(','end','once') + 1;
header_end = regexp(erw_source,'\)\r?\n\r?\n%%%% INPUT','start','once') - 1;
names = strtrim(strsplit(regexprep(erw_source(header_start:header_end),'\.\.\.',''),','));
args = cell(1,numel(names));
for k = 1:numel(names)
    if ~strcmp(names{k},'~'), args{k} = F.inputs.(names{k}); end
end
iAn = find(strcmp(names,'An_tot'));
ibp = find(strcmp(names,'bulk_pH_params'));
bp = args{ibp};
bp.suppress_no_converge_warning = true;
bp.skip_erw_data = true;
args{ibp} = bp;

modes = {'matlab','fast'};
if exist('ERW_Appelo_Seed_Kernel_mex','file') == 3 && ...
        exist('ERW_Charge_Root_Kernel_mex','file') == 3
    modes{end+1} = 'mex';
end
fprintf('ERW_function step time, median of %d interleaved rounds (%s, %s, %d computational thread(s))\n', ...
    n_rounds,version('-release'),computer,maxNumCompThreads);
fprintf('%-6s %-7s %-6s %-9s','An x','pH','passes','lsqnonlin');
time_heads = strcat(modes,'_ms');
gain_heads = strcat(modes(2:end),'_x');
fprintf(' %9s',time_heads{:});
fprintf(' %8s',gain_heads{:});
fprintf('\n');
total = zeros(1,numel(modes));
for m = [1.0 1.15 1.5 2 3]
    a = args;
    a{iAn} = args{iAn}*m;
    t = zeros(n_rounds,numel(modes));
    for q = 1:numel(modes)
        setenv('TC_ERW_KERNEL',modes{q});
        [~,~,B,~,~,V] = ERW_function(a{:});          %% warm-up
    end
    for r = 1:n_rounds
        for q = 1:numel(modes)
            setenv('TC_ERW_KERNEL',modes{q});
            t0 = tic;
            ERW_function(a{:});
            t(r,q) = toc(t0);
        end
    end
    med = median(t,1);
    total = total + med;
    used_trr = contains(V.candidate_summary,'trr');
    fprintf('%-6.2f %-7.3f %-6d %-9s',m,B.pH,B.activity_passes,mat2str(used_trr));
    fprintf(' %9.2f',1e3*med);
    fprintf(' %8.1f',med(1)./med(2:end));
    fprintf('\n');
end
fprintf('%-30s','sum of medians');
fprintf(' %9.2f',1e3*total);
fprintf(' %8.1f',total(1)./total(2:end));
fprintf('\n');
end
