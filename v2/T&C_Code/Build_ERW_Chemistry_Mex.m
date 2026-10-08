function Build_ERW_Chemistry_Mex(build_dir)
%%% Compile the ERW chemistry kernels into MEX files for TC_ERW_KERNEL=mex.
%%%
%%% Builds ERW_Appelo_Seed_Kernel (with ERW_Al_Equilibrium_Target) and
%%% ERW_Charge_Root_Kernel with MATLAB Coder into this folder, as
%%% ERW_Appelo_Seed_Kernel_mex and ERW_Charge_Root_Kernel_mex for this
%%% platform (.mexa64 on Bouchet, .mexw64 on Windows).  A MEX runs only on
%%% the platform and MATLAB release it was built with, so every deployed
%%% root builds its own (build_erw_mex.sbatch at the package root).  The
%%% campaign fingerprints hash the *.m sources, not the MEX files.
%%%
%%% Coder's MEX code calls MATLAB's own scalar math library (mwmathutil:
%%% muDoubleScalarPower, ...) and the build keeps IEEE arithmetic (no
%%% fast-math, no FMA contraction), so the compiled kernels can reproduce
%%% the MATLAB ones bit for bit; ERW_Chemistry_Mex_Check verifies that at
%%% the end of the build and again the first time a run asks for 'mex'.
%%% A manifest beside the MEX files records the SHA-256 of the compiled
%%% sources, the MATLAB release and the platform.
%%%
%%% build_dir   Coder's working folder (default: a new folder in tempdir);
%%%             only the MEX files and the manifest are written here.

code_dir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(build_dir)
    build_dir = fullfile(tempdir,sprintf('erw_mex_build_%s', ...
        char(java.util.UUID.randomUUID)));
end
if ~exist(build_dir,'dir')
    mkdir(build_dir);
end
sources = {'ERW_Appelo_Seed_Kernel.m','ERW_Al_Equilibrium_Target.m', ...
    'ERW_Charge_Root_Kernel.m'};
for k = 1:numel(sources)
    if ~exist(fullfile(code_dir,sources{k}),'file')
        error('ERW:mexBuild','Kernel source %s is missing.',sources{k});
    end
end

cfg = coder.config('mex');
cfg.GenerateReport = false;
cfg.EnableOpenMP = false;              %% one column per worker, as now
cfg.ResponsivenessChecks = false;      %% no Ctrl-C polling in the loops
cfg.IntegrityChecks = true;            %% bounds checks stay on
cfg.EnableJIT = false;                 %% real C, compiled
set_if_present(cfg,'SIMDAcceleration','None');
set_if_present(cfg,'InstructionSetExtensions','None');
set_if_present(cfg,'OptimizeReductions',false);

real_row = coder.typeof(0,[1 Inf],[false true]);
seed_args = {real_row,0,0,coder.typeof(0,[1 4]),0,0,0, ...
    coder.typeof(0,[1 5]),0,coder.typeof(0,[1 6]),coder.typeof(0,[1 3]), ...
    0,0,coder.typeof(0,[1 7]),false};
scan_args = {real_row,0,0,0,0,0,0,0,0,0,0,0,0,0};

old_dir = pwd;
restore = onCleanup(@() cd(old_dir));
cd(build_dir);
addpath(code_dir);
fprintf('Building ERW chemistry MEX files in %s\n  (Coder working folder %s)\n', ...
    code_dir,build_dir);
t0 = tic;
codegen('-config',cfg,'ERW_Appelo_Seed_Kernel','-args',seed_args, ...
    '-d',fullfile(build_dir,'seed'), ...
    '-o',fullfile(code_dir,'ERW_Appelo_Seed_Kernel_mex'));
codegen('-config',cfg,'ERW_Charge_Root_Kernel','-args',scan_args, ...
    '-d',fullfile(build_dir,'scan'), ...
    '-o',fullfile(code_dir,'ERW_Charge_Root_Kernel_mex'));
fprintf('  built in %.0f s\n',toc(t0));
clear ERW_Appelo_Seed_Kernel_mex ERW_Charge_Root_Kernel_mex
rehash path

manifest = fullfile(code_dir,'ERW_Chemistry_Kernels_mex_manifest.txt');
fid = fopen(manifest,'w');
if fid < 0
    error('ERW:mexBuild','Cannot write %s.',manifest);
end
closer = onCleanup(@() fclose(fid));
fprintf(fid,'built %s\nmatlab %s (%s)\ncomputer %s\nmexext %s\n', ...
    char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')),version, ...
    version('-release'),computer,mexext);
compiler = mex.getCompilerConfigurations('C','Selected');
if ~isempty(compiler)
    fprintf(fid,'compiler %s %s\n',compiler(1).Name,compiler(1).Version);
end
for k = 1:numel(sources)
    fprintf(fid,'%s  %s\n',file_sha256(fullfile(code_dir,sources{k})),sources{k});
end
clear closer

ERW_Chemistry_Mex_Check();
fprintf('PASS: ERW chemistry MEX files built and bit-identical to the MATLAB kernels on the reference steps (%s)\n', ...
    manifest);
end

function set_if_present(cfg,name,value)
if isprop(cfg,name)
    cfg.(name) = value;
end
end

function hex = file_sha256(file_path)
fid = fopen(file_path,'r');
if fid < 0
    error('ERW:mexBuild','Cannot read %s.',file_path);
end
bytes = fread(fid,Inf,'*uint8');
fclose(fid);
digest = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
digest.update(bytes);
hash = typecast(digest.digest(),'uint8');
hex = lower(reshape(dec2hex(hash,2).',1,[]));
end
