function varargout = Field_Campaign_Grid(action,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% The field grid of a 35-year field campaign (Field_SSP_2D and
%% Field_Uncertainty_2D), as runtime/Field_2D checks it for the five-year
%% runs: a file saved by Build_Field_Grid_Topology, beside the deployed root
%% (a validated root never gains files), named by at most 32 letters,
%% digits and underscores, which name the campaign.
%%
%%   [file,name] = Field_Campaign_Grid('read',repo_root,prefix)
%%       TC_FIELD_TOPOLOGY, resolved, outside the root, well named
%%   topo = Field_Campaign_Grid('load',file,prefix)
%%       the grid; at least two squares and its build record
%%   Field_Campaign_Grid('provenance',p,file,topo,prefix)
%%       a saved run_provenance names this file, its SHA-256, its
%%       topography switch and its square count
%%   runs = Field_Campaign_Grid('runs',text,calcite_case,prefix)
%%       a run list the controller passed (TC_FIELD_RUNS, comma
%%       separated): runs 1-9, or 3, 6 and 9 in a calcite-case campaign
%%   gib = Field_Campaign_Grid('memory',n_squares,n_days)
%%       the memory a run certainly holds by its end [GiB, Slurm's G]
%%   Field_Campaign_Grid('check_memory',n_squares,n_days,prefix)
%%       refuse a job whose memory (SLURM_MEM_PER_NODE) is below that
%%   kernel = Field_Campaign_Grid('kernel',prefix)
%%       the ERW chemistry kernel the job runs (TC_ERW_KERNEL): fast, the
%%       only one a field run runs (unset reads as fast)
%%   Field_Campaign_Grid('kernel_provenance',p,kernel,prefix)
%%       a saved run_provenance records that kernel
%%
%% prefix is the calling campaign's error identifier prefix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch action
    case 'read'
        [varargout{1},varargout{2}] = read_grid_file(varargin{:});
    case 'load'
        varargout{1} = load_grid(varargin{:});
    case 'provenance'
        assert_grid_provenance(varargin{:});
    case 'runs'
        varargout{1} = run_list(varargin{:});
    case 'memory'
        varargout{1} = memory_floor_gib(varargin{:});
    case 'check_memory'
        check_memory(varargin{:});
    case 'kernel'
        varargout{1} = erw_kernel(varargin{:});
    case 'kernel_provenance'
        assert_kernel_provenance(varargin{:});
    otherwise
        error('Field_Campaign_Grid:UnknownAction','Unknown action %s.',action);
end
end

function [file,name] = read_grid_file(repo_root,prefix)
file = strtrim(getenv('TC_FIELD_TOPOLOGY'));
if isempty(file) || ~isfile(file)
    error([prefix ':MissingTopology'], ...
        'TC_FIELD_TOPOLOGY must name an existing grid file, not ''%s''.',file);
end
[~,attributes] = fileattrib(file);
file = attributes.Name;
[~,root_attributes] = fileattrib(repo_root);
if startsWith(file,[root_attributes.Name filesep])
    error([prefix ':TopologyInsideRoot'], ...
        'Keep the grid file beside the deployed root, not inside it.');
end
[~,name,extension] = fileparts(file);
if ~strcmp(extension,'.mat') || isempty(regexp(name,'^[A-Za-z0-9_]{1,32}$','once'))
    error([prefix ':BadTopologyName'], ...
        ['The grid file must be <name>.mat with a name of at most 32 ' ...
         'letters, digits and underscores; it names the campaign.']);
end
end

function topo = load_grid(file,prefix)
S = struct();
contents = whos('-file',file);
if any(strcmp({contents.name},'topo'))
    S = load(file,'topo');
end
if ~isfield(S,'topo') || ~isstruct(S.topo) || ~isscalar(S.topo) || ...
        ~isfield(S.topo,'field_grid_build') || ~isfield(S.topo,'N_cells') || ...
        S.topo.N_cells < 2
    error([prefix ':BadTopology'], ...
        '%s is not a field grid saved by Build_Field_Grid_Topology.',file);
end
topo = S.topo;
end

function assert_grid_provenance(p,file,topo,prefix)
%% The result names the grid file it ran on and that file's content.
fid = fopen(file,'rb');
if fid < 0
    error([prefix ':TopologyUnreadable'],'Cannot read %s.',file);
end
bytes = fread(fid,Inf,'*uint8');
fclose(fid);
md = javaMethod('getInstance','java.security.MessageDigest','SHA-256');
md.update(typecast(bytes,'int8'));
hex = lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
attached = isfield(topo,'topography_attached') && topo.topography_attached;
ok = isfield(p,'field_topology_file') && strcmp(p.field_topology_file,file) && ...
    isfield(p,'field_topology_file_sha256') && strcmp(p.field_topology_file_sha256,hex) && ...
    isfield(p,'opt_dem') && p.opt_dem == attached && ...
    isfield(p,'n_columns') && p.n_columns == topo.N_cells;
if ~ok
    error([prefix ':WrongGridProvenance'], ...
        'The result does not carry the grid it was validated on: %s.',file);
end
end

function runs = run_list(text,calcite_case,prefix)
text = strtrim(char(text));
if isempty(text) || isempty(regexp(text,'^[1-9](,[1-9])*$','once'))
    error([prefix ':BadRunList'], ...
        'TC_FIELD_RUNS must reach the task as runs separated by commas, not ''%s''.',text);
end
runs = str2double(strsplit(text,','));
allowed = 1:9;
if ~isempty(calcite_case); allowed = [3 6 9]; end
if numel(unique(runs)) ~= numel(runs) || any(~ismember(runs,allowed))
    error([prefix ':BadRunList'], ...
        'The runs %s are not distinct runs of %s.',text,mat2str(allowed));
end
end

function gib = memory_floor_gib(n_squares,n_days)
%% Package y keeps every daily array whole for each square: 672 values a
%% square and day, plus the 12 hourly-derived daily stores (O_d over six
%% active layers, ZWT_d, Lk_d, f_d, q_runon_d, Rh_d, Rd_d) and the 79
%% values of an ERW day's summary (FIELD_RECORDS_DAILY_2026-10-08.md,
%% section 2). The arrays are allocated whole at the start and filled as
%% the run goes, so they are all resident by its end. 4 GiB stands for
%% MATLAB, the thread pool and the per-square state. The first ERW day
%% builds each square's whole-run timetable on its worker, a passing peak
%% this does not count; the pilot's MaxRSS measures it.
n_daily_rows = n_days+1;
gib = 4 + 8*(672+12+79)*n_daily_rows*n_squares/2^30;
end

function kernel = erw_kernel(prefix)
%% Package z's fast ERW chemistry kernels compute the original code's
%% numbers bit for bit (ERW_KERNELS_2026-10-08.md), and a field run runs
%% them only: unset reads as fast, anything else is refused (MATLAB runs no
%% MEX function on the thread pool every field run has). The kernel is in
%% neither the stamp nor the checkpoint identity, so every element checks
%% it here and every result must record it.
kernel = lower(strtrim(getenv('TC_ERW_KERNEL')));
if isempty(kernel); kernel = 'fast'; end
if ~strcmp(kernel,'fast')
    error([prefix ':BadKernel'], ...
        'Field runs use the fast ERW kernels: TC_ERW_KERNEL must be fast or unset, not ''%s''.', ...
        kernel);
end
end

function assert_kernel_provenance(p,kernel,prefix)
if ~isfield(p,'erw_kernel') || ~strcmp(char(p.erw_kernel),kernel)
    error([prefix ':WrongKernel'], ...
        'The result does not record the ERW kernel %s it ran with.',kernel);
end
end

function check_memory(n_squares,n_days,prefix)
floor_gib = memory_floor_gib(n_squares,n_days);
allocated_mib = str2double(strtrim(getenv('SLURM_MEM_PER_NODE')));
if ~(isfinite(allocated_mib) && allocated_mib > 0)
    fprintf('Memory: SLURM_MEM_PER_NODE unset; the end of the run needs at least %.0f GiB.\n', ...
        floor_gib);
    return
end
fprintf('Memory: %.0f GiB reserved; %d squares over %d days hold at least %.0f GiB by the end.\n', ...
    allocated_mib/1024,n_squares,n_days,floor_gib);
if allocated_mib/1024 < floor_gib
    error([prefix ':InsufficientMemory'], ...
        ['%d squares over %d days hold at least %.0f GiB by the end of the ' ...
         'run; %.0f GiB is reserved. Raise TC_FIELD_MEM.'], ...
        n_squares,n_days,floor_gib,allocated_mib/1024);
end
end
