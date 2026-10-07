function deltas = Extract_CMIP6_Monthly_Deltas( ...
    nexgddp_root,output_file,Lat,Lon,scenario,models,varargin)
% Extract monthly change factors from NEX-GDDP-CMIP6 at one grid cell.
% Returns additive deltas for temperature and relative humidity and
% multiplicative factors for precipitation, shortwave and wind speed.
% nexgddp_root is a local directory tree or an OPeNDAP base URL holding
% files named <var>_day_<model>_<scenario>_<variant>_<grid>_<year>*.nc.
% Deltas are formed between a historical baseline window and a future
% window of the named scenario, both drawn from the same model, so the
% observational baseline period need not match BaselineYears exactly.
%
% Name-value options:
%   BaselineYears    [1985 2014]  CMIP6 historical reference window
%   FutureYears      [2036 2065]  window centred on the ramp end year
%   PrecipQuantiles  []           probabilities for QDM ratios, e.g. 0.05:0.05:0.95
%   Verbose          true

narginchk(6,inf);
opt = parse_options(varargin);

nexgddp_root = char(nexgddp_root);
output_file = char(output_file);
scenario = char(scenario);
if ischar(models)
    models = {models};
end
validateattributes(Lat,{'numeric'},{'scalar','real','finite','>=',-90,'<=',90});
validateattributes(Lon,{'numeric'},{'scalar','real','finite','>=',-180,'<=',180});
if isempty(models) || ~iscellstr(models)
    error('Extract_CMIP6_Monthly_Deltas:BadModels', ...
        'models must be a non-empty cell array of model names.');
end
[output_dir,~,output_ext] = fileparts(output_file);
if ~isempty(output_dir) && ~isfolder(output_dir)
    error('Extract_CMIP6_Monthly_Deltas:MissingOutputDirectory', ...
        'Output directory not found: %s',output_dir);
end
if ~strcmpi(output_ext,'.mat')
    error('Extract_CMIP6_Monthly_Deltas:BadOutputExtension', ...
        'Output file must have a .mat extension.');
end

% NEX-GDDP longitude runs 0-360.
lon_query = mod(Lon,360);

vars = struct( ...
    'name',  {'tas','hurs','pr','rsds','sfcWind'}, ...
    'field', {'dT','dRH','fP','fSW','fWs'}, ...
    'mode',  {'additive','additive','ratio','ratio','ratio'});

nm = numel(models);
nv = numel(vars);
raw = nan(nm,nv,12);
for k = 1:nm
    for v = 1:nv
        base = monthly_climatology(nexgddp_root,models{k},'historical', ...
            vars(v).name,opt.BaselineYears,Lat,lon_query,opt.Verbose);
        futu = monthly_climatology(nexgddp_root,models{k},scenario, ...
            vars(v).name,opt.FutureYears,Lat,lon_query,opt.Verbose);
        switch vars(v).mode
            case 'additive'
                raw(k,v,:) = futu - base;
            case 'ratio'
                if any(base<=0)
                    error('Extract_CMIP6_Monthly_Deltas:NonPositiveBaseline', ...
                        ['Non-positive baseline climatology for %s in model %s; ' ...
                        'a ratio delta is undefined.'],vars(v).name,models{k});
                end
                raw(k,v,:) = futu ./ base;
        end
    end
end

if any(~isfinite(raw(:)))
    error('Extract_CMIP6_Monthly_Deltas:NonFiniteDelta', ...
        'Delta extraction produced NaN or infinite values.');
end

deltas = struct();
deltas.scenario = scenario;
deltas.models = models(:);
deltas.Lat = Lat;
deltas.Lon = Lon;
deltas.BaselineYears = opt.BaselineYears;
deltas.FutureYears = opt.FutureYears;
deltas.created = datestr(now,'yyyy-mm-ddTHH:MM:SS');

for v = 1:nv
    across = squeeze(raw(:,v,:));          % models x 12
    if nm==1
        across = reshape(across,1,12);
    end
    deltas.(vars(v).field) = mean(across,1);
    deltas.([vars(v).field '_members']) = across;
    deltas.([vars(v).field '_p10']) = prctile_local(across,10);
    deltas.([vars(v).field '_p90']) = prctile_local(across,90);
end

% Sanity bounds. These catch unit or window mistakes, not physical surprises.
if any(deltas.dT<-1) || any(deltas.dT>12)
    warning('Extract_CMIP6_Monthly_Deltas:ImplausibleWarming', ...
        'Monthly dT outside -1 to 12 K; check the window and scenario.');
end
if any(deltas.fP<0.3) || any(deltas.fP>3)
    warning('Extract_CMIP6_Monthly_Deltas:ImplausiblePrecipFactor', ...
        'Monthly fP outside 0.3 to 3.0; check the baseline window.');
end

if ~isempty(opt.PrecipQuantiles)
    deltas.PrecipQuantiles = opt.PrecipQuantiles(:)';
    deltas.fP_quantile = precip_quantile_ratios(nexgddp_root,models,scenario, ...
        opt.BaselineYears,opt.FutureYears,Lat,lon_query, ...
        opt.PrecipQuantiles,opt.Verbose);
end

save(output_file,'-struct','deltas','-v7');
fprintf('Wrote %s: %d models, scenario %s, %d-%d vs %d-%d.\n', ...
    output_file,nm,scenario,opt.BaselineYears(1),opt.BaselineYears(2), ...
    opt.FutureYears(1),opt.FutureYears(2));
end

% ---------------------------------------------------------------------

function clim = monthly_climatology(root,model,scenario,vname,years,Lat,lon_query,verbose)
% Mean annual cycle of one variable at the nearest grid cell.
%
% The historical baseline is identical across scenarios, so without a cache
% it would be re-read for each one. Keyed on everything that changes the
% result, so a scenario climatology is never wrongly reused.
persistent clim_cache
if isempty(clim_cache)
    clim_cache = containers.Map('KeyType','char','ValueType','any');
end
key = sprintf('%s|%s|%s|%d-%d|%.6f|%.6f', ...
    model,scenario,vname,years(1),years(2),Lat,lon_query);
if isKey(clim_cache,key)
    if verbose
        fprintf('  %s %s %s %d-%d (cached)\n', ...
            model,scenario,vname,years(1),years(2));
    end
    clim = clim_cache(key);
    return
end

totals = zeros(1,12);
counts = zeros(1,12);
for yr = years(1):years(2)
    f = locate_year_file(root,model,scenario,vname,yr);
    if isempty(f)
        error('Extract_CMIP6_Monthly_Deltas:MissingYear', ...
            'No %s file for %s %s %d under %s',vname,model,scenario,yr,root);
    end
    if verbose
        fprintf('  %s %s %s %d\n',model,scenario,vname,yr);
    end
    [ilat,ilon] = nearest_cell(f,Lat,lon_query);
    info = ncinfo(f,vname);
    nt = info.Size(3);
    series = squeeze(ncread(f,vname,[ilon ilat 1],[1 1 nt]));
    series = double(series(:));
    mo = month_index(f,yr,nt);
    for m = 1:12
        sel = mo==m;
        if any(sel)
            totals(m) = totals(m) + sum(series(sel));
            counts(m) = counts(m) + sum(sel);
        end
    end
end
if any(counts==0)
    error('Extract_CMIP6_Monthly_Deltas:EmptyMonth', ...
        'No days found for at least one month of %s %s.',model,vname);
end
clim = totals ./ counts;
clim_cache(key) = clim;
end

function f = locate_year_file(root,model,scenario,vname,yr)
pat = sprintf('%s_day_%s_%s_*%d*.nc',vname,model,scenario,yr);
if isfolder(root)
    hits = dir(fullfile(root,'**',pat));
    if isempty(hits)
        f = '';
    else
        f = fullfile(hits(1).folder,hits(1).name);
    end
else
    f = resolve_remote_file(root,model,scenario,vname,yr);
end
end

function f = resolve_remote_file(root,model,scenario,vname,yr)
% NEX-GDDP filenames carry a variant label, a grid label and a version suffix
% (_v1.1, _v2.0, or none) that vary by model and variable and cannot be
% derived from the model name. Read the THREDDS catalogue once per
% model/scenario/variable and cache the listing rather than guessing.
persistent cache
if isempty(cache)
    cache = containers.Map('KeyType','char','ValueType','any');
end
key = sprintf('%s|%s|%s',model,scenario,vname);
if ~isKey(cache,key)
    cache(key) = discover_remote_listing(root,model,scenario,vname);
end
listing = cache(key);

hit = listing.files(~cellfun(@isempty, ...
    regexp(listing.files,sprintf('_%d(_v[\\d.]+)?\\.nc$',yr),'once')));
if isempty(hit)
    error('Extract_CMIP6_Monthly_Deltas:RemoteYearNotFound', ...
        'No %s file for %s %s %d in the catalogue (%d files listed).', ...
        vname,model,scenario,yr,numel(listing.files));
end
f = sprintf('%s/%s/%s/%s/%s/%s', ...
    root,model,scenario,listing.variant,vname,hit{1});
end

function listing = discover_remote_listing(root,model,scenario,vname)
% THREDDS serves OPeNDAP under /dodsC/ and catalogues under /catalog/.
cat_root = strrep(root,'/dodsC/','/catalog/');
opts = weboptions('ContentType','text','Timeout',90);

fprintf('  resolving %s %s %s from the catalogue\n',model,scenario,vname);

variant_url = sprintf('%s/%s/%s/catalog.xml',cat_root,model,scenario);
try
    xml = webread(variant_url,opts);
catch err
    error('Extract_CMIP6_Monthly_Deltas:CatalogUnreachable', ...
        ['Could not read %s (%s). Check the model name and that %s is ' ...
        'published for it.'],variant_url,err.message,scenario);
end
variants = unique(regexp(xml,'r\d+i\d+p\d+f\d+','match'));
if isempty(variants)
    error('Extract_CMIP6_Monthly_Deltas:NoVariantFound', ...
        'No variant label found in %s',variant_url);
end
listing.variant = variants{1};

file_url = sprintf('%s/%s/%s/%s/%s/catalog.xml', ...
    cat_root,model,scenario,listing.variant,vname);
try
    xml = webread(file_url,opts);
catch err
    error('Extract_CMIP6_Monthly_Deltas:CatalogUnreachable', ...
        'Could not read %s (%s)',file_url,err.message);
end
files = unique(regexp(xml,sprintf('%s_day_[^"<>\\s]*?\\.nc',vname),'match'));
if isempty(files)
    error('Extract_CMIP6_Monthly_Deltas:NoFilesFound', ...
        'No %s NetCDF files listed in %s',vname,file_url);
end
listing.files = files(:);
fprintf('    variant %s, %d files listed\n', ...
    listing.variant,numel(listing.files));
end

function [ilat,ilon] = nearest_cell(f,Lat,lon_query)
lat_name = axis_name(f,{'lat','latitude'});
lon_name = axis_name(f,{'lon','longitude'});
lats = double(ncread(f,lat_name));
lons = double(ncread(f,lon_name));
[~,ilat] = min(abs(lats-Lat));
dlon = abs(mod(lons-lon_query+180,360)-180);
[~,ilon] = min(dlon);
end

function name = axis_name(f,candidates)
available = {ncinfo(f).Variables.Name};
for k = 1:numel(candidates)
    idx = find(strcmpi(available,candidates{k}),1);
    if ~isempty(idx)
        name = available{idx};
        return
    end
end
error('Extract_CMIP6_Monthly_Deltas:MissingAxis', ...
    'No coordinate axis found. Expected one of: %s',strjoin(candidates,', '));
end

function mo = month_index(f,yr,nt)
% CMIP6 models use 365_day, 360_day or standard calendars. Deriving the
% month from the day index per calendar avoids datenum errors that would
% silently shift the annual cycle.
cal = 'standard';
try
    tname = axis_name(f,{'time'});
    cal = lower(ncreadatt(f,tname,'calendar'));
catch
end
switch cal
    case {'360_day','360day'}
        lens = repmat(30,1,12);
    case {'365_day','365day','noleap'}
        lens = [31 28 31 30 31 30 31 31 30 31 30 31];
    case {'366_day','366day','all_leap'}
        lens = [31 29 31 30 31 30 31 31 30 31 30 31];
    otherwise
        leap = (mod(yr,4)==0 && mod(yr,100)~=0) || mod(yr,400)==0;
        lens = [31 28+leap 31 30 31 30 31 31 30 31 30 31];
end
if sum(lens)~=nt
    error('Extract_CMIP6_Monthly_Deltas:CalendarMismatch', ...
        'Year %d has %d records but calendar "%s" implies %d.', ...
        yr,nt,cal,sum(lens));
end
mo = zeros(nt,1);
p = 0;
for m = 1:12
    mo(p+(1:lens(m))) = m;
    p = p + lens(m);
end
end

function ratios = precip_quantile_ratios(root,models,scenario,base_years, ...
    fut_years,Lat,lon_query,probs,verbose)
% Wet-day quantile ratios by month, averaged across models. Used by the
% QDM precipitation option in Apply_CMIP6_Deltas_To_ERA5_Land.
probs = probs(:)';
wet_threshold = 1/86400;   % 1 mm/day expressed as kg m-2 s-1
acc = nan(numel(models),12,numel(probs));
for k = 1:numel(models)
    if verbose
        fprintf('  quantiles: %s\n',models{k});
    end
    base = collect_daily(root,models{k},'historical','pr',base_years,Lat,lon_query);
    futu = collect_daily(root,models{k},scenario,'pr',fut_years,Lat,lon_query);
    for m = 1:12
        b = base.values(base.month==m);
        f = futu.values(futu.month==m);
        b = b(b>wet_threshold);
        f = f(f>wet_threshold);
        if numel(b)<30 || numel(f)<30
            error('Extract_CMIP6_Monthly_Deltas:TooFewWetDays', ...
                'Month %d has too few wet days for stable quantiles.',m);
        end
        qb = quantile_local(b,probs);
        qf = quantile_local(f,probs);
        acc(k,m,:) = qf ./ qb;
    end
end
ratios = squeeze(mean(acc,1));
if isscalar(models)
    ratios = reshape(ratios,12,numel(probs));
end
end

function out = collect_daily(root,model,scenario,vname,years,Lat,lon_query)
% Cached on the same basis as monthly_climatology: the historical wet-day
% sample is shared by every scenario.
persistent daily_cache
if isempty(daily_cache)
    daily_cache = containers.Map('KeyType','char','ValueType','any');
end
key = sprintf('%s|%s|%s|%d-%d|%.6f|%.6f', ...
    model,scenario,vname,years(1),years(2),Lat,lon_query);
if isKey(daily_cache,key)
    out = daily_cache(key);
    return
end

vals = [];
mos = [];
for yr = years(1):years(2)
    f = locate_year_file(root,model,scenario,vname,yr);
    if isempty(f)
        error('Extract_CMIP6_Monthly_Deltas:MissingYear', ...
            'No %s file for %s %s %d.',vname,model,scenario,yr);
    end
    [ilat,ilon] = nearest_cell(f,Lat,lon_query);
    nt = ncinfo(f,vname).Size(3);
    s = double(squeeze(ncread(f,vname,[ilon ilat 1],[1 1 nt])));
    vals = [vals; s(:)];
    mos  = [mos; month_index(f,yr,nt)];
end
out.values = vals;
out.month = mos;
daily_cache(key) = out;
end

function q = quantile_local(x,probs)
% Linear interpolation of order statistics; avoids a Statistics Toolbox call.
x = sort(x(:));
n = numel(x);
pos = probs(:)*n + 0.5;
lo = floor(pos);
hi = ceil(pos);
lo = min(max(lo,1),n);
hi = min(max(hi,1),n);
w = pos - floor(pos);
q = (1-w).*x(lo) + w.*x(hi);
q = reshape(q,1,[]);
end

function p = prctile_local(a,pct)
p = zeros(1,size(a,2));
for j = 1:size(a,2)
    p(j) = quantile_local(a(:,j),pct/100);
end
end

function opt = parse_options(args)
opt.BaselineYears = [1985 2014];
opt.FutureYears = [2036 2065];
opt.PrecipQuantiles = [];
opt.Verbose = true;
if mod(numel(args),2)~=0
    error('Extract_CMIP6_Monthly_Deltas:BadOptions', ...
        'Options must be supplied as name-value pairs.');
end
for k = 1:2:numel(args)
    name = validatestring(args{k}, ...
        {'BaselineYears','FutureYears','PrecipQuantiles','Verbose'});
    opt.(name) = args{k+1};
end
validateattributes(opt.BaselineYears,{'numeric'},{'vector','numel',2,'integer','increasing'});
validateattributes(opt.FutureYears,{'numeric'},{'vector','numel',2,'integer','increasing'});
if ~isempty(opt.PrecipQuantiles)
    validateattributes(opt.PrecipQuantiles,{'numeric'}, ...
        {'vector','>',0,'<',1,'increasing'});
end
end
