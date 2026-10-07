function out = Fit_Empirical_Variogram_Range(x_m, y_m, values, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit_Empirical_Variogram_Range
%%
%% Build the isotropic empirical (binned) semivariogram of a scattered
%% field and least-squares-fit an exponential model to recover its
%% effective correlation range.  This is the measured length scale that
%% Run_Topographic_Grid_Convergence uses to size the process grid -- the
%% fixed 80 m "notebook recipe" range in Fit_Ordinary_Kriging is a default,
%% not a site measurement.
%%
%% INPUTS
%%   x_m, y_m   [N x 1] projected sample coordinates [m]
%%   values     [N x 1] the sampled attribute (any units; detrend first for
%%              a surface with a regional gradient, e.g. elevation)
%%
%% Name/value options
%%   'MaxLagM'   longest separation binned [m].  Default: 0.4 * the
%%               sample-cloud diameter.
%%   'NBins'     number of lag bins (default 15)
%%   'MinPairs'  bins with fewer pairs than this are dropped (default 20)
%%   'Model'     'exponential' (default) or 'gaussian' or 'spherical'
%%               (uses Variogram_Semivariance for the model form)
%%
%% OUTPUT struct
%%   .range_m       fitted effective range [m]  (exponential: gamma reaches
%%                  ~95% sill at range_m; the decay length is range_m/3, the
%%                  pykrige convention Variogram_Semivariance keeps)
%%   .sill .nugget  fitted partial sill + nugget (value units^2)
%%   .lag_m .gamma .n_pairs   the empirical variogram that was fitted
%%   .rmse          fit residual RMSE
%%   .model         the model name used
%%   .n_samples     number of finite samples used
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
addpath(here);

opt = struct('MaxLagM',[],'NBins',15,'MinPairs',20,'Model','exponential');
for k = 1:2:numel(varargin); opt.(varargin{k}) = varargin{k+1}; end

x = double(x_m(:)); y = double(y_m(:)); v = double(values(:));
ok = isfinite(x) & isfinite(y) & isfinite(v);
x = x(ok); y = y(ok); v = v(ok);
n = numel(v);
if n < 10
    error('Fit_Empirical_Variogram_Range:TooFewSamples', ...
        'Need at least 10 located finite samples, got %d.',n);
end

%% All unique pairs.
[ii,jj] = find(triu(true(n),1));
d  = hypot(x(ii)-x(jj), y(ii)-y(jj));
sv = 0.5*(v(ii)-v(jj)).^2;

diameter = max(hypot(x-mean(x),y-mean(y)))*2;
max_lag = opt.MaxLagM;
if isempty(max_lag); max_lag = 0.4*diameter; end
max_lag = min(max_lag, max(d));

edges = linspace(0, max_lag, opt.NBins+1);
lag = zeros(opt.NBins,1); gamma = zeros(opt.NBins,1); npair = zeros(opt.NBins,1);
for b = 1:opt.NBins
    in = d > edges(b) & d <= edges(b+1);
    npair(b) = nnz(in);
    if npair(b) > 0
        lag(b)   = mean(d(in));
        gamma(b) = mean(sv(in));
    end
end
keep = npair >= opt.MinPairs & lag > 0;
lag = lag(keep); gamma = gamma(keep); npair = npair(keep);
if numel(lag) < 4
    error('Fit_Empirical_Variogram_Range:SparseVariogram', ...
        ['Only %d usable lag bins (>= %d pairs). Widen MaxLagM, lower ' ...
         'MinPairs, or accept that the sampling cannot resolve a range.'], ...
        numel(lag),opt.MinPairs);
end

%% Least-squares fit of [partial_sill, range, nugget], pair-count weighted.
model = lower(char(opt.Model));
var_v = var(v,1);
p0 = [max(var_v,eps); max(lag)/2; 0.05*var_v];
lb = [0; edges(2); 0];
ub = [10*var_v + eps; 3*max_lag; var_v];
w  = sqrt(npair);
resid = @(p) w.*(Variogram_Semivariance(model,[p(1) p(2) p(3)],lag) - gamma);

if exist('lsqnonlin','file') == 2
    o = optimoptions('lsqnonlin','Display','off');
    p = lsqnonlin(resid,p0,lb,ub,o);
else
    %% Toolbox-free fallback: fine grid search over the range, with the
    %% sill and nugget solved in closed form (they are linear) at each
    %% trial range.
    best = inf; p = p0;
    for r = linspace(edges(2), 3*max_lag, 400)
        [ps,ng,e] = linear_sill_nugget(model,r,lag,gamma,w,var_v);
        if e < best; best = e; p = [ps; r; ng]; end
    end
end

fit = Variogram_Semivariance(model,[p(1) p(2) p(3)],lag);
out = struct();
out.range_m   = p(2);
out.sill      = p(1);
out.nugget    = p(3);
out.lag_m     = lag;
out.gamma     = gamma;
out.n_pairs   = npair;
out.model     = model;
out.n_samples = n;
out.rmse      = sqrt(mean((fit-gamma).^2));
out.decay_length_m = p(2)/3;    % exponential e-folding (pykrige convention)
%% A fitted range below the shortest sampled lag, or a nugget that is most
%% of the sill, means the sampling cannot resolve any spatial structure --
%% the attribute imposes no grid-resolution constraint.
out.shortest_lag_m = min(lag);
out.nugget_fraction = p(3)/max(p(1)+p(3),eps);
out.resolved = out.range_m > out.shortest_lag_m && out.nugget_fraction < 0.75;
if ~out.resolved
    warning('Fit_Empirical_Variogram_Range:Unresolved', ...
        ['Fitted range %.1f m is below the shortest sampled lag %.1f m ' ...
         '(nugget %.0f%% of sill): the sampling cannot resolve spatial ' ...
         'structure in this attribute.'], out.range_m, out.shortest_lag_m, ...
         100*out.nugget_fraction);
end
end

%% =======================================================================
function [psill,nugget,err] = linear_sill_nugget(model,range_m,lag,gamma,w,var_v)
%% For a fixed range, the exponential/gaussian/spherical semivariance is
%% linear in (partial sill, nugget): gamma = psill*shape(lag;range) + nugget.
shape = Variogram_Semivariance(model,[1 range_m 0],lag);   % psill=1, nugget=0
A = [w.*shape, w];
b = w.*gamma;
c = lsqnonneg_local(A,b);
psill = c(1); nugget = min(c(2), var_v);
err = norm(A*[psill;nugget] - b);
end

%% =======================================================================
function x = lsqnonneg_local(A,b)
if exist('lsqnonneg','file') == 2
    x = lsqnonneg(A,b);
else
    x = max(A\b, 0);
end
end

