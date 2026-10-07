function gamma = Variogram_Semivariance(variogram_model,model_parameters,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Semivariance of the PyKrige variogram models, in their own parameterisation.
%%
%% The range parameter is the effective range, not the exponential decay
%% length: exponential decays over range/3 and gaussian over range*4/7,
%% exactly as pykrige.variogram_models does. Keeping that convention is what
%% lets a MATLAB surface reproduce the notebook value for value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~(isnumeric(d) && isreal(d)) || any(d(:) < 0)
    error('Variogram_Semivariance:BadDistance', ...
        'd must hold real non-negative separation distances.');
end
model_parameters = double(model_parameters(:)).';

switch lower(char(variogram_model))
    case 'linear'
        expect_parameters(model_parameters,2,'linear');
        slope = model_parameters(1);
        nugget = model_parameters(2);
        gamma = slope*d+nugget;
    case 'power'
        expect_parameters(model_parameters,3,'power');
        scale = model_parameters(1);
        exponent = model_parameters(2);
        nugget = model_parameters(3);
        gamma = scale*d.^exponent+nugget;
    case 'gaussian'
        expect_parameters(model_parameters,3,'gaussian');
        [psill,range_m,nugget] = deal(model_parameters(1), ...
            model_parameters(2),model_parameters(3));
        gamma = psill*(1-exp(-(d.^2)/(range_m*4/7)^2))+nugget;
    case 'exponential'
        expect_parameters(model_parameters,3,'exponential');
        [psill,range_m,nugget] = deal(model_parameters(1), ...
            model_parameters(2),model_parameters(3));
        gamma = psill*(1-exp(-d/(range_m/3)))+nugget;
    case 'spherical'
        expect_parameters(model_parameters,3,'spherical');
        [psill,range_m,nugget] = deal(model_parameters(1), ...
            model_parameters(2),model_parameters(3));
        gamma = psill+nugget+zeros(size(d));
        within = d <= range_m;
        gamma(within) = psill*(3*d(within)/(2*range_m)- ...
            d(within).^3/(2*range_m^3))+nugget;
    otherwise
        error('Variogram_Semivariance:UnsupportedModel', ...
            ['variogram_model must be linear, power, gaussian, ' ...
             'exponential or spherical.']);
end

end

function expect_parameters(model_parameters,expected_count,model_name)
if numel(model_parameters) ~= expected_count
    error('Variogram_Semivariance:BadParameterCount', ...
        'The %s model takes %d parameters, not %d.', ...
        model_name,expected_count,numel(model_parameters));
end
end
