function test_erw_particle_distribution()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A cumulative particle-size distribution must retain its class masses.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

d_um = [1 250 500 1000 2000 3000 4000];
cdf_percent = [0 30 55 78 92 98 100];
expected = [0 0.30 0.25 0.23 0.14 0.06 0.02];

mass_fraction = ERW_Particle_Mass_Fractions(d_um,cdf_percent);

assert(max(abs(mass_fraction-expected)) < 1e-14, ...
    'CDF conversion changed the measured particle-class masses.');
assert(abs(sum(mass_fraction)-1) < 1e-14, ...
    'Particle-class mass fractions do not sum to one.');
assert_error(@() ERW_Particle_Mass_Fractions( ...
    [0 250 500 1000 2000 3000 4000],cdf_percent), ...
    'ERW:invalidParticleDistribution');

fprintf('PASS: cumulative particle sizes retain measured class masses.\n');

end

function assert_error(fn,id)
try
    fn();
catch ME
    assert(strcmp(ME.identifier,id), ...
        'Expected %s, received %s.',id,ME.identifier);
    return
end
error('Expected error %s was not raised.',id);
end
