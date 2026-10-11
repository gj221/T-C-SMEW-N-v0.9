function fraction = Plant_Reserve_Fraction(reserve,capacity,element_name)
% Reserve fill fraction, reproducing the legacy TC-SMEW clamped ratio.
%
% Legacy Plant_Exports computed this inline and clamped the *result*, never
% the inputs:
%
%   fll = Nreserve./Nsto_tm1; fll(fll>1)=1; fll(fll<0)=0;
%
% A negative flexible storage therefore produced a negative ratio that
% clamped to zero, and the run continued. Validating the inputs instead and
% erroring on a negative capacity was a behavioural change, not a port: it
% turns a benign excursion at an exhausted reserve into a hard stop. The
% clamped-ratio behaviour is restored here.
%
% Two departures from legacy are kept deliberately. A non-finite input still
% raises, because legacy propagated NaN silently and that is a genuine defect
% rather than a tolerated state. And the 0/0 case, which legacy left as NaN,
% returns the saturated limit so the zero-capacity boundary stays defined.

if ~isfinite(reserve) || ~isfinite(capacity)
    error('Plant_Reserve_Fraction:InvalidState', ...
        '%s reserve and storage capacity must be finite.',element_name);
end

tol = 100*eps(max([1,reserve,capacity]));

%% Zero-capacity boundary. Legacy divided by zero here and clamped, giving
%% the saturated limit for a positive reserve. Only a capacity that is zero
%% to within tolerance takes this branch: a capacity that is genuinely
%% negative falls through to the ratio, where legacy clamped it to zero.
if abs(capacity) <= tol
    fraction = double(reserve > tol);
    return
end

fraction = min(max(reserve/capacity,0),1);
end
