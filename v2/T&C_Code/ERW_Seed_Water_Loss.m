function Q = ERW_Seed_Water_Loss(ERW_TC_const_input,TT_m_input)
Q = NaN;
if isstruct(ERW_TC_const_input) && ...
        isfield(ERW_TC_const_input,'seed_water_loss_l_m2_d')
    q_input = ERW_TC_const_input.seed_water_loss_l_m2_d;
    if isscalar(q_input) && isfinite(q_input) && q_input > 0
        Q = q_input;
        return
    end
end

if isstruct(TT_m_input) && isfield(TT_m_input,'Tw') && ...
        isfield(TT_m_input,'Lw')
    Tw = TT_m_input.Tw(:);
    Lw = TT_m_input.Lw(:);
    if ~isempty(Tw) && isequal(size(Tw),size(Lw)) && ...
            all(isfinite(Tw)) && all(isfinite(Lw))
        q_forcing = mean(Tw+Lw)*1000;
        if isfinite(q_forcing) && q_forcing > 0
            Q = q_forcing;
        end
    end
end
end
