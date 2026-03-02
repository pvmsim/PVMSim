function k = metrics_kpis(runout)
%METRICS_KPIS  Unified KPIs for GUI and CLI (explicit I–V/P–V errors).
%
% Definitions (explicit, pointwise unless stated):
%   eI = I_meas - I_calc         [A]
%   eP = P_meas - P_calc         [W]
%   RMSE = sqrt(mean(e.^2))      (no weighting)
%   MAE  = mean(abs(e))          (no weighting)
%   IAE  = trapz(V_sorted, abs(e_sorted))  (integral over V)
%   RE_Pmpp = 100*|Pmpp_calc - Pmpp_meas|/|Pmpp_meas|
%
% Requires runout fields:
%   V_meas, I_meas, I_calc
% Optional:
%   P_meas, P_calc (if not present, computed from V and I)

    % ---- required ----
    mustHave = {'V_meas','I_meas','I_calc'};
    for i = 1:numel(mustHave)
        if ~isfield(runout, mustHave{i}) || isempty(runout.(mustHave{i}))
            error('metrics_kpis:MissingField', 'runout.%s is required.', mustHave{i});
        end
    end

    V  = runout.V_meas(:);
    Im = runout.I_meas(:);
    Ic = runout.I_calc(:);

    % ---- sanitize finite points ----
    ok = isfinite(V) & isfinite(Im) & isfinite(Ic);
    if nnz(ok) < 2
        error('metrics_kpis:NotEnoughPoints', 'Need at least 2 finite points for KPIs.');
    end
    V  = V(ok);
    Im = Im(ok);
    Ic = Ic(ok);

    % ---- compute P if needed ----
    if isfield(runout,'P_meas') && ~isempty(runout.P_meas)
        Pm = runout.P_meas(:); Pm = Pm(ok);
    else
        Pm = V .* Im;
    end

    if isfield(runout,'P_calc') && ~isempty(runout.P_calc)
        Pc = runout.P_calc(:); Pc = Pc(ok);
    else
        Pc = V .* Ic;
    end

    % ---- errors ----
    eI = Im - Ic;     % [A]
    eP = Pm - Pc;     % [W]

    % ---- pointwise RMSE/MAE ----
    k.RMSE_I = sqrt(mean(eI.^2));
    k.MAE_I  = mean(abs(eI));

    k.RMSE_P = sqrt(mean(eP.^2));
    k.MAE_P  = mean(abs(eP));

    % ---- IAE as integral over V (trapz) ----
    [Vs, idx] = sort(V, 'ascend');
    eIs = eI(idx);
    ePs = eP(idx);

    k.IAE_I = trapz(Vs, abs(eIs));  % [A·V]
    k.IAE_P = trapz(Vs, abs(ePs));  % [W·V]

    % ---- Relative error at Pmpp ----
    Pmpp_me = max(Pm);
    Pmpp_ca = max(Pc);

    if isfinite(Pmpp_me) && abs(Pmpp_me) > 1e-12
        k.RE_Pmpp = 100 * abs(Pmpp_ca - Pmpp_me) / abs(Pmpp_me); % [%]
    else
        k.RE_Pmpp = NaN;
    end

    % ---- convenience (optional) ----
    k.Pmpp_meas = Pmpp_me;
    k.Pmpp_calc = Pmpp_ca;
    k.N_points  = numel(V);
end