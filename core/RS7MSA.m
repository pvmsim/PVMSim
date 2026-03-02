function runout = RS7MSA(panel, V_meas, I_meas, meta, logFcn)
%RS7MSA  DDM + RS6 + multi-start SA in one call (source of truth).

%% 0) Inputs
if nargin < 4 || isempty(meta)
    meta = struct('ID','','G',NaN,'T',NaN,'Date','');
end
if nargin < 5
    logFcn = [];
end

V_meas = V_meas(:);
I_meas = I_meas(:);
curves = {V_meas, I_meas};


% --- 0b) Controls from meta (Seed, MaxIter) ---
DEFAULT_SEED    = 42;
DEFAULT_MAXITER = 80000;

% Seed (defensive normalization; utils_rng will also enforce)
if ~isstruct(meta) || isempty(meta)
    meta = struct('ID','','G',NaN,'T',NaN,'Date','');
end

if isfield(meta,'Seed') && isfinite(meta.Seed)
    meta.Seed = double(round(meta.Seed));
else
    meta.Seed = DEFAULT_SEED;
end
meta.Seed = min(max(meta.Seed, 1), 4294967295);

% MaxIter (drive both Phase 1 and Phase 2 lengths)
if isfield(meta,'MaxIter') && isfinite(meta.MaxIter)
    meta.MaxIter = double(round(meta.MaxIter));
else
    meta.MaxIter = DEFAULT_MAXITER;
end
meta.MaxIter = min(max(meta.MaxIter, 1000), 500000);

localLog('INFO', sprintf('Meta controls | Seed=%u | MaxIter=%d', uint32(meta.Seed), meta.MaxIter));



localLog('INFO', sprintf( ...
    'Starting RS7MSA for curve "%s" (%d points, G=%.1f W/m², T=%.1f °C).', ...
    meta.ID, numel(V_meas), meta.G, meta.T));

%% 3) Constants
q = 1.602176634e-19;
k = 1.380649e-23;
epsv = 1e-12;

T_K = panel.Tref_C + 273.15;
Vt  = panel.Ns * (k * T_K / q);

fprintf('T [K] = %.2f\n', T_K);
fprintf('Vt (módulo) = %.12f\n', Vt);

%% 4) Bounds/steps (C)
C = struct();

C.IPH_LB = 0.9 * panel.Isc;
C.IPH_UB = 1.1 * panel.Isc;

C.IO1_LB = 1e-12;  C.IO1_UB = 50e-6;
C.IO2_LB = 1e-12;  C.IO2_UB = 50e-6;

C.RS_LB  = 0.0;    C.RS_UB  = 0.8;
C.RSH_LB = 5.0;    C.RSH_UB = 2e3;

C.N1_LB = 1.0;     C.N1_UB = 2.0;
C.N2_LB = 1.0;     C.N2_UB = 2.0;

C.D_RS_COARSE      = 0.10;
C.D_RSH_COARSE     = 10.0;
C.D_N_COARSE       = 0.10;
C.IO_FACTOR_COARSE = 2.0;

C.D_RS_FINE        = 0.01;
C.D_RSH_FINE       = 2.0;
C.D_N_FINE         = 0.02;
C.IO_FACTOR_FINE   = 1.25;

C.TOPK = 20;
C.EPS  = epsv;

%% 5) RNG
[seed, rng_state0] = utils_rng(meta, DEFAULT_SEED);
localLog('INFO', sprintf('RNG initialized | Seed=%u | Generator=twister', seed));

%% 6) Phase 1: RS6 + diversity
paramsF1 = struct();
paramsF1.iters         = meta.MaxIter;
paramsF1.p_global      = 0.75;
paramsF1.p_kick        = 0.15;
paramsF1.topk          = C.TOPK;
paramsF1.k_per_bucket  = 2;
paramsF1.verbose_every = 5000;

[best_theta7_f1, best_rmse_f1, topk_phase1] = random_search_phase1(curves, Vt, C, paramsF1);

fprintf('\nRMSE residual (fase 1): %.15f\n', best_rmse_f1);
fprintf('theta fase 1 [Iph, Io1, Io2, Rs, Rsh, n1, n2]:\n');
disp(best_theta7_f1);

% Sort top-k by RMSE before Phase 2
if ~isempty(topk_phase1)
    rm = arrayfun(@(s) s.rmse, topk_phase1(:));
    ok = isfinite(rm);
    nBad = sum(~ok);

    if ~any(ok)
        error('Phase 1 produced only non-finite RMSE seeds. Aborting Phase 2.');
    end

    topk_phase1 = topk_phase1(ok);
    rm = rm(ok);

    [rm, idx] = sort(rm, 'ascend');
    topk_phase1 = topk_phase1(idx);

    localLog('INFO', sprintf(['Top-k seeds sorted for Phase 2. ' ...
        'Best RMSE=%.6g, Worst RMSE=%.6g, k=%d, dropped=%d'], ...
        rm(1), rm(end), numel(rm), nBad));
end

%% 7) Phase 2: multi-start SA
paramsF2 = struct();
paramsF2.iters_per_start = meta.MaxIter;
paramsF2.starts          = 10;
paramsF2.T0              = 0.02;
paramsF2.Tend            = 1e-5;
paramsF2.p_restart_best  = 0.02;
paramsF2.p_global_jump   = 0.01;
paramsF2.verbose_every   = 5000;

[best_theta7_f2, best_rmse_f2, hist2] = phase2_simulated_annealing(topk_phase1, curves, Vt, C, paramsF2);

fprintf('\nRMSE residual (fase 2): %.15f\n', best_rmse_f2);
fprintf('theta fase 2 [Iph, Io1, Io2, Rs, Rsh, n1, n2]:\n');
format shortG
disp(best_theta7_f2)

%% 8) I–V calculated (Newton)
I_calc = RS7MSA_model_current(best_theta7_f2, V_meas, Vt, C);

rmse_I = sqrt(mean((I_meas - I_calc).^2));
mae_I  = mean(abs(I_meas - I_calc));

fprintf('\nRMSE(I_meas vs I_calc): %.6f A\n', rmse_I);
fprintf('MAE(I_meas vs I_calc):  %.6f A\n', mae_I);

%% 9) Pack runout
rng_stateF = rng;

P_meas = V_meas .* I_meas;
P_calc = V_meas .* I_calc;

runout = struct();
runout.panel       = panel;
runout.V_meas      = V_meas;
runout.I_meas      = I_meas;
runout.I_calc      = I_calc;
runout.P_meas      = P_meas;
runout.P_calc      = P_calc;
runout.rmse_I      = rmse_I;
runout.mae_I       = mae_I;

runout.rmse_obj    = best_rmse_f2;
runout.rmse_IV     = rmse_I;

runout.meta        = meta;
runout.theta7      = best_theta7_f2;
runout.hist2       = hist2;
runout.rmse_phase2 = best_rmse_f2;
runout.res_I       = I_meas - I_calc;

runout.seed     = seed;
runout.rng0     = rng_state0;
runout.rngFinal = rng_stateF;

localLog('INFO', sprintf('Finished RS7MSA | %s | RMSE_I = %.4f A', ...
    string(panel.modelo), rmse_I));

%% 10) Summary
fprintf('\n=== RESUMEN FINAL ===\n');
fprintf('Modelo: %s\n', panel.modelo);
fprintf('Fase 1 RMSE residual: %.8f\n', best_rmse_f1);
fprintf('Fase 2 RMSE residual: %.8f\n', best_rmse_f2);
fprintf('RMSE curva I-V (corriente): %.8f A\n', rmse_I);
fprintf('Theta final [Iph, Io1, Io2, Rs, Rsh, n1, n2]:\n');
format shortG
disp(best_theta7_f2)

    function localLog(level,msg)
        if ~isempty(logFcn)
            try
                logFcn(level,msg);
            catch
            end
        end
    end

end

% ---------------------------------------------------------------------
function [best_theta7, best_rmse, topk_list] = random_search_phase1(curves, Vt, C, P)

    buckets = containers.Map('KeyType','char','ValueType','any');

    best_rmse   = inf;
    best_theta7 = [];
    best_theta6 = [];
    best_iph    = NaN;

    for t = 1:P.iters

        if isempty(best_theta6) || (rand < P.p_global)
            cand6 = sample_global_theta6(C);
        else
            cand6 = mutate_RS7MSA_phase1(best_theta6, C, 0.50);
        end

        cand6 = clamp_theta6(cand6, C);

        if rand < P.p_kick
            Rs  = cand6(1);
            Rsh = cand6(2);
            Io1 = cand6(3);
            Io2 = cand6(4);
            n1  = cand6(5);
            n2  = cand6(6);

            rsh_sigma = 0.15 * max(C.RSH_UB - C.RSH_LB, 1.0);
            Rsh = clamp_val(Rsh + randn * rsh_sigma, C.RSH_LB, C.RSH_UB);

            logIo2 = log10(max(Io2, C.IO2_LB)) + randn * 2.0;
            logIo2 = min(max(logIo2, log10(C.IO2_LB)), log10(C.IO2_UB));
            Io2 = 10^logIo2;

            n2_sigma = 0.20 * max(C.N2_UB - C.N2_LB, 1e-6);
            n2 = clamp_val(n2 + randn * n2_sigma, C.N2_LB, C.N2_UB);

            cand6 = clamp_theta6([Rs, Rsh, Io1, Io2, n1, n2], C);
        end

        [rmse_c, iph_c, theta7_c] = objective_autoph(cand6, curves, Vt, C);

        buckets = push_topk_stratified(buckets, rmse_c, iph_c, cand6, C, P.k_per_bucket);

        if isfinite(rmse_c) && (rmse_c < best_rmse)
            best_rmse   = rmse_c;
            best_theta6 = cand6;
            best_theta7 = theta7_c;
            best_iph    = iph_c;
        end

        if P.verbose_every > 0 && mod(t, P.verbose_every) == 0 && ~isempty(best_theta6)
            Rs  = best_theta6(1);
            Rsh = best_theta6(2);
            Io1 = best_theta6(3);
            Io2 = best_theta6(4);
            n1  = best_theta6(5);
            n2  = best_theta6(6);
            fprintf(['[RS6-div] iter=%d best=%.6f Iph*=%.4f Rs=%.4f Rsh=%.2f ' ...
                     'Io1=%.3e Io2=%.3e n1=%.3f n2=%.3f buckets=%d\n'], ...
                     t, best_rmse, best_iph, Rs, Rsh, Io1, Io2, n1, n2, buckets.Count);
        end
    end

    all_items = {};
    keys_ = buckets.keys;
    for i = 1:numel(keys_)
        lst = buckets(keys_{i});
        for j = 1:numel(lst)
            all_items{end+1,1} = lst{j};
        end
    end

    if isempty(all_items)
        topk_list = struct('rmse',{},'iph',{},'theta6',{});
        return;
    end

    rmse_vec = cellfun(@(s) s.rmse, all_items);
    [~, idx] = sort(rmse_vec, 'ascend');
    all_items = all_items(idx);

    nkeep = min(P.topk, numel(all_items));
    topk_list = [all_items{1:nkeep}];
end

% ---------------------------------------------------------------------
function [global_best_theta7, global_best_rmse, history] = phase2_simulated_annealing(topk_phase1, curves, Vt, C, P)

    if isempty(topk_phase1)
        error('topk_phase1 está vacío.');
    end

    n_starts = min(P.starts, numel(topk_phase1));
    global_best_rmse = inf;
    global_best_theta7 = [];
    history = zeros(n_starts * P.iters_per_start, 5);
    hidx = 0;

    if P.iters_per_start <= 1
        alpha = 1.0;
    else
        alpha = (P.Tend / P.T0)^(1 / (P.iters_per_start - 1));
    end

    for s = 1:n_starts
        seed6 = topk_phase1(s).theta6;
        curr6 = jitter_seed(seed6, C);

        [curr_rmse, curr_iph, curr7] = objective_autoph(curr6, curves, Vt, C);
        if ~isfinite(curr_rmse), curr_rmse = 1e9; end

        best_rmse = curr_rmse;
        best_iph  = curr_iph;
        best6     = curr6;
        best7     = curr7;

        T = P.T0;

        for t = 1:P.iters_per_start

            if rand < P.p_restart_best
                curr6 = best6;
                curr_rmse = best_rmse;
                curr_iph  = best_iph;
                curr7     = best7;
            end

            if rand < P.p_global_jump
                cand6 = sample_global_theta6(C);
            else
                cand6 = mutate_RS7MSA_phase2(curr6, C, 0.12);
            end

            [cand_rmse, cand_iph, cand7] = objective_autoph(cand6, curves, Vt, C);
            if ~isfinite(cand_rmse), cand_rmse = 1e9; end

            delta = cand_rmse - curr_rmse;

            if delta <= 0
                accept = true;
            else
                x = -delta / max(T, 1e-12);
                x = min(max(x, -700), 700);
                accept = (rand < exp(x));
            end

            if accept
                curr6 = cand6;
                curr_rmse = cand_rmse;
                curr_iph  = cand_iph;
                curr7     = cand7;

                if curr_rmse < best_rmse
                    best_rmse = curr_rmse;
                    best_iph  = curr_iph;
                    best6     = curr6;
                    best7     = curr7;
                end
            end

            if P.verbose_every > 0 && mod(t, P.verbose_every) == 0
                Rs  = best6(1);
                Rsh = best6(2);
                Io1 = best6(3);
                Io2 = best6(4);
                n1  = best6(5);
                n2  = best6(6);
                fprintf(['[SA-F2] start=%d/%d iter=%d best=%.6f curr=%.6f T=%.2e Iph*=%.4f ' ...
                         'Rs=%.4f Rsh=%.2f Io1=%.3e Io2=%.3e n1=%.3f n2=%.3f\n'], ...
                         s, n_starts, t, best_rmse, curr_rmse, T, best_iph, ...
                         Rs, Rsh, Io1, Io2, n1, n2);
            end

            hidx = hidx + 1;
            history(hidx,:) = [s-1, t, best_rmse, curr_rmse, T];
            T = T * alpha;
        end

        if best_rmse < global_best_rmse
            global_best_rmse = best_rmse;
            global_best_theta7 = best7;
        end
    end

    history = history(1:hidx,:);
end

% ---------------------------------------------------------------------
function theta6 = sample_global_theta6(C)
    Rs  = C.RS_LB  + rand * (C.RS_UB  - C.RS_LB);
    Rsh = C.RSH_LB + rand * (C.RSH_UB - C.RSH_LB);
    n1  = C.N1_LB  + rand * (C.N1_UB  - C.N1_LB);
    n2  = C.N2_LB  + rand * (C.N2_UB  - C.N2_LB);

    logIo1 = log10(C.IO1_LB) + rand * (log10(C.IO1_UB) - log10(C.IO1_LB));
    logIo2 = log10(C.IO2_LB) + rand * (log10(C.IO2_UB) - log10(C.IO2_LB));
    Io1 = 10^logIo1;
    Io2 = 10^logIo2;

    theta6 = [Rs, Rsh, Io1, Io2, n1, n2];
end

% ---------------------------------------------------------------------
function theta6_out = jitter_seed(theta6_in, C)
    theta6_out = theta6_in;
    for i = 1:4
        theta6_out = mutate_RS7MSA_phase1(theta6_out, C, 0.35);
    end
end

% ---------------------------------------------------------------------
function theta6c = clamp_theta6(theta6, C)

    Rs  = theta6(1);
    Rsh = theta6(2);
    Io1 = theta6(3);
    Io2 = theta6(4);
    n1  = theta6(5);
    n2  = theta6(6);

    if ~isfinite(Rs),  Rs  = C.RS_LB;  end
    if ~isfinite(Rsh), Rsh = C.RSH_LB; end
    if ~isfinite(Io1), Io1 = C.IO1_LB; end
    if ~isfinite(Io2), Io2 = C.IO2_LB; end
    if ~isfinite(n1),  n1  = C.N1_LB;  end
    if ~isfinite(n2),  n2  = C.N2_LB;  end

    Rs  = clamp_val(Rs,  C.RS_LB,  C.RS_UB);
    Rsh = clamp_val(Rsh, C.RSH_LB, C.RSH_UB);
    Io1 = clamp_val(Io1, C.IO1_LB, C.IO1_UB);
    Io2 = clamp_val(Io2, C.IO2_LB, C.IO2_UB);
    n1  = clamp_val(n1,  C.N1_LB,  C.N1_UB);
    n2  = clamp_val(n2,  C.N2_LB,  C.N2_UB);

    theta6c = [Rs, Rsh, Io1, Io2, n1, n2];
end

% ---------------------------------------------------------------------
function [rmse_val, Iph_star, theta7] = objective_autoph(theta6, curves, Vt, C)

    theta6 = clamp_theta6(theta6, C);
    Iph_star = iph_opt_closed(theta6, curves, Vt, C);

    Rs  = theta6(1);
    Rsh = theta6(2);
    Io1 = theta6(3);
    Io2 = theta6(4);
    n1  = theta6(5);
    n2  = theta6(6);

    theta7 = [Iph_star, Io1, Io2, Rs, Rsh, n1, n2];
    rmse_val = objective_rmse_obj(theta7, curves, Vt, C);
end

% ---------------------------------------------------------------------
function Iph_star = iph_opt_closed(theta6, curves, Vt, C)

    Rs  = theta6(1);
    Rsh = theta6(2);
    Io1 = theta6(3);
    Io2 = theta6(4);
    n1  = theta6(5);
    n2  = theta6(6);

    denom1 = max(n1 * Vt, C.EPS);
    denom2 = max(n2 * Vt, C.EPS);
    Rsh_d  = max(Rsh, C.EPS);

    V = curves{1};
    I = curves{2};

    u1 = (V + I .* Rs) ./ denom1;
    u2 = (V + I .* Rs) ./ denom2;
    u1 = min(max(u1, -100), 100);
    u2 = min(max(u2, -100), 100);

    term = Io1 .* (exp(u1) - 1.0) + Io2 .* (exp(u2) - 1.0) + ...
           (V + I .* Rs) ./ Rsh_d + I;

    Iph_star = mean(term);
    Iph_star = clamp_val(Iph_star, C.IPH_LB, C.IPH_UB);
end

% ---------------------------------------------------------------------
function buckets = push_topk_stratified(buckets, rmse, iph, theta6, C, k_per_bucket)

    if isempty(rmse) || ~isfinite(rmse)
        return
    end

    theta6 = clamp_theta6(theta6, C);

    Rsh = theta6(2);
    Io2 = theta6(4);
    n2  = theta6(6);

    if Rsh < 50
        b_rsh = 0;
    elseif Rsh < 500
        b_rsh = 1;
    else
        b_rsh = 2;
    end

    logIo2 = log10(max(Io2, C.IO2_LB));
    if logIo2 < -10
        b_io2 = 0;
    elseif logIo2 < -8
        b_io2 = 1;
    else
        b_io2 = 2;
    end

    if n2 < 1.3
        b_n2 = 0;
    elseif n2 < 1.7
        b_n2 = 1;
    else
        b_n2 = 2;
    end

    key = sprintf('%d_%d_%d', b_rsh, b_io2, b_n2);

    item = struct('rmse',double(rmse),'iph',double(iph),'theta6',theta6);

    if ~isKey(buckets, key)
        buckets(key) = {item};
    else
        lst = buckets(key);
        lst{end+1} = item;

        rmse_list = cellfun(@(s) s.rmse, lst);
        [~, idx] = sort(rmse_list, 'ascend');
        lst = lst(idx);

        if numel(lst) > k_per_bucket
            lst = lst(1:k_per_bucket);
        end
        buckets(key) = lst;
    end
end

% ---------------------------------------------------------------------
function y = clamp_val(x, lb, ub)
    y = min(max(x, lb), ub);
end