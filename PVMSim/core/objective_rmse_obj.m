function rmse = objective_rmse_obj(theta7, curves, Vt, C)
%OBJECTIVE_RMSE_OBJ  RMSE of implicit residual using measured I

    Iph = theta7(1);
    Io1 = theta7(2);
    Io2 = theta7(3);
    Rs  = theta7(4);
    Rsh = theta7(5);
    n1  = theta7(6);
    n2  = theta7(7);

    % Defensive clamps (same as before)
    Io1 = min(max(Io1, C.IO1_LB), C.IO1_UB);
    Io2 = min(max(Io2, C.IO2_LB), C.IO2_UB);
    Rs  = min(max(Rs,  C.RS_LB),  C.RS_UB);
    Rsh = min(max(Rsh, C.RSH_LB), C.RSH_UB);
    n1  = min(max(n1,  C.N1_LB),  C.N1_UB);
    n2  = min(max(n2,  C.N2_LB),  C.N2_UB);

    denom1 = max(n1 * Vt, C.EPS);
    denom2 = max(n2 * Vt, C.EPS);
    Rsh_d  = max(Rsh, C.EPS);

    V = curves{1};
    I = curves{2};

    u1 = (V + I .* Rs) ./ denom1;
    u2 = (V + I .* Rs) ./ denom2;
    u1 = min(max(u1, -100), 100);
    u2 = min(max(u2, -100), 100);

    e1 = exp(u1);
    e2 = exp(u2);

    r = Iph - Io1 .* (e1 - 1.0) - Io2 .* (e2 - 1.0) - ...
        (V + I .* Rs) ./ Rsh_d - I;

    if any(~isfinite(r))
        rmse = 1e9;
        return
    end

    rmse = sqrt(mean(r.^2));
end