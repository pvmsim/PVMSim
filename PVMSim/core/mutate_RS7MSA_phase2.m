function theta6_out = mutate_RS7MSA_phase2(theta6_in, C, p_big)
%MUTATE_RS7MSA_PHASE2  Phase-2 mutation
% theta6 = [Rs, Rsh, Io1, Io2, n1, n2]

    if nargin < 3, p_big = 0.12; end

    Rs  = theta6_in(1);
    Rsh = theta6_in(2);
    Io1 = theta6_in(3);
    Io2 = theta6_in(4);
    n1  = theta6_in(5);
    n2  = theta6_in(6);

    rs_range  = max(C.RS_UB  - C.RS_LB,  1e-12);
    rsh_range = max(C.RSH_UB - C.RSH_LB, 1e-12);
    n_range   = max(C.N1_UB  - C.N1_LB,  1e-12);

    if rand < p_big
        d_rs  = 0.005 * rs_range;
        d_rsh = 0.01  * rsh_range;
        d_n   = 0.002 * n_range;
        d_log = 1.2;
    else
        d_rs  = 0.001 * rs_range;
        d_rsh = 0.01  * rsh_range;
        d_n   = 0.001 * n_range;
        d_log = 0.2;
    end

    m = @() (randi(3) - 2); % {-1,0,1}

    Rs2  = clamp_val(Rs  + m()*d_rs,  C.RS_LB,  C.RS_UB);
    Rsh2 = clamp_val(Rsh + m()*d_rsh, C.RSH_LB, C.RSH_UB);
    n1_2 = clamp_val(n1  + m()*d_n,   C.N1_LB,  C.N1_UB);
    n2_2 = clamp_val(n2  + m()*d_n,   C.N2_LB,  C.N2_UB);

    if (~isfinite(Io1)) || (Io1 <= 0), Io1 = C.IO1_LB; end
    if (~isfinite(Io2)) || (Io2 <= 0), Io2 = C.IO2_LB; end

    logIo1 = log10(max(Io1, C.IO1_LB));
    logIo2 = log10(max(Io2, C.IO2_LB));

    logIo1_2 = min(max(logIo1 + m()*d_log, log10(C.IO1_LB)), log10(C.IO1_UB));
    logIo2_2 = min(max(logIo2 + m()*d_log, log10(C.IO2_LB)), log10(C.IO2_UB));

    Io1_2 = 10^(logIo1_2);
    Io2_2 = 10^(logIo2_2);

    theta6_out = [Rs2, Rsh2, Io1_2, Io2_2, n1_2, n2_2];

end

function y = clamp_val(x, lb, ub)
    y = min(max(x, lb), ub);
end