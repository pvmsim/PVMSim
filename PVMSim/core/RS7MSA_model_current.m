function Icalc = RS7MSA_model_current(theta7, Vvec, Vt, C)
%RS7MSA_MODEL_CURRENT  Solve I(V) using Newton-Raphson

    Iph = theta7(1);
    Io1 = theta7(2);
    Io2 = theta7(3);
    Rs  = theta7(4);
    Rsh = max(theta7(5), C.EPS);
    n1  = theta7(6);
    n2  = theta7(7);

    denom1 = max(n1 * Vt, C.EPS);
    denom2 = max(n2 * Vt, C.EPS);

    N = numel(Vvec);
    Icalc = zeros(size(Vvec));

    Iprev = Iph;

    for i = 1:N
        V = Vvec(i);

        I = Iprev;

        for it = 1:80
            u1 = (V + I*Rs) / denom1;
            u2 = (V + I*Rs) / denom2;
            u1 = min(max(u1, -100), 100);
            u2 = min(max(u2, -100), 100);

            e1 = exp(u1);
            e2 = exp(u2);

            f  = Iph - Io1*(e1 - 1) - Io2*(e2 - 1) - (V + I*Rs)/Rsh - I;
            df = -Io1*e1*(Rs/denom1) - Io2*e2*(Rs/denom2) - (Rs/Rsh) - 1.0;

            if ~isfinite(f) || ~isfinite(df) || abs(df) < 1e-14
                break
            end

            Inew = I - f/df;

            if ~isfinite(Inew)
                break
            end

            if abs(Inew - I) < 1e-10
                I = Inew;
                break
            end

            I = Inew;
        end

        Icalc(i) = I;
        Iprev = I;
    end
end