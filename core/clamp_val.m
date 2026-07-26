function y = clamp_val(x, lb, ub)
    y = min(max(x, lb), ub);
end