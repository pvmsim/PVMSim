function [seed, rng_state0] = utils_rng(meta, seed_default)
%UTILS_RNG  Reproducibility-aware RNG initialization.
%   - If meta.Seed exists and is finite scalar -> use it.
%   - Else use seed_default (default 42).
%   Returns:
%     seed        : effective seed used (1..2^32-1)
%     rng_state0  : rng() struct after initialization

    if nargin < 2 || isempty(seed_default) || ~isfinite(seed_default)
        seed_default = 42;
    end
    seed_default = double(round(seed_default));
    seed_default = min(max(seed_default, 1), 4294967295);

    seed = seed_default;

    if nargin >= 1 && ~isempty(meta) && isstruct(meta) && isfield(meta,'Seed')
        s = meta.Seed;

        if ischar(s) || isstring(s)
            s = str2double(string(s));
        end

        if isfinite(s) && isscalar(s)
            seed = floor(abs(double(s)));
        end
    end

    % Normalize to valid MT19937 range (1..2^32-1)
    seed = mod(seed, 2^32 - 1);
    if seed == 0
        seed = seed_default;
    end

    rng(seed,'twister');
    rng_state0 = rng;
end