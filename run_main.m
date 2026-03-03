function out = run_main(varargin)
%RUN_MAIN Reproducible CLI entry-point for PV_DDM_Project (no GUI).
%
% Examples:
%   run_main("module","config/modules/STM6_40_36.txt", ...
%            "iv","data/iv/example_STM6_40_36.txt", ...
%            "seed",42);
%
%   run_main("module","config/modules/STM6_40_36.txt", ...
%            "iv","data/iv/example_STM6_40_36.txt", ...
%            "batchN",30,"consecutiveSeeds",true,"seed",42);

% -----------------------------
% 0) Repo root and defaults
% -----------------------------
repoRoot = fileparts(mfilename("fullpath"));

% --- Ensure repo is on MATLAB path (robust) ---
if exist('repoRoot','var') && isfolder(repoRoot)
    addpath(repoRoot);
    addpath(genpath(fullfile(repoRoot,'core')));
    rehash;
else
    error('repoRoot is not a valid folder.');
end

cfg = default_cfg(repoRoot);
if exist(fullfile(repoRoot,"config","cfg_default.m"),"file") == 2
    try
        cfg2 = feval("cfg_default"); %#ok<FVAL>
        if isstruct(cfg2), cfg = merge_struct(cfg, cfg2); end
    catch
        % keep fallback cfg
    end
end

% -----------------------------
% 1) Parse inputs
% -----------------------------
ip = inputParser;
ip.FunctionName = "run_main";

addParameter(ip,"module", cfg.module, @(x)ischar(x)||isstring(x));
addParameter(ip,"iv",     cfg.iv,     @(x)ischar(x)||isstring(x));

addParameter(ip,"outRoot", cfg.outRoot, @(x)ischar(x)||isstring(x));
addParameter(ip,"seed",    cfg.seed,    @(x)isnumeric(x)&&isscalar(x)&&isfinite(x));
addParameter(ip,"batchN",  cfg.batchN,  @(x)isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>=1);
addParameter(ip,"consecutiveSeeds", cfg.consecutiveSeeds, @(x)islogical(x)||ismember(x,[0 1]));
addParameter(ip,"appVersion", cfg.appVersion, @(x)ischar(x)||isstring(x));
addParameter(ip,"dryRun", false, @(x)islogical(x)||ismember(x,[0 1]));

parse(ip,varargin{:});
opt = ip.Results;

modulePath = normalize_path(repoRoot, string(opt.module));
ivPath     = normalize_path(repoRoot, string(opt.iv));
outRoot    = normalize_path(repoRoot, string(opt.outRoot));

if ~isfile(modulePath)
    error("run_main:FileNotFound","Module file not found: %s", modulePath);
end
if ~isfile(ivPath)
    error("run_main:FileNotFound","I-V file not found: %s", ivPath);
end

batchN = round(double(opt.batchN));
seed0  = double(opt.seed);
useConsecutive = logical(opt.consecutiveSeeds);

% -----------------------------
% 2) Load inputs (CLI does not depend on GUI)
% -----------------------------
[panel, moduleMeta] = load_module_txt(modulePath);
[meta, V_meas, I_meas, ivMeta] = load_iv_txt(ivPath);

meta.Seed = seed0; % requested seed (effective seed handled by core/pipeline)
meta.ModuleFile = string(modulePath);
meta.IVFile     = string(ivPath);

% -----------------------------
% 3) Prepare output folder + logger
% -----------------------------
ts = char(datetime("now","Format","yyyy-MM-dd_HH-mm-ss"));
modName = string(moduleMeta.name);
if strlength(modName)==0
    [~,modName] = fileparts(modulePath);
end
ivID = string(meta.ID);
if strlength(ivID)==0
    [~,ivID] = fileparts(ivPath);
end

tag = ts + "_" + sanitize(modName) + "_" + sanitize(ivID);
if batchN > 1
    tag = tag + "_batch" + string(batchN);
else
    tag = tag + "_seed" + string(seed0);
end

outDir = fullfile(outRoot,"runs",char(tag));
if ~exist(outDir,"dir")
    mkdir(outDir);
end

logFile = fullfile(outDir,"log.txt");
logFcn  = make_logger(logFile);

logFcn("INFO", sprintf("Repo root: %s", repoRoot));
logFcn("INFO", sprintf("Module: %s", modulePath));
logFcn("INFO", sprintf("IV: %s", ivPath));
logFcn("INFO", sprintf("IV ID: %s | Points: %d | G=%.3f | T=%.3f", meta.ID, numel(V_meas), meta.G, meta.T));
logFcn("INFO", sprintf("Requested seed: %.0f | batchN=%d | consecutive=%d", seed0, batchN, useConsecutive));

if opt.dryRun
    logFcn("WARN","Dry run enabled. Exiting without running the solver.");
    out = struct("outDir",outDir,"dryRun",true);
    return
end

% -----------------------------
% 4) Run (prefer pipeline if available)
% -----------------------------
t0 = tic;

if batchN <= 1
    runout = run_single(repoRoot, panel, V_meas, I_meas, meta, logFcn);
    elapsed = toc(t0);

    % Save results
    save(fullfile(outDir,"run_results.mat"),"runout","-v7.3");

    % Save summary CSV
    summary = summarize_single(runout, elapsed);
    writetable(struct2table(summary), fullfile(outDir,"run_summary.csv"));

    % Build run_config.json
    runCfg = build_run_config(repoRoot, opt, modulePath, ivPath, meta, panel, runout);
    write_json(fullfile(outDir,"run_config.json"), runCfg);

    % Checksums
    write_checksums(fullfile(outDir,"checksums.sha256"), repoRoot, ...
        {modulePath, ivPath, ...
         fullfile(repoRoot,"core","RS7MSA.m"), ...
         fullfile(repoRoot,"core","mutate_theta6_phase1.m"), ...
         fullfile(repoRoot,"core","mutate_theta6_phase2.m"), ...
         fullfile(repoRoot,"core","utils_rng.m")}, ...
        logFcn);

    logFcn("INFO", sprintf("Done. Elapsed: %.3f s", elapsed));

    out = struct("outDir",outDir,"runout",runout,"summary",summary);

else
    batch = run_batch(repoRoot, panel, V_meas, I_meas, meta, batchN, seed0, useConsecutive, logFcn);
    elapsed = toc(t0);

    save(fullfile(outDir,"run_results.mat"),"batch","-v7.3");

    % Per-run table
    writetable(batch.runsTable, fullfile(outDir,"runs_table.csv"));

    % Batch summary CSV
    writetable(struct2table(batch.summary), fullfile(outDir,"run_summary.csv"));

    runCfg = build_run_config(repoRoot, opt, modulePath, ivPath, meta, panel, batch);
    write_json(fullfile(outDir,"run_config.json"), runCfg);

    write_checksums(fullfile(outDir,"checksums.sha256"), repoRoot, ...
        {modulePath, ivPath, ...
         fullfile(repoRoot,"core","RS7MSA.m"), ...
         fullfile(repoRoot,"core","mutate_theta6_phase1.m"), ...
         fullfile(repoRoot,"core","mutate_theta6_phase2.m"), ...
         fullfile(repoRoot,"core","utils_rng.m")}, ...
        logFcn);

    logFcn("INFO", sprintf("Done (batch). Elapsed: %.3f s", elapsed));

    out = struct("outDir",outDir,"batch",batch);

end

end

% =========================================================================
% Local helpers
% =========================================================================

function cfg = default_cfg(repoRoot)
cfg = struct();
cfg.module = fullfile(repoRoot,"config","modules","STM6_40.txt");
cfg.iv     = fullfile(repoRoot,"data","iv","CurvasIV_STM6_40.txt");
cfg.outRoot = fullfile(repoRoot,"outputs");
cfg.seed = 42;
cfg.batchN = 1;
cfg.consecutiveSeeds = true;
cfg.appVersion = "dev";
end

function s = merge_struct(a,b)
s = a;
fn = fieldnames(b);
for i = 1:numel(fn)
    s.(fn{i}) = b.(fn{i});
end
end

function p = normalize_path(repoRoot, pIn)
pIn = string(pIn);
if strlength(pIn)==0
    p = "";
    return
end
% allow relative to repo
if ~isfile(pIn) && ~isfolder(pIn)
    pTry = fullfile(repoRoot, char(pIn));
    p = string(pTry);
else
    p = pIn;
end
end

function name = sanitize(s)
s = string(s);
s = regexprep(s,"[^A-Za-z0-9_\-]+","_");
s = regexprep(s,"_+","_");
name = s;
end

function logFcn = make_logger(logFile)
fid = fopen(logFile,"a");
if fid < 0
    error("run_main:Log","Cannot open log file: %s", logFile);
end
cleanupObj = onCleanup(@() fclose(fid)); %#ok<NASGU>

    function write(level,msg)
        ts = char(datetime("now","Format","yyyy-MM-dd HH:mm:ss"));
        line = sprintf("[%s] %-5s %s", ts, upper(string(level)), string(msg));
        fprintf(fid,"%s\n",line);
        % CLI feedback is ok
        fprintf("%s\n",line);
    end

logFcn = @write;
end

function [panel, metaOut] = load_module_txt(fpath)
% Parse simple "name = value;" TXT (MATLAB-like assignments).
lines = string(readlines(fpath));
M = struct();
for i = 1:numel(lines)
    s = strtrim(lines(i));
    if s=="" || startsWith(s,"%") || startsWith(s,"%%")
        continue
    end
    pct = strfind(s,"%");
    if ~isempty(pct)
        s = strtrim(extractBefore(s,pct(1)));
    end
    if s==""; continue; end

    tok = regexp(s,'^(?<name>[A-Za-z]\w*)\s*=\s*(?<rhs>.+?)\s*;\s*$','names');
    if isempty(tok); continue; end

    key = tok.name;
    rhs = strtrim(string(tok.rhs));

    if startsWith(rhs,"'") && endsWith(rhs,"'")
        M.(key) = char(extractBetween(rhs,2,strlength(rhs)-1));
        continue
    end
    if strcmpi(rhs,"NaN")
        M.(key) = NaN;
        continue
    end
    num = str2double(rhs);
    if ~isnan(num)
        M.(key) = num;
    else
        M.(key) = char(rhs);
    end
end

% Minimal mapping to the core "panel" struct
panel = struct();
panel.modelo = "";
metaOut = struct();
[~,metaOut.name] = fileparts(fpath);

% Common names in your module files
panel.Pmax_W   = pick_num(M, {"Pmax_W","Pmax","PmaxE","Pmpp"}, NaN);
panel.Voc_ref  = pick_num(M, {"Voc_val","Voc","Vocn"}, NaN);
panel.Isc      = pick_num(M, {"Isc_val","Isc","Iscn"}, NaN);
panel.Vmpp_ref = pick_num(M, {"Vmpp_val","Vmpp","Vmp"}, NaN);
panel.Impp_ref = pick_num(M, {"Impp_val","Impp","Imp"}, NaN);
panel.Ns       = pick_num(M, {"Ns"}, NaN);
panel.N_pv     = pick_num(M, {"Np","N_pv"}, NaN);

panel.Gref     = pick_num(M, {"Gref"}, 1000);
panel.Tref_C   = pick_num(M, {"Tref"}, 25);

panel.kIsc     = pick_num(M, {"Ki","kIsc"}, NaN);
panel.kVoc_pct = pick_num(M, {"Kv","kVoc_pct"}, NaN);

if isfield(M,"modelo"), panel.modelo = string(M.modelo); end
if strlength(panel.modelo)==0, panel.modelo = string(metaOut.name); end

end

function v = pick_num(M, names, defaultVal)
v = defaultVal;
for k = 1:numel(names)
    nm = string(names{k});
    if isfield(M, nm)
        raw = M.(nm);
        if isnumeric(raw) && isscalar(raw) && isfinite(raw)
            v = double(raw); return
        end
        if ischar(raw) || isstring(raw)
            x = str2double(string(raw));
            if isfinite(x)
                v = x; return
            end
        end
    end
end
end

function [meta,V,I,metaOut] = load_iv_txt(filename)
% Format:
%   % ID=...
%   % G=...
%   % T=...
%   % Date=...
%   V,I lines as "V,I"
meta = struct("ID","","G",NaN,"T",NaN,"Date","");
fid = fopen(filename,"r");
if fid < 0
    error("run_main:IV","Cannot open file: %s", filename);
end
c = onCleanup(@() fclose(fid));

data = zeros(0,2);
while ~feof(fid)
    line = strtrim(fgetl(fid));
    if startsWith(line,"%")
        if contains(line,"ID="),   meta.ID   = strtrim(extractAfter(line,"=")); end
        if contains(line,"G="),    meta.G    = str2double(extractAfter(line,"=")); end
        if contains(line,"T="),    meta.T    = str2double(extractAfter(line,"=")); end
        if contains(line,"Date="), meta.Date = strtrim(extractAfter(line,"=")); end
        continue
    end
    if line=="" || contains(lower(line),"v_meas")
        continue
    end
    vals = sscanf(line,"%f,%f");
    if numel(vals)==2
        data(end+1,:) = vals(:).'; %#ok<AGROW>
    end
end

V = data(:,1);
I = data(:,2);

if strlength(string(meta.ID))==0
    [~,name] = fileparts(filename);
    meta.ID = name;
end

metaOut = meta;
end

function runout = run_single(repoRoot, panel, V_meas, I_meas, meta, logFcn)
% Prefer pipeline if present
if exist(fullfile(repoRoot,"pipeline","pipeline_single_run.m"),"file") == 2 || exist("pipeline_single_run","file")==2
    try
        runout = pipeline_single_run(panel, V_meas, I_meas, meta, logFcn); %#ok<NASGU>
        return
    catch
        % fallback to core
    end
end
runout = RS7MSA(panel, V_meas, I_meas, meta, logFcn);
end

function batch = run_batch(repoRoot, panel, V_meas, I_meas, meta, N, seed0, useConsecutive, logFcn)
% Prefer pipeline if present
if exist(fullfile(repoRoot,"pipeline","pipeline_batch_run.m"),"file") == 2 || exist("pipeline_batch_run","file")==2
    try
        batch = pipeline_batch_run(panel, V_meas, I_meas, meta, N, seed0, useConsecutive, logFcn); %#ok<NASGU>
        return
    catch
        % fallback to local implementation
    end
end

baseSeed = seed0;
if useConsecutive
    seeds = baseSeed + (0:N-1);
else
    rng(baseSeed,"twister");
    seeds = randi(1e9,[N,1]);
end

rmse_obj = nan(N,1);
tRun     = nan(N,1);
status   = strings(N,1);

for i = 1:N
    meta_i = meta;
    meta_i.Seed = double(seeds(i));
    t0 = tic;
    try
        runout_i = RS7MSA(panel, V_meas, I_meas, meta_i, []);
        tRun(i) = toc(t0);
        if isfield(runout_i,"rmse_obj") && isfinite(runout_i.rmse_obj)
            rmse_obj(i) = runout_i.rmse_obj;
            status(i) = "OK";
        else
            status(i) = "FAIL";
        end
    catch
        tRun(i) = toc(t0);
        status(i) = "FAIL";
    end

    if mod(i,max(1,ceil(N/10)))==0 || i==N
        logFcn("INFO", sprintf("Batch %d/%d | last status=%s", i, N, status(i)));
    end
end

ok = isfinite(rmse_obj);
summary = struct();
summary.model = string(panel.modelo);
summary.N = sum(ok);
summary.best  = min(rmse_obj(ok));
summary.worst = max(rmse_obj(ok));
summary.mean  = mean(rmse_obj(ok));
summary.std   = std(rmse_obj(ok),0);
summary.time_s_run = mean(tRun,"omitnan");

runsTable = table((1:N).', seeds(:), rmse_obj(:), tRun(:), status(:), ...
    "VariableNames", {"run","seed","rmse_obj","time_s","status"});

batch = struct();
batch.summary = summary;
batch.runsTable = runsTable;
batch.seeds = seeds(:);
batch.rmse_obj = rmse_obj(:);
batch.tRun = tRun(:);
end

function summary = summarize_single(runout, elapsed)
summary = struct();
summary.elapsed_s = elapsed;

% prefer rmse_obj if present
summary.rmse_obj = getfield_def(runout,"rmse_obj", NaN);
summary.rmse_I   = getfield_def(runout,"rmse_I",   NaN);

% include seed if stored by core
summary.seed_effective = getfield_def(runout,"seed", NaN);
end

function v = getfield_def(s, field, defaultVal)
if isstruct(s) && isfield(s, field)
    vv = s.(field);
    if isnumeric(vv) && isscalar(vv) && isfinite(vv)
        v = double(vv);
    else
        v = defaultVal;
    end
else
    v = defaultVal;
end
end

function cfg = build_run_config(repoRoot, opt, modulePath, ivPath, meta, panel, resultStruct)
cfg = struct();

cfg.timestamp = char(datetime("now","Format","yyyy-MM-dd HH:mm:ss"));
cfg.app_version = char(string(opt.appVersion));
cfg.matlab_version = version;
cfg.repo_root = char(repoRoot);

cfg.inputs = struct();
cfg.inputs.module_file = char(modulePath);
cfg.inputs.iv_file     = char(ivPath);
cfg.inputs.iv_id       = char(string(meta.ID));
cfg.inputs.G_Wm2       = double(meta.G);
cfg.inputs.T_C         = double(meta.T);

cfg.requested = struct();
cfg.requested.seed  = double(opt.seed);
cfg.requested.batchN = double(opt.batchN);
cfg.requested.consecutiveSeeds = logical(opt.consecutiveSeeds);

cfg.panel = panel;

% result summary
cfg.output = struct();
if isfield(resultStruct,"rmse_obj"), cfg.output.rmse_obj = resultStruct.rmse_obj; end
if isfield(resultStruct,"rmse_I"),   cfg.output.rmse_I   = resultStruct.rmse_I;   end
if isfield(resultStruct,"seed"),     cfg.output.seed_effective = resultStruct.seed; end

cfg.git_commit = try_git_commit(repoRoot);

end

function commit = try_git_commit(repoRoot)
commit = "";
try
    cwd = pwd;
    c = onCleanup(@() cd(cwd));
    cd(repoRoot);
    [st, out] = system("git rev-parse HEAD");
    if st==0
        commit = strtrim(string(out));
    end
catch
    commit = "";
end
end

function write_json(filePath, s)
txt = jsonencode(s);
txt = pretty_json(txt);
fid = fopen(filePath,"w");
if fid < 0
    error("run_main:JSON","Cannot write: %s", filePath);
end
c = onCleanup(@() fclose(fid));
fprintf(fid,"%s\n",txt);
end

function txt = pretty_json(txt)
% Minimal pretty printer: keeps jsonencode output readable enough
% without external dependencies.
txt = string(txt);
txt = replace(txt, "},{", "},\n{");
txt = replace(txt, ',"', sprintf(',\n  "'));
txt = replace(txt, "{", "{\n  ");
txt = replace(txt, "}", "\n}");
txt = char(txt);
end

function write_checksums(outFile, repoRoot, fileList, logFcn)
fileList = fileList(~cellfun(@isempty,fileList));
fid = fopen(outFile,"w");
if fid < 0
    logFcn("WARN", sprintf("Cannot write checksums file: %s", outFile));
    return
end
c = onCleanup(@() fclose(fid));

for i = 1:numel(fileList)
    f = string(fileList{i});
    if ~isfile(f)
        continue
    end
    h = sha256_file(f);
    rel = make_relative(repoRoot, f);
    fprintf(fid,"%s  %s\n", h, rel);
end
end

function rel = make_relative(root, p)
root = char(string(root));
p    = char(string(p));
rel = p;
if startsWith(p, root)
    rel = p(numel(root)+2:end);
end
end

function h = sha256_file(filePath)
% No toolboxes. Uses Java MessageDigest.
md = java.security.MessageDigest.getInstance("SHA-256");
fis = java.io.FileInputStream(java.io.File(char(filePath)));
dis = java.security.DigestInputStream(fis, md);
buf = zeros(1,8192,"int8");
while dis.read(buf) ~= -1
end
dis.close();
hash = typecast(md.digest(),"uint8");
h = lower(reshape(dec2hex(hash,2).',1,[]));
end