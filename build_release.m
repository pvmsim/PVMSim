%% Package PVMSim v0.1.0 without a .prj file
clear;
clc;

repoRoot = string(pwd);

%% 1. Verify repository root
assert(isfile(fullfile(repoRoot, "run_main.m")), ...
    "run_main.m was not found. Set the Current Folder to the PVMSim repository root.");

assert(isfile(fullfile(repoRoot, "app", "PVMSim.mlapp")), ...
    "app/PVMSim.mlapp was not found.");

assert(isfolder(fullfile(repoRoot, "core")), ...
    "The core folder was not found.");

assert(isfolder(fullfile(repoRoot, "config")), ...
    "The config folder was not found.");

assert(isfolder(fullfile(repoRoot, "data")), ...
    "The data folder was not found.");

assert(isfolder(fullfile(repoRoot, "tests")), ...
    "The tests folder was not found.");

%% 2. Obtain or create a permanent toolbox identifier
identifierFile = fullfile(repoRoot, "toolbox_identifier.txt");
identifier = "";

installed = matlab.addons.toolbox.installedToolboxes;

if ~isempty(installed)
    names = string({installed.Name});
    idx = find(strcmpi(names, "PVMSim"), 1);

    if ~isempty(idx)
        identifier = string(installed(idx).Guid);
        fprintf("Using identifier from installed PVMSim toolbox.\n");
    end
end

if strlength(identifier) == 0 && isfile(identifierFile)
    identifier = strtrim(string(fileread(identifierFile)));
    fprintf("Using identifier stored in toolbox_identifier.txt.\n");
end

if strlength(identifier) == 0
    identifier = string(char(java.util.UUID.randomUUID));
    fprintf("A new permanent toolbox identifier was created.\n");
end

% Always store the identifier in the repository
fid = fopen(identifierFile, "w");

if fid < 0
    error("Cannot create toolbox_identifier.txt: %s", identifierFile);
end

fprintf(fid, "%s\n", identifier);
fclose(fid);

assert(isfile(identifierFile), ...
    "toolbox_identifier.txt could not be created.");

fprintf("Toolbox identifier: %s\n", identifier);

%% 3. Select files and folders to package
toolboxFiles = [
    fullfile(repoRoot, "run_main.m")
    fullfile(repoRoot, "README.md")
    fullfile(repoRoot, "LICENSE")
    fullfile(repoRoot, "CITATION.cff")
    fullfile(repoRoot, "toolbox_identifier.txt")
    fullfile(repoRoot, "app")
    fullfile(repoRoot, "config")
    fullfile(repoRoot, "core")
    fullfile(repoRoot, "data")
    fullfile(repoRoot, "docs")
    fullfile(repoRoot, "tests")
    fullfile(repoRoot, "third_party")
];

% Optional repository resources
optionalResources = [
    fullfile(repoRoot, "THIRD_PARTY_NOTICES.md")
    fullfile(repoRoot, "LICENSE-CC0.txt")
    fullfile(repoRoot, "figures")
    fullfile(repoRoot, "Documentation")
];

for k = 1:numel(optionalResources)
    resource = optionalResources(k);

    if isfile(resource) || isfolder(resource)
        toolboxFiles(end + 1, 1) = resource; %#ok<SAGROW>
    end
end

% Verify required files and folders.
for k = 1:numel(toolboxFiles)
    item = toolboxFiles(k);

    assert(isfile(item) || isfolder(item), ...
        "Required toolbox resource not found: %s", item);
end

%% 4. Prepare output folder
releaseDir = fullfile(repoRoot, "Release");

if ~isfolder(releaseDir)
    mkdir(releaseDir);
end

outputFile = fullfile(releaseDir, "PVMSim.mltbx");

if isfile(outputFile)
    delete(outputFile);
end

%% 5. Configure toolbox metadata
opts = matlab.addons.toolbox.ToolboxOptions( ...
    repoRoot, identifier);

opts.ToolboxName = "PVMSim";
opts.ToolboxVersion = "0.1.0";

opts.Summary = ...
    "Reproducible extraction of double-diode photovoltaic model parameters.";

opts.Description = ...
    "PVMSim is a MATLAB application for reproducible parameter " + ...
    "extraction of the double-diode photovoltaic model from measured " + ...
    "current-voltage curves. It includes graphical and command-line " + ...
    "interfaces, benchmark data, controlled random seeds, structured " + ...
    "run artifacts, integrity hashes, batch execution, and an automated " + ...
    "smoke test.";

opts.AuthorName = "Liomnis Osorio et al.";
opts.AuthorEmail = "pvmsim.matlab@gmail.com";

opts.MinimumMatlabRelease = "R2025b";
opts.MaximumMatlabRelease = "";

opts.SupportedPlatforms.Win64 = true;
opts.SupportedPlatforms.Mac = true;
opts.SupportedPlatforms.Glnxa64 = true;
opts.SupportedPlatforms.MatlabOnline = false;

%% 6. Define packaged files and installation paths
opts.ToolboxFiles = toolboxFiles;

opts.ToolboxMatlabPath = [
    repoRoot
    fullfile(repoRoot, "app")
    fullfile(repoRoot, "core")
];

opts.AppGalleryFiles = ...
    fullfile(repoRoot, "app", "PVMSim.mlapp");

opts.OutputFile = outputFile;

%% 7. Package toolbox
matlab.addons.toolbox.packageToolbox(opts);

%% 8. Verify generated package
assert(isfile(outputFile), ...
    "PVMSim.mltbx was not generated.");

generatedVersion = ...
    matlab.addons.toolbox.toolboxVersion(outputFile);

fprintf("\nToolbox generated successfully:\n");
fprintf("%s\n", outputFile);
fprintf("Version: %s\n", generatedVersion);
fprintf("Identifier: %s\n", identifier);