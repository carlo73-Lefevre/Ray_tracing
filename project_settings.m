function settings = project_settings(varargin)
%PROJECT_SETTINGS  Load project directory configuration.
%   SETTINGS = PROJECT_SETTINGS() returns a structure with commonly used
%   directories for the ray tracing project. The function looks for a
%   'config.json' file located in the same folder and reads the field
%   'root' to determine the root directory. The root can also be passed as
%   a name-value pair argument:
%       SETTINGS = PROJECT_SETTINGS('root', ROOTDIR)
%
%   The returned structure contains the following fields:
%       .root      - project root directory
%       .functions - path to the 'Funzioni' folder
%       .data      - path to 'Dati Caomsol'
%       .class     - path to 'Class'
%       .results   - path to 'Results'
%       .boxWingResults - path to 'Box_wing/Confronto Box Wing Surf .../Massimo_risultati'
%
%   Using this helper avoids hard-coded absolute paths inside the start
%   scripts.

p = inputParser;
p.addParameter('root','',@ischar);
p.parse(varargin{:});
settings.root = p.Results.root;

if isempty(settings.root)
    cfgFile = fullfile(fileparts(mfilename('fullpath')),'config.json');
    if exist(cfgFile,'file')
        cfg = jsondecode(fileread(cfgFile));
        if isfield(cfg,'root')
            settings.root = cfg.root;
        end
    end
end

if isempty(settings.root)
    settings.root = pwd; % fallback
end

settings.functions = fullfile(settings.root,'Funzioni');
settings.data      = fullfile(settings.root,'Dati Caomsol');
settings.class     = fullfile(settings.root,'Class');
settings.results   = fullfile(settings.root,'Results');
settings.boxWingResults = fullfile(settings.root,'Box_wing', ...
    'Confronto Box Wing Surf Surf_Visco', 'Script di Confronto', ...
    'Massimo_risultati');
end
