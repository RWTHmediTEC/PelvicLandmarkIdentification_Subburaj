clearvars; close all; opengl hardware

addpath('src/external/matGeom/matGeom/meshes3d')

dbPath='../../Database/VSD/';

% Mode
mode='subburaj';

subjects=dir([dbPath 'Bones/z*.mat']);
NoS=length(subjects);

autoLM=cell(14,NoS);
autoLM=cellfun(@(x) nan(1,3), autoLM,'uni',0);
for s=1%:NoS
    % Load meshes and build pelvis
    load(fullfile(subjects(s).folder, subjects(s).name),'B','M')
    [~, Hip_R] = ismember('Hip_R',{B.name});
    [~, Hip_L] = ismember('Hip_L',{B.name});
    [~, Sacrum] = ismember('Sacrum',{B.name});
    [pelvis.vertices, pelvis.faces] = concatenateMeshes([B(Hip_R).mesh, B(Hip_L).mesh]);
    % Load APP coordinate system transformation
    load(fullfile([dbPath 'APP_CS'], subjects(s).name),'APP_TFM')
    Landmarks = pelvicLandmarkIdSubburaj(pelvis, 'initalTransform', APP_TFM, 'visu',1, 'mode', mode, 'sym', 1);
    
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';