clearvars; close all; opengl hardware

% Clone example data
if ~exist('VSD', 'dir')
    try
    !git clone https://github.com/RWTHmediTEC/VSDFullBodyBoneModels VSD
    rmdir('VSD/.git', 's')
    catch
        warning([newline 'Clone (or copy) the example data from: ' ...
            'https://github.com/RWTHmediTEC/VSDFullBodyBoneModels' newline 'to: ' ...
            fileparts([mfilename('fullpath'), '.m']) '\VSD' ...
            ' and try again!' newline])
        return
    end
end

addpath(genpath('src'))

% Load subject names
load('VSD\MATLAB\res\VSD_Subjects.mat', 'Subjects')

% Mode
mode='subburaj';

for s=1%:size(Subjects, 1)
    load(['VSD\Bones\' Subjects.Number{s} '.mat'],'B');
    % Construct the pelvic bone
    [pelvis.vertices, pelvis.faces] = concatenateMeshes(B(1:3).mesh);
    % Get APP coordinate system transformation
    APP_TFM = pelvicLandmarkID(pelvis, 'visu',0);
    % Combine the hip bones
    [hipBones.vertices, hipBones.faces] = concatenateMeshes(B(2:3).mesh);
    % ID landmarks
    Landmarks = pelvicLandmarkID_Subburaj(hipBones, 'initalTransform', APP_TFM, 'visu',1, 'mode', mode, 'sym', 1);
    set(gcf, 'Name',['Subject: ' Subjects.Number{s}], 'NumberTitle', 'Off')
end

% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';