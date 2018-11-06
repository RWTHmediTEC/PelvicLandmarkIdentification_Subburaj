clearvars; close all; opengl hardware
% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';

subjects=dir('data\*.mat');
NoS=length(subjects);

for s=1%:NoS
    
load(fullfile(subjects(s).folder, subjects(s).name), 'pelvis')
Landmarks = PelvicLandmarkIdSubburaj(pelvis, 'visu',1);

end