clearvars; close all; opengl hardware
% [List.f, List.p] = matlab.codetools.requiredFilesAndProducts([mfilename '.m']); 
% List.f = List.f'; List.p = List.p';

subjects=dir('data/*.mat');
NoS=length(subjects);

autoLM=cell(14,NoS);
for s=1:NoS
    
load(fullfile(subjects(s).folder, subjects(s).name), 'pelvis')
Landmarks = pelvicLandmarkIdSubburaj(pelvis, 'visu',1, 'curv', 40, 'mode', 'full', 'sym', 1);

autoLM{1,s} =nan(1,3);
autoLM{2,s} =nan(1,3);
autoLM{3,s} =Landmarks(1).centroids(2,:); % ASIS_L
autoLM{4,s} =Landmarks(1).centroids(1,:); % ASIS_R
autoLM{5,s} =Landmarks(2).centroids(2,:); % AIIS_L
autoLM{6,s} =Landmarks(2).centroids(1,:); % AIIS_R
autoLM{7,s} =Landmarks(6).centroids(2,:); % PT_L
autoLM{8,s} =Landmarks(6).centroids(1,:); % PT_R
autoLM{9,s} =Landmarks(3).centroids(2,:); % PSIS_L
autoLM{10,s}=Landmarks(3).centroids(1,:); % PSIS_R
autoLM{11,s}=Landmarks(4).centroids(2,:); % PIIS_L
autoLM{12,s}=Landmarks(4).centroids(1,:); % PIIS_R
autoLM{13,s}=Landmarks(5).centroids(2,:); % IS_L
autoLM{14,s}=Landmarks(5).centroids(1,:); % IS_R

end

detected=cellfun(@(x) ~all(isnan(x)), autoLM);