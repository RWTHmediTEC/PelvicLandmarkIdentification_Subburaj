function Landmarks = pelvicLandmarkIdSubburaj(pelvis, varargin)
%CURVATUREANALYSIS detects boney landmarks of the pelvis
%
% REQUIRED INPUT:
%   pelvis: A mesh of the pelvis (hip bones and sacrum) consisting of one
%       component with the fields: pelvis.vertices, pelvis.faces
%       ATTENTION: The mesh has to be transformed into the automatic pelvic
%       coordiante system [Kai 2014]. Use the function: automaticPelvicCS.m
%       Then the pubic symphisis (PS) is the origin of the CS.
% OPTIONAL INPUT:
%   'visualization': true (default) or false
%
% OUTPUT: A struct with the following fields:
%   SubMesh: two meshes of the surface regions. 
%   Centroids: centroids of the surface regions.
%   Area: area of the surface regions.
%   ID: ids of the corresponding centroid vertices.
%   Name: The short name of the landmark.
% The first row is the right side, the second row is the left side.
%
% Depending on the spatial relationship matrix following landmarks can be
% identified:
%       ASIS = Anterior Superior Iliac Spine
%       AIIS = Anterior Inferior Iliac Spine
%       PSIS = Posterior Superior Iliac Spine
%       PIIS = Posterior Inferior Iliac Spine
%       IS   = Ischial Spine
%       PT   = Pubic Tubercle
%       IPY  = Iliac Pubic Symphysis
%       IIT  = Iliac Ischial Tuberosity
%       IT   = Iliac Tubercle
%
% REFERENCES:
%   2008 - Subburaj et al. - 3D Shape Reasoning for Identifying Anatomical Landmarks
%   2009 - Subburaj et al. - Automated identification of anatomical landmarks on 3D bone models reconstructed from CT scan images
%
%   TO-DO / IDEAS:
%   Clean up & standardize code
%
% AUTHOR: Felix Krooss
% 	RWTH Aachen University
% VERSION: 1.0
% DATE: 2017-06-14

addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'src']))
addpath(genpath([fileparts([mfilename('fullpath'), '.m']) '\' 'rsc']))

% settings
expectedMode = {'subburaj','small','full', 'rm'};
defaultMode = {'small'};
p = inputParser;
logParValidFunc=@(x) (islogical(x) || isequal(x,1) || isequal(x,0));
addParameter(p,'visualization', false, logParValidFunc);
% Choose spatial adjacency matrix
addParameter(p,'mode', defaultMode, @(x) any(validatestring(lower(x),expectedMode)));
% Choose if identical landmarks have to be on opposite sides
addParameter(p,'symmetry', false, logParValidFunc); 
addParameter(p,'curvatureThreshold', 40, @(x) validateattributes(x, {'numeric'},{'scalar','>', 0,'<', 50}));
addParameter(p,'initalTransform', nan, @isTransform3d)
addParameter(p,'meanTransform', false, logParValidFunc)
parse(p,varargin{:});

visu = p.Results.visualization;
curvThreshold = p.Results.curvatureThreshold;
mode = p.Results.mode;
sym = p.Results.symmetry;
iTFM = p.Results.initalTransform;
meanTFM = p.Results.meanTransform;

%% Pre processing
if ~any(isnan(iTFM))
    pelvis = transformPoint3d(pelvis, iTFM);
end

if meanTFM
    pelvis = transformPoint3d(pelvis, createTranslation3d(-mean(pelvis.vertices)));
end

%% Curvature Analysis

% Compute Gaussian- and Mean Curvature
curvatureOptions.curvature_smoothing = 0;
curvatureOptions.verb = 0;
[~,~,~,~,cMean,cGauss,~] = ...
    compute_curvature(pelvis.vertices,pelvis.faces,curvatureOptions);

% Group vertices according to mean curvature
C_l = cMean<=prctile(cMean,  0 + curvThreshold);
C_h = cMean>=prctile(cMean,100 - curvThreshold);

% Second level grouping according to Gaussian curvature
% 1,4 negative Gaussian, 2,5 Gaussian=0, 3,6 positive Gaussian
G_1 = C_l & (cGauss<0);
G_2 = C_l & (cGauss==0);
G_3 = C_l & (cGauss>0);

G_4 = C_h & (cGauss<0);
G_5 = C_h & (cGauss==0);
G_6 = C_h & (cGauss>0);

% Create seed-triangle list
% Get triangles with identical grouping of vertices
S_1 = sum(G_1(pelvis.faces),2) == 3;
S_2 = sum(G_2(pelvis.faces),2) == 3;
S_3 = sum(G_3(pelvis.faces),2) == 3;
S_4 = sum(G_4(pelvis.faces),2) == 3;
S_5 = sum(G_5(pelvis.faces),2) == 3;
S_6 = sum(G_6(pelvis.faces),2) == 3;

% Get faces of corresponding triangles
S_1f = pelvis.faces(S_1,:);
S_2f = pelvis.faces(S_2,:);
S_3f = pelvis.faces(S_3,:);
S_4f = pelvis.faces(S_4,:);
S_5f = pelvis.faces(S_5,:);
S_6f = pelvis.faces(S_6,:);

if visu
    MonitorsPos = get(0,'MonitorPositions');
    figHandle = figure('Units','pixels','renderer','opengl', 'Color', 'w');
    if     size(MonitorsPos,1) == 1
        set(figHandle,'OuterPosition',[1 50 MonitorsPos(1,3)-1 MonitorsPos(1,4)-50]);
    elseif size(MonitorsPos,1) == 2
        set(figHandle,'OuterPosition',[1+MonitorsPos(1,3) 50 MonitorsPos(2,3)-1 MonitorsPos(2,4)-50]);
    end
    H_Light(1) = light; light('Position', -1*(get(H_Light(1),'Position')));
    axis on; axis equal; hold on
    xlabel x; ylabel y; zlabel z;
    
%     % Color Code Mesh according to curvature
%     cSorted = sort(cMean);
%     idxs = arrayfun(@(x)find(cSorted==x,1),cMean);
%     cmap = jet(length(cMean));
%     cmapS = cmap(idxs,:);
%     patch(pelvis,'FaceColor', 'interp', 'FaceVertexCData', cmapS, 'EdgeColor', 'none');
    
    % Alternative visualization
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [0.75 0.75 0.75];
    patchProps.FaceAlpha = 0.75;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    patch(pelvis, patchProps)
end

%% Create surface regions from seed list

% Pit region
% Counter for number of surface regions
rid=1;
pitRegion = struct('faces',{});
% Create surface regions while seed list is not empty
while ~isempty(S_3f)
    % Counter for number of triangles in surface region
    counter = 1;
    % Use the first triangle of the seed list as reference
    % Get all edge connected faces of the reference triangle
    tempList=sum(ismember(S_3f,S_3f(1,:)),2)>=2;
    % Create surface region
    % Insert reference triangle and edge connected faces
    pitRegion(rid).faces=S_3f(tempList,:);
    % Recursive search for edge connected triangles of neighboring faces
    % while the number of triangles in surface region is bigger than
    % counter variable (number in previous step), repeat
    while sum(tempList) > counter
        counter = sum(tempList);
        % Get all edge connected faces of triangles in surface region
        tempList=sum(ismember(S_3f,pitRegion(rid).faces),2)>=2;
        % Insert newly found triangles in surface region
        pitRegion(rid).faces=S_3f(tempList,:);
    end
    rid=rid+1;
    % Delete all surface triangles from seed list
    S_3f(tempList,:)=[];
end

for i=1:length(pitRegion)
    rFaces=sum(ismember(pelvis.faces,pitRegion(i).faces),2) == 3;
    subMesh = removeMeshFaces(pelvis, ~rFaces);
    % Extract data from region
    pitRegion(i).vertices = subMesh.vertices;
    pitRegion(i).area =  meshSurfaceArea(subMesh.vertices,subMesh.faces);
    pitRegion(i).centroids = mean(subMesh.vertices);
end

% Peak region
% Counter for number of surface regions
rid=1;
peakRegion = struct('faces',{});
% Create surface regions while seed list is not empty
while ~isempty(S_6f)
    % Counter for number of triangles in surface region
    counter = 1;
    % Use the first triangle of the seed list as reference
    % Get all edge connected faces of the reference triangle
    tempList=sum(ismember(S_6f,S_6f(1,:)),2)>=2;
    % Create surface region
    % Insert reference triangle and edge connected faces
    peakRegion(rid).faces=S_6f(tempList,:);
    % Recursive search for edge connected triangles of neighboring faces
    % while the number of triangles in surface region is bigger than
    % counter variable (number in previous step), repeat
    while sum(tempList) > counter
        counter = sum(tempList);
        % Get all edge connected faces of triangles in surface region
        tempList=sum(ismember(S_6f,peakRegion(rid).faces),2)>=2;
        % Insert newly found triangles in surface region
        peakRegion(rid).faces=S_6f(tempList,:);
    end
    rid=rid+1;
    % Delete all surface triangles from seed list
    S_6f(tempList,:)=[];
end

for i=1:length(peakRegion)
    rFaces=sum(ismember(pelvis.faces,peakRegion(i).faces),2) == 3;
    subMesh = removeMeshFaces(pelvis, ~rFaces);
    % extract data from region
    peakRegion(i).vertices = subMesh.vertices;
    peakRegion(i).area =  meshSurfaceArea(subMesh.vertices,subMesh.faces);
    peakRegion(i).centroids = mean(subMesh.vertices);
end

% Ravine region
% Counter for number of surface regions
rid=1;
ravineRegion = struct('faces',{});
% Create surface regions while seed list is not empty
while ~isempty(S_1f)
    % Counter for number of triangles in surface region
    counter = 1;
    % Use the first triangle of the seed list as reference
    % Get all edge connected faces of the reference triangle
    tempList=sum(ismember(S_1f,S_1f(1,:)),2)>=2;
    % Create surface region
    % Insert reference triangle and edge connected faces
    ravineRegion(rid).faces=S_1f(tempList,:);
    % Recursive search for edge connected triangles of neighboring faces
    % while the number of triangles in surface region is bigger than
    % counter variable (number in previous step), repeat
    while sum(tempList) > counter
        counter = sum(tempList);
        % Get all edge connected faces of triangles in surface region
        tempList=sum(ismember(S_1f,ravineRegion(rid).faces),2)>=2;
        % Insert newly found triangles in surface region
        ravineRegion(rid).faces=S_1f(tempList,:);
    end
    rid=rid+1;
    % Delete all surface triangles from seed list
    S_1f(tempList,:)=[];
end

for i=1:length(ravineRegion)
    rFaces=sum(ismember(pelvis.faces,ravineRegion(i).faces),2) == 3;
    subMesh = removeMeshFaces(pelvis, ~rFaces);
    % extract data from region
    ravineRegion(i).vertices = subMesh.vertices;
    ravineRegion(i).area =  meshSurfaceArea(subMesh.vertices,subMesh.faces);
    ravineRegion(i).centroids = mean(subMesh.vertices);
end

% Ridge region
% Counter for number of surface regions
rid=1;
ridgeRegion = struct('faces',{});
% Create surface regions while seed list is not empty
while ~isempty(S_4f)
    % Counter for number of triangles in surface region
    counter = 1;
    % Use the first triangle of the seed list as reference
    % Get all edge connected faces of the reference triangle
    tempList=sum(ismember(S_4f,S_4f(1,:)),2)>=2;
    % Create surface region
    % Insert reference triangle and edge connected faces
    ridgeRegion(rid).faces=S_4f(tempList,:);
    % Recursive search for edge connected triangles of neighboring faces
    % while the number of triangles in surface region is bigger than
    % counter variable (number in previous step), repeat
    while sum(tempList) > counter
        counter = sum(tempList);
        % Get all edge connected faces of triangles in surface region
        tempList=sum(ismember(S_4f,ridgeRegion(rid).faces),2)>=2;
        % Insert newly found triangles in surface region
        ridgeRegion(rid).faces=S_4f(tempList,:);
    end
    rid=rid+1;
    % Delete all surface triangles from seed list
    S_4f(tempList,:)=[];
end

for i=1:length(ridgeRegion)
    rFaces=sum(ismember(pelvis.faces,ridgeRegion(i).faces),2) == 3;
    subMesh = removeMeshFaces(pelvis, ~rFaces);
    % extract data from region
    ridgeRegion(i).vertices = subMesh.vertices;
    ridgeRegion(i).area =  meshSurfaceArea(subMesh.vertices,subMesh.faces);
    ridgeRegion(i).centroids = mean(subMesh.vertices);
end

% Flat region
% Counter for number of surface regions
rid=1;
flatRegion = struct('faces',{});
% Create surface regions while seed list is not empty
while ~isempty(S_2f)
    % Counter for number of triangles in surface region
    counter = 1;
    % Use the first triangle of the seed list as reference
    % Get all edge connected faces of the reference triangle
    tempList=sum(ismember(S_2f,S_4f(1,:)),2)>=2;
    % Create surface region
    % Insert reference triangle and edge connected faces
    flatRegion(rid).faces=S_2f(tempList,:);
    % Recursive search for edge connected triangles of neighboring faces
    % while the number of triangles in surface region is bigger than
    % counter variable (number in previous step), repeat
    while sum(tempList) > counter
        counter = sum(tempList);
        % Get all edge connected faces of triangles in surface region
        tempList=sum(ismember(S_2f,flatRegion(rid).faces),2)>=2;
        % Insert newly found triangles in surface region
        flatRegion(rid).faces=S_2f(tempList,:);
    end
    rid=rid+1;
    % Delete all surface triangles from seed list
    S_4f(tempList,:)=[];
end

for i=1:length(flatRegion)
    rFaces=sum(ismember(pelvis.faces,flatRegion(i).faces),2) == 3;
    subMesh = removeMeshFaces(pelvis, ~rFaces);
    % extract data from region
    flatRegion(i).vertices = subMesh.vertices;
    flatRegion(i).area =  meshSurfaceArea(subMesh.vertices,subMesh.faces);
    flatRegion(i).centroids = mean(subMesh.vertices);
end

%Create one struct P containing all surface shapes with the fields:
% P.type: surface shape of region, P.faces: faces of region, P.vertices: vertices of
% region, P.area: area of surface region. P.centroids: centroids of surface
% region
P = struct('type',{});
for i=1:length(pitRegion)
    P(i).type = 'Pit';
    P(i).faces = pitRegion(i).faces;
    P(i).vertices = pitRegion(i).vertices;
    P(i).area = pitRegion(i).area;
    P(i).centroids = pitRegion(i).centroids;
end
rLength = length(P);
for i=1:length(ravineRegion)
    P(i + rLength).type = 'Ravine';
    P(i + rLength).faces = ravineRegion(i).faces;
    P(i + rLength).vertices = ravineRegion(i).vertices;
    P(i + rLength).area = ravineRegion(i).area;
    P(i + rLength).centroids = ravineRegion(i).centroids;
end
rLength = length(P);
for i=1:length(flatRegion)
    P(i + rLength).type = 'Flat';
    P(i + rLength).faces = flatRegion(i).faces;
    P(i + rLength).vertices = flatRegion(i).vertices;
    P(i + rLength).area = flatRegion(i).area;
    P(i + rLength).centroids = flatRegion(i).centroids;
end
rLength = length(P);
for i=1:length(ridgeRegion)
    P(i + rLength).type = 'Ridge';
    P(i + rLength).faces = ridgeRegion(i).faces;
    P(i + rLength).vertices = ridgeRegion(i).vertices;
    P(i + rLength).area = ridgeRegion(i).area;
    P(i + rLength).centroids = ridgeRegion(i).centroids;
end
rLength = length(P);
for i=1:length(peakRegion)
    P(i + rLength).type = 'Peak';
    P(i + rLength).faces = peakRegion(i).faces;
    P(i + rLength).vertices = peakRegion(i).vertices;
    P(i + rLength).area = peakRegion(i).area;
    P(i + rLength).centroids = peakRegion(i).centroids;
end

% Filter region according to surface area
% Regions with area below threshold can be ignored
i=1;
while i <= length(P)
    if P(i).area < 10.0
        P(i)=[];
        i=i-1;
    end
    i=i+1;
end
clear i;

% Create spatial adjancy matrix from region
Q_R=zeros(length(P)+1,length(P)+1,3);
% Add PS as reference point: If the pelvis was transformed into the the  
% automatic coordiante system [Kai 2014], PS is the origin: [0 0 0]
PS = [0 0 0]; 
% Calculate deviations between PS and all centroids of surface regions
for i=1:length(P)
    Q_R(1,i+1,1) = PS(1) - abs(P(i).centroids(1));
    Q_R(i+1,1,1) = abs(P(i).centroids(1)) - PS(1);
    Q_R(1,i+1,2) = PS(2) - P(i).centroids(2);
    Q_R(i+1,1,2) = P(i).centroids(2) - PS(2);
    Q_R(1,i+1,3) = PS(3) - P(i).centroids(3);
    Q_R(i+1,1,3) = P(i).centroids(3) - PS(3);
end

% Calculate deviations between all centroids of surface regions
for i=1:length(P)
    for j=1:length(P)
        Q_R(i+1,j+1,1) = abs(P(i).centroids(1))-abs(P(j).centroids(1));
        Q_R(i+1,j+1,2) = P(i).centroids(2)-P(j).centroids(2);
        Q_R(i+1,j+1,3) = P(i).centroids(3)-P(j).centroids(3);
    end
end

% Transform adjacency matrix to matrix coded according to anatomical terms
% of location. L: Lateral, M: Medial, A: Anterior, P: Posterior, S:
% Superior, I: Inferior, -: identical location in one plane (thresholded)
Q_A= cell(length(Q_R),length(Q_R));
for i=1:length(Q_A)
    for j=1:length(Q_A)
        Q_A{i,j} = '';
    end
end

for i=1:length(Q_R)
    for j=i+1:length(Q_R)-1
        
        if abs(Q_R(i,j,1)) < 5
            a ='-';
        elseif Q_R(i,j,1) < 0
            a ='L';
        elseif Q_R(i,j,1) > 0
            a ='M';
        end
        
        if abs(Q_R(i,j,2)) < 5
            b =':-';
        elseif Q_R(i,j,2) < 0
            b =':A';
        elseif Q_R(i,j,2) > 0
            b =':P';
        end
        
        if abs(Q_R(i,j,3)) < 5
            c = ':-';
        elseif Q_R(i,j,3) < 0
            c =':S';
        elseif Q_R(i,j,3) > 0
            c =':I';
        end
        Q_A{i,j} = [a b c];
    end
    Q_A{i,i} = '-:-:-';
end
clear a b c;

% Find landmarks in surface regions
% Load spatial adjacency matrix R_A
load(strcat(mode,'.mat'),'R_A');

% Initial search for reference point PS
% All candidates (surface region) fulfilling spatial relationship found in
% R_A are potential landmarks
C_PS = searchAdjM(R_A(4,5:end),Q_A(1,:));
possibleCandidate = 1;
for i=1:length(C_PS)
    % Create list containing candidates for each landmark
    CL(i).indices = C_PS(i).indices;
end

% For each landmark subject to detection (amount of CL found in R_A),
% repeat
i=1;
while i<length(CL)
    d_CL = [];
    for j=1:length(CL)-i
        d_t(j).indices =[];
    end
    
    %Determine landmark type and area
    lmType = R_A{i+4,2};
    lmArea = str2double(R_A{i+4,3});
    
    % Iterate over all candidates
    for k=1:length(CL(i).indices)
        
        %Determine candidates Type and Area
        type = P(CL(i).indices(k)-1).type;
        lmFaces=sum(ismember(pelvis.faces,P(CL(i).indices(k)-1).faces),2) == 3;
        subMesh = removeMeshFaces(pelvis, ~lmFaces);
        area = meshSurfaceArea(subMesh.vertices,subMesh.faces);
        
        %Check if type and area match landmark
        %         if ~(area >= lmArea*0.25 & area <= lmArea*1.75 & strcmp(type,lmType))
        if ~(strcmp(type,lmType))
            % if not, continue and skip candidate
            continue;
        end
        
        possibleCandidate = 1;
        % Determine all regions fulfilling spatial relationship to
        % candidate
        d = searchAdjM(R_A(i+4,i+5:end),Q_A(CL(i).indices(k),:));
        d_h = searchAdjM(R_A(i+5:end,i+4).',Q_A(:,CL(i).indices(k)).');
        for j=1:length(d)
            % create list with potential regions in spatial relation to
            % candidate
            d(j).indices = [d(j).indices d_h(j).indices];
            d(j).indices = unique(d(j).indices);
        end
        % If no potential regions for every candidate is found, candidate
        % can't be the desired landmark. Skip Candidate.
        for j=1:length(d)-1
            if isempty(intersect(CL(j+i).indices, d(j).indices))
                possibleCandidate = 0;
                break;
            end
        end
        % If potential regions are found for every candidate, create new
        % list including candidates for all landmarks. List is filled
        % whenever candidate fulfills requirements
        if possibleCandidate == 1
            d_CL = [d_CL CL(i).indices(k)];
            for j=1:length(d)
                d_t(j).indices = [d_t(j).indices d(j).indices];
                d_t(j).indices = unique(d_t(j).indices);
            end
        end
    end
    % Update initial list of candidates. Candidates found in old list, and
    % are not part of new list, are deleted.
    CL(i).indices = d_CL;
    for j=1:length(d_t)
        CL(j+i).indices = intersect(CL(j+i).indices, d_t(j).indices);
    end
    
    i=i+1;
    clear d_t d_CL;
end
clear d d_h i;




% Repeat search for landmarks in reverse order
possibleCandidate = 1;
i=length(CL);
% For each landmark subject to detection (amount of CL found in R_A),
% repeat
while i>=2
    d_CL = [];
    for j=i-1:-1:1
        d_t(j).indices =[];
    end
    
    %Determine landmark type and area
    lmType = R_A{i+4,2};
    lmArea = str2double(R_A{i+4,3});
    
    
    for k=1:length(CL(i).indices)
        
        %Determine Candidate Type and Area
        type = P(CL(i).indices(k)-1).type;
        lmFaces=sum(ismember(pelvis.faces,P(CL(i).indices(k)-1).faces),2) == 3;
        subMesh = removeMeshFaces(pelvis, ~lmFaces);
        area = meshSurfaceArea(subMesh.vertices,subMesh.faces);
        
        %Check if type and area match landmark
        %         if ~(area >= lmArea*0.25 & area <= lmArea*1.75 & strcmp(type,lmType))
        if ~(strcmp(type,lmType))
            continue;
        end
        
        possibleCandidate = 1;
        % Determine all regions fulfilling spatial relationship to
        % candidate
        d = searchAdjM(R_A(i+4,5:i+3),Q_A(CL(i).indices(k),:));
        d_h = searchAdjM(R_A(5:i+3,i+4).',Q_A(:,CL(i).indices(k)).');
        for j=1:length(d)
            % create list with potential regions in spatial relation to
            % candidate
            d(j).indices = [d(j).indices d_h(j).indices];
            d(j).indices = unique(d(j).indices);
        end
        % If no potential regions for every candidate is found, candidate
        % can't be the desired landmark. Skip Candidate.
        for j=1:length(d)
            if isempty(intersect(CL(j).indices, d(j).indices))
                possibleCandidate = 0;
                break;
            end
        end
        % If potential regions are found for every candidate, create new
        % list including candidates for all landmarks. List is filled
        % whenever candidate fulfills requirements
        if possibleCandidate == 1
            d_CL = [d_CL CL(i).indices(k)];
            for j=1:length(d)
                d_t(j).indices = [d_t(j).indices d(j).indices];
                d_t(j).indices = unique(d_t(j).indices);
            end
        end
    end
    % Update initial list of candidates. Candidates found in old list, and
    % are not part of new list, are deleted.
    CL(i).indices = d_CL;
    for j=1:length(d_t)
        CL(j).indices = intersect(CL(j).indices, d_t(j).indices);
    end
    clear d_t d_CL;
    
    i=i-1;
end
clear d d_h i;

% Store Candidates in output
for i=1:length(CL)
    lmArea = str2double(R_A{i+4,3});
    lmType = R_A{i+4,2};
    lmName = R_A{i+4,1};
    Landmarks(i).subMesh =[];
    Landmarks(i).centroids =[];
    Landmarks(i).area =[];
    Landmarks(i).id =[];
    Landmarks(i).name = lmName;
    
    for j=1:length(CL(i).indices)
        lmFaces=sum(ismember(pelvis.faces,P(CL(i).indices(j)-1).faces),2) == 3;
        subMesh = removeMeshFaces(pelvis, ~lmFaces);
        area = meshSurfaceArea(subMesh.vertices,subMesh.faces);
        type = P(CL(i).indices(j)-1).type;
        
        % Check if candidate type and area are matching landmark
        if strcmp(type,lmType)
            %         if area >= lmArea*0.7 & area <= lmArea*1.3 & strcmp(type,lmType)
            Landmarks(i).subMesh = [Landmarks(i).subMesh subMesh];
            Landmarks(i).centroids = [Landmarks(i).centroids; P(CL(i).indices(j)-1).centroids];
            Landmarks(i).area = [Landmarks(i).area; area];
            Landmarks(i).id = [Landmarks(i).id CL(i).indices(j)];
            
        end
    end
end

% If more than 2 candidates are detected for one landmark, keep only those
% with the highest area and delete the rest.

% NEW: Group candidates in positive- and negative side first.
% If two landmarks are detected, choose the one with the highest area

if sym
    for i=1:length(Landmarks)
        if isempty(Landmarks(i).centroids)
            pos =[]; neg=[];
        else
            pos = Landmarks(i).area(sign(Landmarks(i).centroids(:,1))==1);
            neg = Landmarks(i).area(sign(Landmarks(i).centroids(:,1))==-1);
        end
        if ~isempty(pos)
            idxPos = find(Landmarks(i).area==max(pos));
            posArea =  Landmarks(i).area(idxPos);
            posId = Landmarks(i).id(idxPos);
            posCentroids = Landmarks(i).centroids(idxPos,:);
            posSubMesh = Landmarks(i).subMesh(idxPos);
        else
            posArea =  nan;
            posId = nan;
            posCentroids = nan(1,3);
            posSubMesh.vertices = [];
            posSubMesh.faces = [];
        end
        
        if ~isempty(neg)
            idxNeg = find(Landmarks(i).area==max(neg));
            negArea =  Landmarks(i).area(idxNeg);
            negId = Landmarks(i).id(idxNeg);
            negCentroids = Landmarks(i).centroids(idxNeg,:);
            negSubMesh = Landmarks(i).subMesh(idxNeg);
        else
            negArea =  nan;
            negId = nan;
            negCentroids = nan(1,3);
            negSubMesh.vertices = [];
            negSubMesh.faces = [];
        end
        
        Landmarks(i).area = [posArea; negArea];
        Landmarks(i).id = [posId; negId];
        Landmarks(i).centroids = [posCentroids; negCentroids];
        Landmarks(i).subMesh = [posSubMesh; negSubMesh];
    end
else
    for i=1:length(Landmarks)
        if length(Landmarks(i).area) > 2
            j=length(Landmarks(i).area);
            while j > 2
                idx = find(Landmarks(i).area==min(Landmarks(i).area));
                Landmarks(i).area(idx) = [];
                Landmarks(i).id(idx) = [];
                Landmarks(i).centroids(idx,:) = [];
                Landmarks(i).subMesh(idx) = [];
                j=j-1;
            end
        end
    end
end

%visualization
if visu == true
    patchProps.EdgeColor = 'none';
    patchProps.FaceColor = [0.75 0.75 0.75];
    % patchProps.FaceAlpha = 0.5;
    patchProps.EdgeLighting = 'gouraud';
    patchProps.FaceLighting = 'gouraud';
    
    pointProps.Linestyle = 'none';
    pointProps.Marker = 'o';
    pointProps.MarkerEdgeColor = 'k';
    pointProps.MarkerFaceColor = 'r';
    pointProps.Marker = 'd';
    
    for i=1:length(Landmarks)
        patchProps.FaceColor = rand(1,3);
        for j=1:length(Landmarks(i).subMesh)
            if ~isempty(Landmarks(i).subMesh(j).vertices)
            % patch(Landmarks(i).subMesh(j), patchProps);
            drawPoint3d(mean(Landmarks(i).subMesh(j).vertices), pointProps);
            text(mean(Landmarks(i).subMesh(j).vertices(:,1)), ...
                mean(Landmarks(i).subMesh(j).vertices(:,2)), ...
                mean(Landmarks(i).subMesh(j).vertices(:,3)), ...
                R_A(i+4,1),...
                'FontWeight','bold','FontSize',16,...
                'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
            end
        end
    end
    mouseControl3d
    medicalViewButtons('RAS')
end


end