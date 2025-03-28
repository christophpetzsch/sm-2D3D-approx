clear all 
addpath('data', genpath('tools'), genpath('external'));

%% Data Loading
load data/faust/faust_041.mat
load data/faust/faust_sketch4.mat


%% Data Preprocessing (downsampling and computation of features)
Norig = N; 
Morig = M;


% downsampling settings for our experiments
reductionFactor2d = 1.0;
reductionFactor3d = 0.5;

varWKS = 6;
numEV = 25;
sizeWKS = 100;

[N, NNN] = computeShapeStruct(N, 1, numEV, sizeWKS, varWKS, true, sum(logical(N.boundary)));
N.segmentation = Norig.segmentation(NNN);
[M, MNN] = computeShapeStruct(M, 1, numEV, sizeWKS, varWKS, false, 0);
M.segmentation = Morig.segmentation(MNN);
% reduction
[M, N] = unidistMeshes(M, N, reductionFactor3d, reductionFactor2d);



%% Settings
graphType = 2; % 0 -> Laehner et al (Product Graph), 1 -> Roetzer et al (Conjugate Product Graph), 2 -> Petzsch et al (Blended Product Graph)
retrieval_normalization = false;
landmarks = false;


%% Compute results
tic
if(~landmarks)
    matching = sm_2d3d_dijkstra(M, N, graphType, retrieval_normalization);
    %% Visualize result
    boundary = find(N.boundary ~= 0);
    boundary = [boundary; boundary(1)];

    plot_path(M, N.VERT(boundary, 1:2), matching(:, 1), matching(:, 2));
else
    matching = sm_2d3d_dijkstra_landmarks(M,N, graphType, retrieval_normalization);
end
toc



