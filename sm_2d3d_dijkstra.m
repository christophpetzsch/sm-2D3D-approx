function [vertexMatching, e] = sm_2d3d_dijkstra(M, N, conjugateGraph)
%sm_2d3d_dijkstra wrapper for generated mex functions to solve 2D3D matching problem
%
% Note: you can use conjugateGraph = false to compute matching without the conjugate graph.
%       For that you have to download respective source files and put them
%       into dijkstra folder.

%% handle inputs
if nargin < 3
    conjugateGraph = true;
end

%% generate cpp code if not existent
test = dir('dijkstraCG.mex*');
wrnMsg = 'Building relevant source files. Optimistation time might be wrong.';
if isempty(test)
    warning(wrnMsg)
    buildDijkstraCG
end
if ~conjugateGraph
    test = dir('dijkstra.mex*');
    if ismepty(test)
        warning(wrnMsg)
        mex dijkstra/dijkstra.cpp dijkstra/MinHeap.cpp
    end
end



%% duplicate first layer
boundary = find(N.boundary ~= 0);
boundary = [boundary; boundary(1)];
if conjugateGraph
    % boundary = boundary boundary(1:2)
    boundary = [boundary; boundary(2)];
end


%% descriptors

% segmentation is not known => all to segment 0
segmentM = zeros(size(M.VERT, 1), 1);
segmentN = zeros(size(N.VERT, 1), 1);
if isfield(M, 'segmentation') && isfield(N, 'segmentation')
    % only makes sense when segmentation for both shapes known
    segmentM = M.segmentation;
    segmentN = N.segmentation;
end

% best working descriptors are taken from original paper
descr_M = [M.VERT M.rn_wks M.gn_shks segmentM];
descr_N = [N.VERT(boundary, 1:2) zeros(length(boundary),1)...
            N.rn_wks(boundary, :) N.gn_shks(boundary, :) ...
            segmentN(boundary)];

%% OPTIMIZE

if conjugateGraph
    % workaround
    if ~isfield(N, 'BoundCurv')
        N.BoundCurv = 0 * N.Dist2Other;
    end
    if ~isfield(M, 'Cmin')
        M.Cmin = 0 * M.Dist2Other;
        M.Cmax  = M.Cmin;
        M.Umin = zeros(length(M.Cmax), 3);
        M.Umax = M.Umin;
    end
    % end workaround
    med = median(N.Dist2Other); 
    if (med == 0)
        med = 1;
    end
    descr_N = [descr_N(:, 1:3) N.Normal(boundary, :) N.BoundCurv(boundary) N.Dist2Other(boundary)/med];
    descr_M = [descr_M(:, 1:3) M.Cmin M.Cmax M.Umin M.Umax M.Normal M.Dist2Other/med];
    [e, path, corr] = dijkstraCG(descr_M', M.TRIV'-1, descr_N');
else
    if costMode == 6
        % make sure segementation is not used
        segmentM = zeros(size(M.VERT, 1), 1);
        segmentN = zeros(size(N.VERT, 1), 1);
        descr_N = [descr_N(:, 1:3) N.Dist2Other(boundary) segmentN(boundary)];
        descr_M = [descr_M(:, 1:3) M.Dist2Other segmentM];
    elseif costMode > 2
        warning('Curvature Cost not supported for non line-graph calls');
    end
    [e, path, corr] = dijkstra(descr_M', M.TRIV'-1, descr_N');
end


fprintf('Final energy: %.4f\n', e);

% path and corr come out as double and 0-index based
path = uint64(path) + 1;
corr = uint64(corr) + 1;

vertexMatching = [path, corr];
