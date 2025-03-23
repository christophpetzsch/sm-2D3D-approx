function [vertexMatching, e] = sm_2d3d_dijkstra_landmarks(M, N, graphType, retrieval)
%sm_2d3d_dijkstra wrapper for generated mex functions to solve 2D3D matching problem
%
% Note: you can use conjugateGraph = false to compute matching without the conjugate graph.
%       For that you have to download respective source files and put them
%       into dijkstra folder.

%% handle inputs
%if nargin < 3
%    graphType = 1;
%end

%% generate cpp code if not existent
test = dir('dijkstraCG.mex*');
wrnMsg = 'Building relevant source files. Optimisation time might be wrong.';

if graphType == 0
    test = dir('dijkstra.mex*');
    if isempty(test)
        warning(wrnMsg)
        buildDijkstra
    end
end

if graphType == 1
    test = dir('dijkstraCG.mex*');
    if isempty(test)
        warning(wrnMsg)
        buildDijkstraCG
    end
end

if graphType == 2
    test = dir('dijkstraCGapprox.mex*');
    if isempty(test)
        warning(wrnMsg)
        buildDijkstraCGapprox
    end
end




%% Precompute Rotations

test = dir('rotations.mex*')
if isempty(test)
    warning(wrnMsg)
    buildRotations
end

boundary = find(N.boundary ~= 0);
boundary = [boundary; boundary(1)];
if ~ (graphType == 0)
    % boundary = boundary boundary(1:2)
    boundary = [boundary; boundary(2)];
end


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



if ~ (graphType == 0)
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
    fprintf('Precomputing rotations...')
    quaternions = rotations(descr_M', M.TRIV'-1, descr_N');
    disp(' done.')
end



boundary = find(N.boundary ~= 0);
boundary = [boundary; boundary(1)];
plot_nopath_landmarks(M, N.VERT(boundary, 1:2))
uiwait
new_2d_landmarks = evalin('base', 'new_2d_landmarks');
new_3d_landmarks = evalin('base', 'new_3d_landmarks');
landmarks_2d = int32(new_2d_landmarks);
landmarks_3d = int32(new_3d_landmarks);



path = [];
corr = [];

while(~isempty(new_2d_landmarks))

    %% duplicate first layer
    boundary = find(N.boundary ~= 0);
    


    [i, j] = find(ismember(abs(landmarks_2d' - landmarks_2d), [1,length(descr_N')-3]));
    pairs = i < j;
    indices = [i(pairs), j(pairs)];
    if(~isempty(indices))
        if((landmarks_2d(i(1))<landmarks_2d(j(1)) && (landmarks_2d(i(1)) ~= 0)) || ((landmarks_2d(i(1))==length(descr_N')-3) && (landmarks_2d(j(1))==0)))
            prev_shift = landmarks_2d(i(1));
        else
            prev_shift = landmarks_2d(j(1));
        end

        disp(length(descr_N')-2)
        disp(prev_shift)
        disp(landmarks_2d)
        landmarks_2d = mod(landmarks_2d - prev_shift, length(descr_N')-2);
        disp(landmarks_2d)

        [landmarks_2d,sortIdx] = sort(landmarks_2d);
        landmarks_3d = landmarks_3d(sortIdx);
        boundary = circshift(boundary, -prev_shift);
        landmarks_2d(1) = prev_shift;
        disp(landmarks_2d)
    else
        [landmarks_2d,sortIdx] = sort(landmarks_2d);
        landmarks_3d = landmarks_3d(sortIdx);
        boundary = circshift(boundary, -landmarks_2d(1));
        prev_shift = landmarks_2d(1);
        landmarks_2d = landmarks_2d - landmarks_2d(1);
        landmarks_2d(1) = prev_shift;
    end

    landmarks_3d = landmarks_3d - 1;


    boundary = [boundary; boundary(1)];
    if ~ (graphType == 0)
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
    
    if ~ (graphType == 0)
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
        if graphType == 1
            [e, path, corr] = dijkstraCG(descr_M', M.TRIV'-1, descr_N', landmarks_2d, landmarks_3d, quaternions);
        else 
            [e, path, corr] = dijkstraCGapprox(descr_M', M.TRIV'-1, descr_N', landmarks_2d, landmarks_3d, quaternions);
        end
    else
        costMode = 6
        if costMode == 6
            % make sure segementation is not used
            segmentM = zeros(size(M.VERT, 1), 1);
            segmentN = zeros(size(N.VERT, 1), 1);
            descr_N = [descr_N(:, 1:3) N.Dist2Other(boundary) segmentN(boundary)];
            descr_M = [descr_M(:, 1:3) M.Dist2Other segmentM];
        elseif costMode > 2
            warning('Curvature Cost not supported for non line-graph calls');
        end
        [e, path, corr] = dijkstra(descr_M', M.TRIV'-1, descr_N', landmarks_2d, landmarks_3d);
    end
    
    fprintf('Final energy: %.4f\n', e);
    
    % path and corr come out as double and 0-index based
    path = uint64(path) + 1;
    corr = uint64(corr) + 1;
    

    
    if ~(isempty(landmarks_2d))
        corr(end) = [];
        path(end) = [];
        corr = mod(corr + uint64(prev_shift), length(descr_N')-2)+1;
        last_one_index = find(corr == 1, 1, 'first');
        corr = [corr(last_one_index:end); corr(1:last_one_index-1); corr(last_one_index)];
        path = [path(last_one_index:end); path(1:last_one_index-1); path(last_one_index)];

    end
    

    landmarks_2d(1) = 0;
    landmarks_2d = mod(landmarks_2d + prev_shift, length(descr_N')-2);
    landmarks_3d = landmarks_3d + 1;
    
    boundary = find(N.boundary ~= 0);
    boundary = [boundary; boundary(1)];
    new_2d_landmarks = [];
    new_3d_landmarks = [];
    assignin('base', 'new_2d_landmarks', new_2d_landmarks);
    assignin('base', 'new_3d_landmarks', new_3d_landmarks);    

    plot_path_landmarks(M, N.VERT(boundary, 1:2), path, corr, landmarks_2d, landmarks_3d);
    uiwait
    new_2d_landmarks = evalin('base', 'new_2d_landmarks');
    new_3d_landmarks = evalin('base', 'new_3d_landmarks');
    landmarks_2d = [landmarks_2d, int32(new_2d_landmarks)];
    landmarks_3d = [landmarks_3d, int32(new_3d_landmarks)];
end

if(retrieval)
    path_length = 0;
    for i = 1:length(path)-1
        point1 = [M.VERT(path(i),1), M.VERT(path(i),2), M.VERT(path(i),3),N.VERT(corr(i),1),N.VERT(corr(i),2)];
        point2 = [M.VERT(path(i+1),1), M.VERT(path(i+1),2), M.VERT(path(i+1),3),N.VERT(corr(i+1),1),N.VERT(corr(i+1),2)];

        distance = pdist([point1;point2],'euclidean');
        path_length = path_length + distance;
    end

    point1 = [M.VERT(path(length(path)),1), M.VERT(path(length(path)),2), M.VERT(path(length(path)),3),N.VERT(corr(length(path)),1),N.VERT(corr(length(path)),2)];
    point2 = [M.VERT(path(1),1), M.VERT(path(1),2), M.VERT(path(1),3),N.VERT(corr(1),1),N.VERT(corr(1),2)];

    distance = pdist([point1;point2],'euclidean');
    path_length = path_length + distance;

    e = e/path_length;

    fprintf('Energy after normalization: %.4f\n', e);
end

vertexMatching = [path, corr];

Result.M = M; 
Result.N = N; 
Result.matching = vertexMatching;
Result.energy = e;
