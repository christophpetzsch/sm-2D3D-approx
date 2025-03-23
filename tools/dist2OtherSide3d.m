function [dist2other, dist2otherAll] =  dist2OtherSide3d(M, normalize)
%dist2OtherSide3d Computes distance to other side of triangle mesh
% will return infinity if normal at current vertex is not 
if ~exist('TriangleRayIntersection.m','file')
    p = mfilename('fullpath');
    filepath = fileparts(p);
    addpath(strcat(filepath, '/TriangleRayIntersection'))
    if ~exist('TriangleRayIntersection.m','file')
        error(['Please install Toolbox TriangleRayIntersection: '...
        'https://de.mathworks.com/matlabcentral/fileexchange/33073-triangle-ray-intersection'])
    end
end

normalize_ = false;
if nargin > 1
    normalize_ = normalize;
end

TR = triangulation(M.TRIV, M.VERT);
if ~isfield(M, 'Normal')
    M.Normal = vertexNormal(TR);
end

faceNormals = faceNormal(TR);

% weird way of dealing with triangles
vert1 = M.VERT(M.TRIV(:,1),:);
vert2 = M.VERT(M.TRIV(:,3),:);
vert3 = M.VERT(M.TRIV(:,2),:);
% origin of rays 
orig = M.VERT;
% direction of rays (negative normal)
dir = -M.Normal;

orig = orig + eps * dir;

% weird way of dealing with multiple rays => is slower than for-loop
% nPoints = length(orig);
% orig = repelem(orig, length(vert1), 1);
% dir = repelem(dir, length(vert1), 1);
% vert1 = repmat(vert1, nPoints, 1);
% vert2 = repmat(vert2, nPoints, 1);
% vert3 = repmat(vert3, nPoints, 1);

dist2other = zeros(size(M.VERT(:, 1)));
dist2otherAll = zeros(size(M.VERT(:, 1), 1), 100);
maxSize = 1;
for i = 1:length(orig)
    [triIdx, ~, ~, ~, xcoor] = TriangleRayIntersection (orig(i, :), dir(i, :), vert1, vert2, vert3,  'planeType', 'two sided');
    xcoor = xcoor(~isnan(xcoor(:, 1)), :);
    if isempty(xcoor)
        dist2other(i) = inf;
        dist2otherAll(i, :) = inf;
        continue;
    end
    normalsOfIntersectingFaces = faceNormals(triIdx, :);
    dirRep = repmat(dir(i, :), size(normalsOfIntersectingFaces, 1), 1);
    anglesWNormalsIntersecting = atan2(rowNorm(cross(normalsOfIntersectingFaces, dirRep, 2)),...
        dot(normalsOfIntersectingFaces, dirRep, 2));
    angleCondition = anglesWNormalsIntersecting < (pi/2);
    % angle conditions ensures that we only consider outwards facing
    % triangles (ray comes from inside of mesh)
    xcoor = xcoor(angleCondition, :);
    
    % make sure we do not calc the dist to the current element
    distEachElem = rowNorm(M.VERT(i, :) - xcoor);
    xcoor = xcoor(distEachElem > 1e-10, :);
    if isempty(xcoor)
        dist2other(i) = inf;
        dist2otherAll(i, :) = inf;
        continue;
    end
    distEachElem = distEachElem(distEachElem > 1e-10, :);
    distEachElem = sort(distEachElem);
    dist2other(i) = distEachElem(1);
    dist2otherAll(i, 1:length(distEachElem)) = distEachElem;
    if length(distEachElem) > maxSize
        maxSize = length(distEachElem);
    end
end


% smoothing (averaging over neighbors)
% A = adjacency_matrix(M.TRIV);
% normaldist = (normaldist + sum(normaldist(A), 2)) ./ (sum(A, 2) + 1)

%fix infinities
isInfVec = isinf(dist2other);
idxInf = find(isInfVec);
pointsNotInf = M.VERT(~isInfVec, :);
dist2otherNotInf = dist2other(~isInfVec);
nnInfPoints = knnsearch(pointsNotInf,  M.VERT(idxInf, :), 'k', 2);
dist2other(idxInf) = 0.5 * ( dist2otherNotInf(nnInfPoints(:, 1)) + dist2otherNotInf(nnInfPoints(:, 2)) );
dist2otherAll(idxInf, 1) = dist2other(idxInf);


dist2otherAll = dist2otherAll(:, 1:maxSize);

% normalization
%normaldist = normaldist ./ max(normaldist);
if normalize_
median_ = median(dist2other(~isinf(dist2other)));
dist2other = dist2other ./ median_;
dist2otherAll = dist2otherAll ./ median_;
end

end

