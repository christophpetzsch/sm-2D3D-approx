function dist2OtherSide = dist2OtherSide2d(N, normalize)
normalize_ = false;
if nargin > 1
    normalize_ = normalize;
end

% find normals 
if ~isfield(N, 'Normal')
    N.Normal = normals2D(N.VERT(logical(N.boundary), :));
end

% find diagonal of bounding box 
MINS = min(N.VERT(logical(N.boundary), 1:2));
MAXS = max(N.VERT(logical(N.boundary), 1:2));
ShapeDiag = norm(MINS-MAXS);

% find with help of normals distance of other side of shape 
dist2OtherSide = zeros(sum(logical(N.boundary)), 1);
warning off
poly = polyshape(N.VERT(logical(N.boundary), 1:2));
warning on 
for i = 1: sum(logical(N.boundary))
    idx = i; 
    if i == 230
        disp('')
    end
    line = [ N.VERT(idx, 1:2); N.VERT(idx, 1:2) - ShapeDiag * N.Normal(idx, 1:2)];
    in = intersect(poly,line);
    if isempty(in) 
        dist2OtherSide(i) = inf;
    else
        [~, idxMin] = mink(rowNorm(in - N.VERT(idx, 1:2)), 2);
        % first two in points are ins closest to current vert
        dist2OtherSide(i) = norm(N.VERT(idx, 1:2) - in(idxMin(end),:)); 
    end
end


% do some smoothing peutetre
%dist2OtherSide = smoothdata(dist2OtherSide,'gaussian',10);

% fix infinities
isInfVec = isinf(dist2OtherSide);
idxInf = find(isInfVec);
BoundaryPoints = N.VERT(logical(N.boundary), :);
pointsNotInf = BoundaryPoints(~isInfVec, :);
nnInfPoints = knnsearch(pointsNotInf,  N.VERT(idxInf, :), 'k', 2);
dist2OtherSide(idxInf) = 0.5 * ( dist2OtherSide(nnInfPoints(:, 1)) + ...
    dist2OtherSide(nnInfPoints(:, 2)) );

% normalization
%dist2OtherSide = dist2OtherSide ./ max(dist2OtherSide);

if normalize_
median_ = median(dist2OtherSide(~isinf(dist2OtherSide)));
dist2OtherSide = dist2OtherSide ./ median_;
end

end

