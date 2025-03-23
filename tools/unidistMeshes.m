function [M, N] = unidistMeshes(M, N, meshResamplingFactor, edgeLengthMultiple2d, method)
% this is the important parameter to adust (keep ratio of the elements).
% This reduces the mesh resolution, if that is undesirable, we need to
% upsample beforehand
meshResamplingFactorInt = 0.6;
edgeLengthMultiple2dInt = 1;
methodInt = 'naive';
if nargin > 2
    meshResamplingFactorInt = meshResamplingFactor;
end
if nargin > 3
    edgeLengthMultiple2dInt = edgeLengthMultiple2d;
end
if nargin > 4
    methodInt = method;
end
if ~exist('meshresample')
    p = mfilename('fullpath');
    [filepath,~,~] = fileparts(p);
    %addpath([filepath,'/../external/iso2mesh']);
end

%% 3d
if meshResamplingFactorInt == 1
    VERTnew = M.VERT; 
    TRIVnew = M.TRIV;
    NNintoOriginal = (1:length(M.VERT))';
else
    [VERTnew, TRIVnew, ~, NNintoOriginal] = decimate_libigl(M.VERT, M.TRIV, meshResamplingFactorInt, 'Method', methodInt);
end
% [VERTnew, TRIVnew] = meshresample(M.VERT, M.TRIV, meshResamplingFactorInt);
% [VERTnew, TRIVnew] = meshcheckrepair(VERTnew, TRIVnew);
medEdgeLenghts = median(getEdgeLengths(VERTnew, TRIVnew));
M.VERT = VERTnew;
M.TRIV = TRIVnew;
M.n = length(M.VERT);
M.m = length(M.TRIV);
M = transferQuantities3d(M, NNintoOriginal);


%% 2d
BoundaryPoints = N.VERT(logical(N.boundary), :);
numBoundaryPointsOriginal = sum(logical(N.boundary));
ORIGINAL_SEGMENT = [];
if isfield(N, 'segmentation')
    ORIGINAL_SEGMENT = N.segmentation;
end
BoundaryPoints = [BoundaryPoints; BoundaryPoints(1, :)];
InteriorPoints = N.VERT(sum(logical(N.boundary))+1:end, :);
NewBoundaryPoints = uniformEdgeLengthPolygon(BoundaryPoints, edgeLengthMultiple2dInt * medEdgeLenghts);

NNintoOriginal = knnsearch(BoundaryPoints(:, 1:2), NewBoundaryPoints);

% triangulate interior
ps = polyshape(NewBoundaryPoints(:, 1), NewBoundaryPoints(:, 2), 'Simplify', false);
EdgeConstraint = [(1:length(NewBoundaryPoints))' [(2:length(NewBoundaryPoints))'; 1]];
% sample more points
%mini = min(BoundaryPoints(:, 1:2));
%maxi = max(BoundaryPoints(:, 1:2));
%[X,Y] = meshgrid(mini(1):medEdgeLenghts:maxi(1), mini(2):medEdgeLenghts:maxi(2));
NewInteriorPoints = InteriorPoints(:, 1:2); %[X(:),Y(:)];
TF = isinterior(ps, NewInteriorPoints);
NewInteriorPoints = NewInteriorPoints(TF, :);

DT = delaunayTriangulation([NewBoundaryPoints; NewInteriorPoints], EdgeConstraint);
% remove outside triangles
C = incenter(DT);
TF = isinterior(ps, C);
InsideTriangles = DT.ConnectivityList(TF, :);
if numBoundaryPointsOriginal < length(NewBoundaryPoints)
    N = [];
end
N.VERT = [DT.Points, zeros(length(DT.Points), 1)];
N.boundary = (1:length(NewBoundaryPoints))';
N.TRIV = InsideTriangles;
N.n = length(N.VERT);
N.m = length(N.TRIV);
if numBoundaryPointsOriginal < length(NewBoundaryPoints)
    % remove zero area triangles
    N.S_tri = calc_tri_areas(N);
    tooSmallAreaIdx = N.S_tri < 1e-12;
    if any(tooSmallAreaIdx)
        N.TRIV = N.TRIV(~tooSmallAreaIdx, :);
        N.m = length(N.TRIV);
    end
    
    N = computeShapeStruct(N, 1, 25, 100, 6, true, sum(logical(N.boundary)));
    if ~isempty(ORIGINAL_SEGMENT)
        NNN = knnsearch(BoundaryPoints(: ,1:2), NewBoundaryPoints);
        N.segmentation = ORIGINAL_SEGMENT(NNN);
    end
else  
    N = transferQuantities2d(N, NNintoOriginal);
end

end




function elen = getEdgeLengths(V, F)
    TR = triangulation(F, V);
    E = edges(TR);
    elen = sqrt(sum( (V(E(:, 1), : ) - V(E(:, 2), : ) ).^2, 2));
end

function uniformPoly = uniformEdgeLengthPolygon(poly, edgeLength)
	% code modified from http://de.mathworks.com/matlabcentral/newsreader/view_thread/270021
	% here is a simple polygon, a triangle. Note that
	% I've wrapped the ends, so that the last point is
	% also the first point. This is necessary.
	px = poly(:,1);
	py = poly(:,2);

	% t is the cumulative arclength along the edges of the polygon.
	t = cumsum(sqrt([0,diff(px([1:end 1])').^2] + [0,diff(py([1:end 1])').^2]));

	% The total distance around the polygon is t(end)
	tmax = t(end);
    
    % number of points is tmax / desired edge length
    n = tmax / edgeLength;

	% create a piecewise linear spline for each of px and py,
	% as a function of the cumulative chordwise arclength.
	splx = mkpp(t,[diff(px([1:end, 1]))./diff(t'),px]);
	sply = mkpp(t,[diff(py([1:end, 1]))./diff(t'),py]);

	% now interpolate the polygon splines, splx and sply.
	% Nt is the number of points to generate around the
	% polygon. The first and last points should be replicates
	% at least to within floating point tresh.)
	tint = linspace(0,tmax,n);

	qx = ppval(splx,tint)';
	qy = ppval(sply,tint)';

	uniformPoly = [qx qy];
	% remove duplicate start/end point
	uniformPoly(end,:) = [];
% 	figure
% 	% plot the polygon itself, as well as the generated points.
% 	plot(px,py,'k-v',qx,qy,'ro')
% 	grid on 
end


function N = transferQuantities2d(N, NNintoOriginal)
    if isfield(N, 'wks')
    N.wks = N.wks(NNintoOriginal, :);
    end
    if isfield(N, 'hks')
    N.hks = N.hks(NNintoOriginal, :);
    end
    if isfield(N, 'gn_wks')
    N.gn_wks = N.gn_wks(NNintoOriginal, :);
    end
    if isfield(N, 'rn_wks')
    N.rn_wks = N.rn_wks(NNintoOriginal, :);
    end
    if isfield(N, 'shks')
    N.shks = N.shks(NNintoOriginal, :);
    end
    if isfield(N, 'gn_shks')
    N.gn_shks = N.gn_shks(NNintoOriginal, :);
    end
    if isfield(N, 'rn_shks')
    N.rn_shks = N.rn_shks(NNintoOriginal, :);
    end
    if isfield(N, 'Normal')
    N.Normal = N.Normal(NNintoOriginal, :);
    end
    if isfield(N, 'Dist2Other')
    N.Dist2Other = N.Dist2Other(NNintoOriginal, :);
    end
    if isfield(N, 'segmentation')
    N.segmentation = N.segmentation(NNintoOriginal);
    end
    if isfield(N, 'geoDistToCenter')
    N.geoDistToCenter = N.geoDistToCenter(NNintoOriginal, :);
    end
    if isfield(N, 'BoundCurv')
    N.BoundCurv = N.BoundCurv(NNintoOriginal, :);
    end
end


function M = transferQuantities3d(M, NNintoOriginal)
    if isfield(M, 'wks')
    M.wks = M.wks(NNintoOriginal, :);
    end
    if isfield(M, 'hks')
    M.hks = M.hks(NNintoOriginal, :);
    end
    if isfield(M, 'gn_wks')
    M.gn_wks = M.gn_wks(NNintoOriginal, :);
    end
    if isfield(M, 'rn_wks')
    M.rn_wks = M.rn_wks(NNintoOriginal, :);
    end
    if isfield(M, 'shks')
    M.shks = M.shks(NNintoOriginal, :);
    end
    if isfield(M, 'gn_shks')
    M.gn_shks = M.gn_shks(NNintoOriginal, :);
    end
    if isfield(M, 'rn_shks')
    M.rn_shks = M.rn_shks(NNintoOriginal, :);
    end
    if isfield(M, 'geoDistToCenter')
    M.geoDistToCenter = M.geoDistToCenter(NNintoOriginal, :);
    end
    if isfield(M, 'Umin')
    M.Umin = M.Umin(NNintoOriginal, :);
    end
    if isfield(M, 'Umax')
    M.Umax = M.Umax(NNintoOriginal, :);
    end
    if isfield(M, 'Cmin')
    M.Cmin = M.Cmin(NNintoOriginal, :);
    end
    if isfield(M, 'Cmax')
    M.Cmax = M.Cmax(NNintoOriginal, :);
    end
    if isfield(M, 'Cmean')
    M.Cmean = M.Cmean(NNintoOriginal, :);
    end
    if isfield(M, 'Cgauss')
    M.Cgauss = M.Cgauss(NNintoOriginal, :);
    end
    if isfield(M, 'Dist2OtherAll')
    M.Dist2OtherAll = M.Dist2OtherAll(NNintoOriginal, :);
    end
    if isfield(M, 'Normal')
    M.Normal = M.Normal(NNintoOriginal, :);
    end
    if isfield(M, 'Dist2Other')
    M.Dist2Other = M.Dist2Other(NNintoOriginal, :);
    end
    if isfield(M, 'segmentation')
    M.segmentation = M.segmentation(NNintoOriginal);
    end
end
