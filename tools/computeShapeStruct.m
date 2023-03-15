function [S, NNintoOriginal] = computeShapeStruct(S, reductionFactor, numEV, sizeWKS, varWKS, is2d, boundary2d, useMatlabReduction)
% V: |V| x 3 containing shape vertices
% F: |F| x 3 containing triangles (indices in V)
% reductionFactor \in [0, 1]: discretization percentage
% numWKS: number of eigenvalues of WKS descriptors
% sizeWKS: size of WKS descriptors
% numHKS: number of eigenvalues of HKS descriptors
% is2d: true if V, F are 2d shape and if V = [..., ..., zeros(|V|, 1)];
% boundary2d: integer so that V(1:boundary2d, :) are only the boundary vertices
useMRed = false;
useSmoothing = true;
useNormalization = false;
computeEigenFunctions = false;
computeCurvature = false;
clipDist2Other = true;

if nargin > 8
    useMRed = useMatlabReduction;
end
if reductionFactor > 1
    reductionFactor = 1;
end
if reductionFactor < 0
    reductionFactor = 0;
end

% S = [];
% S.TRIV = F;
% S.VERT = V;
V = S.VERT;
Vbound = [];
if is2d
    Vbound = S.VERT(logical(S.boundary), :);
end

S.n = size(S.VERT, 1);
S.m = size(S.TRIV, 1);
if ~isfield(S, 'Dist2Other') || ~isfield(S, 'Dist2OtherAll')
    if is2d
        S.boundary = (1:boundary2d)';
        S.Dist2Other = dist2OtherSide2d(S, false);

    else
        [S.Dist2Other, S.Dist2OtherAll] = dist2OtherSide3d(S, false);
    end
else
    % median normalize
    %S.Dist2Other = S.Dist2Other ./ median(S.Dist2Other);
end

if clipDist2Other
    medD2O = median(S.Dist2Other);
    idx = S.Dist2Other > 15 * medD2O;
    S.Dist2Other(idx) = 15 * medD2O;
end

%S.normals = normals3D(S);
if ~isfield(S, 'wks') || ~isfield(S, 'hks')  || ~isfield(S, 'evecs') || size(S.wks, 1) ~= size(S.VERT, 1)

    S.S_tri = calc_tri_areas(S);
    [S.evecs, S.evals, S.S, ~] = cotan_LB(S.VERT, S.TRIV, numEV);
    options = struct();
    options.wks_size = sizeWKS;
    options.wks_variance = varWKS;
    options.max_eigen = numEV;

    S = calc_descr(S, options);
end




if ~is2d
    TR = triangulation(S.TRIV, S.VERT);
    if computeCurvature
    if ~exist('principal_curvature')
        warning('principal_curvature not found. Please install gptoolbox https://github.com/alecjacobson/gptoolbox');
    end
    % this function is acutally shit
    % options.verb = 0;
    % options.curvature_smoothing = 1;
    % [S.Umin, S.Umax, S.Cmin, S.Cmax, S.Cmean, S.Cgauss, S.Normal] = compute_curvature(S.VERT, S.TRIV, options);
    % S.Umin = S.Umin';
    % S.Umax = S.Umax';

    % [S.Cmean, S.Cgauss, S.Umin, S.Umax, S.Cmin, S.Cmax] = patchCurvature(TR, false);
    % crashes for less then 6 vertices
    if size(S.VERT, 1) < 7
        F = S.TRIV; V = S.VERT;
        [VV,FF] = upsample(V,F,'Iterations',2);
        [Umax, Umin, Cmax, Cmin] = principal_curvature(VV, FF);
        S.Umax = Umax(1:length(V), :);
        S.Umin = Umin(1:length(V), :);
        S.Cmax = Cmax(1:length(V), :);
        S.Cmin = Cmin(1:length(V), :);
    else
        if useSmoothing
            % upsampling
            F = S.TRIV; V = S.VERT;
            [VV,FF,~] = upsample(V,F,'Iterations', 1);
            % smoothing
            U = laplacian_smooth(VV, FF, 'cotan', [], 0.1, 'implicit', VV, 100);
            % curv
            [Umax, Umin, Cmax, Cmin] = principal_curvature(U, FF);
            % translate to original shape (first vertices are original ones)
            S.Umax = Umax(1:length(V), :);
            S.Umin = Umin(1:length(V), :);
            S.Cmax = Cmax(1:length(V), :);
            S.Cmin = Cmin(1:length(V), :);
        else
            [S.Umax, S.Umin, S.Cmax, S.Cmin] = principal_curvature(S.VERT, S.TRIV);
        end
    end
    S.Cmean = 0.5 * (S.Cmax + S.Cmin);
    S.Cgauss = S.Cmax .* S.Cmin;
    if useNormalization
        S.Cmax = S.Cmax ./ median(S.Cmin);
        S.Cmin = S.Cmin ./ median(S.Cmin);
    end
    end
    S.Normal = vertexNormal(TR)';
    S.Normal = S.Normal';
else
    S.Normal = normals2D(S.VERT((1:boundary2d)', :));
    if computeCurvature
        if useSmoothing
            PTS = S.VERT((1:boundary2d)', :);
            n = size(PTS,1);
            UPSMPLPTS = interp1(1:n,PTS,linspace(1, n, 4 * n)); %4x as many points
            MAPPTS2UPSMPL = knnsearch(UPSMPLPTS, PTS);
            %smooth
            bnd = (1:size(UPSMPLPTS, 1))';
            SMOOTHUPSMPLPTS = curve_smooth(UPSMPLPTS, [bnd, [bnd(2:end); bnd(1)]], 'MaxIters', 10);
    
            tempS.VERT = SMOOTHUPSMPLPTS;
            curv2d = curv2dWrapper(tempS, length(SMOOTHUPSMPLPTS));
            S.BoundCurv = curv2d(MAPPTS2UPSMPL);
        else
            S.BoundCurv = curv2dWrapper(S,boundary2d);
        end
        if useNormalization
            S.BoundCurv = S.BoundCurv ./ median(S.BoundCurv);
        end
    end
end

%% reduce shape
if reductionFactor == 1
    if is2d
        S.boundary = (1:boundary2d)';
    end
    NNintoOriginal = (1:length(V));
    return;
end
NNintoOriginal = []; Vreduced = []; Freduced = [];
if is2d
    [Freduced, Vreduced] = reducepatch(S.TRIV, S.VERT, reductionFactor);
else
    if useMRed
        [Freduced, Vreduced] = reducepatch(S.TRIV, S.VERT, reductionFactor);
    else
        [Vreduced, Freduced, ~, NNintoOriginal] = decimate_libigl(S.VERT, S.TRIV, reductionFactor, 'Method', 'qslim');
        if isempty(Vreduced) % when trimesh contains holes :)
            [Freduced, Vreduced] = reducepatch(S.TRIV, S.VERT, reductionFactor);
        end
    end
end

if is2d
    % reorder reduced vertices such that boundary vertices come first
    boundaryEdges = freeBoundary(triangulation(Freduced, Vreduced));
    permutation = boundaryEdges(:, 1);
    if length(permutation) ~= length(Vreduced)
        idxRestPermutation = length(permutation) + 1;
        permutation = [permutation; zeros(length(Vreduced) - length(permutation), 1)];
        for v = 1:length(Vreduced)
            if all((permutation == v) == 0)
                permutation(idxRestPermutation) = v;
                idxRestPermutation = idxRestPermutation + 1;
            end
        end
    end
    try
    invPermutation(permutation) = 1:length(permutation)';
    catch
        disp('err');
    end

    S.TRIV = invPermutation(Freduced);
    S.VERT = Vreduced(permutation, :);
    S.boundary = (1:length(boundaryEdges))';
else
    S.TRIV = Freduced;
    S.VERT = Vreduced;
end

S.n = size(S.VERT, 1);
S.m = size(S.TRIV, 1);

% map original descriptors to reduced shape
if isempty(NNintoOriginal)
    NNintoOriginal = knnsearch(V, S.VERT);
    if is2d
        NNintoBoundOriginal = knnsearch(Vbound, S.VERT(S.boundary, :));
        NNintoOriginal(S.boundary) = NNintoBoundOriginal;
    end
end
S.wks = S.wks(NNintoOriginal, :);
S.hks = S.hks(NNintoOriginal, :);
S.gn_wks = S.gn_wks(NNintoOriginal, :);
S.rn_wks = S.rn_wks(NNintoOriginal, :);
S.shks = S.shks(NNintoOriginal, :);
S.gn_shks = S.gn_shks(NNintoOriginal, :);
S.rn_shks = S.rn_shks(NNintoOriginal, :);

if ~is2d
    if computeCurvature
        % map original curvatures to reduced shape
        S.Umin = S.Umin(NNintoOriginal, :);
        S.Umax = S.Umax(NNintoOriginal, :);
        S.Cmin = S.Cmin(NNintoOriginal);
        S.Cmax = S.Cmax(NNintoOriginal);
        S.Cmean = S.Cmean(NNintoOriginal);
        S.Cgauss = S.Cgauss(NNintoOriginal);
    end
    S.Dist2Other = S.Dist2Other(NNintoOriginal, :);
    S.Dist2OtherAll = S.Dist2OtherAll(NNintoOriginal, :);
    S.Normal = S.Normal(NNintoOriginal, :);
else
    S.Dist2Other = S.Dist2Other(NNintoOriginal(S.boundary), :);
    S.Normal = S.Normal(NNintoOriginal(S.boundary), :);
    if computeCurvature
        S.BoundCurv = S.BoundCurv(NNintoOriginal(S.boundary));
    end
end
end
