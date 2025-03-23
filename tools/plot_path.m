function plot_path(M, X, path, corr)
figure, 

%% contour
subplot(1,2,1), 
cmp = jet(size(X,1));
temp1 = -min(X(:, 1)) + X(:, 1); temp1 = temp1./max(temp1); temp1 = 0.3 + 0.7 * temp1;
temp2 = -min(X(:, 2)) + X(:, 2); temp2 = temp2./max(temp2); temp2 = 0.2 + 0.6 * temp2;
cmp = [temp1 temp2 (linspace(0.7, 0.3, size(X, 1)))'];
Cplot = [X; X(1,:)];
patch('xdata',[Cplot(:,1); nan], 'ydata',[Cplot(:,2); nan],...
      'facevertexcdata', [cmp; cmp(1:2, :)], 'edgecolor', 'interp', 'linew', 4);
axis equal, axis off

%% mesh
subplot(1,2,2),
if length(unique(corr)) ~= length(Cplot)
    path = [path; path(1)];
    corr = [corr; corr(1)];
    cmp = cmp(corr, :);
end

whites = [1 1 1];
colormap(whites)
render_surf(M);
% add NaN to prevent MATLAB to draw faces   
patch('xdata',[M.VERT(path,1); nan], 'ydata',[M.VERT(path,2); nan], 'zdata',[M.VERT(path,3); nan],...
      'facevertexcdata', [cmp; cmp(1, :)], 'edgecolor', 'interp', 'linew', 4);
axis equal, axis off

hold off
end