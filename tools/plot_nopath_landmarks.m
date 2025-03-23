function plot_nopath(M, X)
figure;
set(gcf, 'Position', [680, 258, 560, 420])


clickCounter = 1;
new_2d_landmarks = [];
new_3d_landmarks = [];


%% contour
ax1 = subplot(1,2,1);
cmp = jet(size(X,1));
temp1 = -min(X(:, 1)) + X(:, 1); temp1 = temp1./max(temp1); temp1 = 0.3 + 0.7 * temp1;
temp2 = -min(X(:, 2)) + X(:, 2); temp2 = temp2./max(temp2); temp2 = 0.2 + 0.6 * temp2;
cmp = [temp1 temp2 (linspace(0.7, 0.3, size(X, 1)))'];
Cplot = [X; X(1,:)];
patch('xdata',[Cplot(:,1); nan], 'ydata',[Cplot(:,2); nan],...
      'facevertexcdata', [cmp; cmp(1:2, :)], 'edgecolor', 'interp', 'linew', 4);
axis equal, axis off
hold(ax1, 'on')

%% mesh
ax2 = subplot(1,2,2);

whites = [1 1 1];
colormap(whites)
render_surf(M);
axis equal, axis off
hold(ax2, 'on')


set(gcf, 'WindowButtonDownFcn', @clickCallback);

function clickCallback(~, ~)
    ax = gca;
    
    % Odd click-> contour, even click-> mesh
    if mod(clickCounter, 2) == 1
        currentPoint = get(ax, 'CurrentPoint');
        clickedPoint = currentPoint(1, 1:2);
        [~, idx] = min(sum((X - clickedPoint).^2, 2));
        new_2d_landmarks = [new_2d_landmarks, idx];
        assignin('base', 'new_2d_landmarks', new_2d_landmarks);
        %disp(['Click ', num2str(clickCounter), ' on contour: Vertex Index ', num2str(idx)]);

        plot(ax1, X(idx,1), X(idx,2), 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
    else
        currentPoint = get(ax, 'CurrentPoint');
        clickedPoint = currentPoint(1, 1:3);  % Extract X, Y, Z coordinates
        [~, idx] = min(sum((M.VERT - clickedPoint).^2, 2));
        new_3d_landmarks = [new_3d_landmarks, idx];  % Store the index
        assignin('base', 'new_3d_landmarks', new_3d_landmarks);
        %disp(['Click ', num2str(clickCounter), ' on 3D mesh: Vertex Index ', num2str(idx)]);

        hold(ax2, 'on');
        sphereRadius = 0.03;
        [x, y, z] = sphere;
        x = sphereRadius * x + M.VERT(idx,1);
        y = sphereRadius * y + M.VERT(idx,2);
        z = sphereRadius * z + M.VERT(idx,3);
        surf(ax2, x, y, z, 'FaceColor', 'k', 'EdgeColor', 'none');
    end

    clickCounter = clickCounter + 1;
end

hold off

disp(new_3d_landmarks)
assignin('base', 'new_2d_landmarks', new_2d_landmarks);
assignin('base', 'new_3d_landmarks', new_3d_landmarks);

end