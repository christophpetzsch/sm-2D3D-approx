function vNormals = normals2D( X )
    % X set of points in Nx3
    
    vecs = [X(end, :); X(1:end-1, :)] - X;
    
    % 90 degree rotation
    rot = [0 1 0; -1 0 0; 0 0 1];

    normals = vecs * rot';
    
    normals = normals ./ sqrt(sum(normals.^2, 2));
  
    
    [~, maxX] = max(X(:, 1));
    %[~, maxY] = max(X(:, 2));
    [~, minX] = min(X(:, 1));
    %[~, minY] = min(X(:, 2));
    XPlusNormals = X + 0.1 * (X(maxX, 1) - X(minX,1)) * normals;
    
    % bounding box must get bigger for outside pointing normals
    if norm(X(maxX,:) - X(minX,:)) > norm(XPlusNormals(maxX,:) - XPlusNormals(minX,:))
        normals = -normals;
    end
    
    
    % transform edge normals to vertex normals by a weighted sum of
    % adjacent edge normals
    vecLengths = sqrt(sum(vecs.^2, 2));
    
    normals = [normals; normals(1, :)];
    vecLengths = [vecLengths; vecLengths(1)];
    
    vNormals = (vecLengths(1:end-1) ./ (vecLengths(1:end-1) + vecLengths(2:end))) .* normals(1:end-1, :)...
             + (vecLengths(2:end)   ./ (vecLengths(1:end-1) + vecLengths(2:end))) .* normals(2:end, :);
    
    vNormals = vNormals ./ sqrt(sum(vNormals.^2, 2)); 
    
%     delta = norm(max(X) - min(X));
%     XPlusVNormals = X + delta * 0.01* vNormals;
%     plot(X(:, 1), X(:, 2))
%     hold on;
%     for i = 1:length(X)
%         plot([X(i, 1) XPlusVNormals(i, 1)], [X(i, 2) XPlusVNormals(i, 2)], 'r')
%     end
   
end