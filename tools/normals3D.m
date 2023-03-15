function [ normals ] = normals3D( M )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here


e21 = M.VERT(M.TRIV(:,2),:) - M.VERT(M.TRIV(:,1),:);
e31 = M.VERT(M.TRIV(:,3),:) - M.VERT(M.TRIV(:,1),:);
tri_n = cross(e21,e31);

normals = zeros(M.n,3);

for i=1:M.n
    ind = (M.TRIV(:,1) == i | M.TRIV(:,2) == i | M.TRIV(:,3) == i);
    normals(i,:) = sum(tri_n(ind,:),1) / sum(ind);
end

normals = normals ./ repmat( sqrt(sum(normals.^2, 2)), 1, 3 );

end
