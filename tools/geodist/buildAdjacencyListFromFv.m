function adjacencyList = buildAdjacencyListFromFv(fv)

    if ( size(fv.TRIV,2) == 4 )
        fv = quadrilateral2triangluar(fv);
    elseif ( size(fv.TRIV,2) > 4 )
        error('Only triangular and quadrilateral meshes are supported');
    end
    vertexFeatures = fv.VERT;
    faces = fv.TRIV;
    
    nFaces = size(faces,1);
    N = size(vertexFeatures,1);

    adjacencyList = cell(1,N);
    
    for f=1:nFaces
        neighs = faces(f,:);

        for n=1:numel(neighs)
            % add all vertices but neighs(n)
            adjacencyList{neighs(n)} = [adjacencyList{neighs(n)} neighs(neighs~=neighs(n))];
        end
    end
end