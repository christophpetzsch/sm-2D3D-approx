function G = buildGraphMatrixFromAdjacencyList(adjacencyList, vertexFeatures)
    N = numel(adjacencyList);
    
    G = sparse(N,N);
    for l=1:N
        neighs = setdiff(unique(adjacencyList{l}),l);
        
        for neighIdx=1:numel(neighs)
            currNeigh = neighs(neighIdx);
            
            % set lower triangle of G
            if ( currNeigh > l )
                rowIdx = currNeigh;
                colIdx = l;
            else
                rowIdx = l;
                colIdx = currNeigh;
            end
            
            if ( G(rowIdx,colIdx) == 0 )
                if ( exist('vertexFeatures', 'var')  )
                    G(rowIdx,colIdx) = norm(vertexFeatures(rowIdx,:)-vertexFeatures(colIdx,:));
                else
                    G(rowIdx,colIdx) = 1;
                end
            end
        end
    end
end