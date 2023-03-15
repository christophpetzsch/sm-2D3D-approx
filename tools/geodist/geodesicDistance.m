function [D, adjListX] = geodesicDistance(X, rowIdx, colIdx)
    if ( isfield(X, 'TRIV') && ~isempty(X.TRIV) )
        adjListX = buildAdjacencyListFromFv(X);
        adjX = buildGraphMatrixFromAdjacencyList(adjListX, X.VERT);
    else
        nn = 3;
        
        adjListX = {};
        nnIdx = knnsearch(X.VERT, X.VERT, 'k', nn+1);
        for i=1:size(X.VERT,1)
            adjListX{i} = setdiff(nnIdx(i,:),i);
        end
        adjX = vertexAdjacency([], X.VERT);
    end
    adjX = 0.5*(adjX+adjX');
    
    [rows,cols,vals] = find(adjX);
    G = graph(rows,cols,vals);
    if ( ~exist('rowIdx', 'var') && ~exist('colId', 'var') )
        % compute full geodesic distances
        D = distances(G);
    else
        % faster for large graphs and small rowIdx/colIdx
        D = distances(G, rowIdx, colIdx);
%         for i=1:numel(rowIdx)
%             for j=1:numel(colIdx)
%                 [~, D(i,j)] = shortestpath(G, rowIdx(i), colIdx(j));
%             end
%         end
    end
end