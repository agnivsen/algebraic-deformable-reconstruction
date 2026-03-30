classdef TriangleIndexer<handle
    %TRIANGLEINDEXER organises triangle indices for a point and its two
    %neighbors
    % 
    
    properties
        triangleMatrix;
        numberOfTriangles;
        numberOfEdges;
        edgeIndexer;
        numUsedEdges;
        uniqueEdgeIndices;
        allUsedEdgeList;
        numMaxNeighbors = 4;
    end
    
    methods
        function obj = TriangleIndexer(nng)
            % NNG adjustment
            kNTemp = size(nng, 2);
            if(kNTemp > obj.numMaxNeighbors)
                nng = nng(:,1:obj.numMaxNeighbors);
            end

            nPts = size(nng,1);
            kN = size(nng, 2);
            obj.edgeIndexer = EdgeIndexer(nng);
            
            triangleMat = zeros(nPts, nPts, nPts, 4);
            triangleIndex = 1;
            allUsedEdgeList = [];
            
            for ii = 1:nPts
                for jj = 1:(kN - 1)
                    p_j = ii; p_q = nng(ii,jj); p_r = nng(ii,jj+1);
                    pos = sort([p_j p_q p_r]);
                    % Checking if: (a) triangle has not already been
                    % indexed, AND, (b) q-th and r-th vertex are neighbors
                    % themselves
                    if( (triangleMat(pos(1), pos(2), pos(3)) == 0) && ( sum(nng(p_q,:) == p_r) ) )
                        % Edges follow the convention: [j-q, j-r, r-q]
                        a = obj.edgeIndexer.getIndex(pos(1), pos(2));
                        b = obj.edgeIndexer.getIndex(pos(1), pos(3));
                        c = obj.edgeIndexer.getIndex(pos(2), pos(3));
                        triangleMat(pos(1), pos(2), pos(3),:) = [triangleIndex a b c];
                        triangleIndex = triangleIndex + 1;
                        allUsedEdgeList = [allUsedEdgeList a b c];
                    end
                end
            end
            
            [uniqueEdges, ~, uniqueEdgeIndices] = unique(allUsedEdgeList);
            obj.triangleMatrix = triangleMat;
            obj.numberOfTriangles = triangleIndex - 1;
            obj.numberOfEdges = obj.edgeIndexer.numberOfEdges;
            obj.numUsedEdges = numel(uniqueEdges);
            obj.uniqueEdgeIndices = uniqueEdgeIndices;
            obj.allUsedEdgeList = allUsedEdgeList;
        end
        
        function [triangleIndex, edgeIndices] = getIndex(obj, p_j, p_q, p_r)
            pos = sort([p_j p_q p_r]);
            triangleIndex = squeeze(obj.triangleMatrix(pos(1), pos(2),pos(3), 1));
            edgeIndices = squeeze(obj.triangleMatrix(pos(1), pos(2),pos(3), 2:4));
        end

        function [utilisedEdgeID, matchCount] = getUtilisedEdgeIndex(obj, edgeID)
            matchCount = 0;
            for ii = 1:numel(obj.allUsedEdgeList)
                if(obj.allUsedEdgeList(ii) == edgeID)
                    matchCount = matchCount + 1;
                    utilisedEdgeID = obj.uniqueEdgeIndices(ii);
                end
            end
            assert(matchCount > 0, 'Incorrect edge ID supplied, no matches found. This has to be a severe bug!');
        end

    end
end

