classdef DeltaTensorIndexer<handle
    %EdgeIndexer Organises edge indices for NNG
    % 
    
    properties
        edgeMatrix;
        numberOfEdges;
    end
    
    methods
        function obj = DeltaTensorIndexer(nng)
            nPts = size(nng,1);
            kN = size(nng, 2);
            
            edgeMat = zeros(nPts, nPts);
            edgeIndex = 1;
            
            for ii = 1:nPts
                for jj = 1:(kN-1)
                    p_j = ii; p_q = nng(ii,jj);
                    pos = sort([p_j p_q]);
                    if(edgeMat(pos(1), pos(2)) == 0)
                        edgeMat(pos(1), pos(2)) = edgeIndex;
                        edgeIndex = edgeIndex + 1;
                    end
                    
                    p_j = nng(ii,jj); p_q = nng(ii,jj+1);
                    pos = sort([p_j p_q]);
                    if(edgeMat(pos(1), pos(2)) == 0)
                        edgeMat(pos(1), pos(2)) = edgeIndex;
                        edgeIndex = edgeIndex + 1;
                    end
                    
                    p_j = ii; p_q = nng(ii,jj+1);
                    pos = sort([p_j p_q]);
                    if(edgeMat(pos(1), pos(2)) == 0)
                        edgeMat(pos(1), pos(2)) = edgeIndex;
                        edgeIndex = edgeIndex + 1;
                    end
                    
                end
            end
            
            obj.edgeMatrix = edgeMat;
            obj.numberOfEdges = edgeIndex - 1;
        end
        
        function [index] = getIndex(obj, p_j, p_q)
            pos = sort([p_j p_q]);
            if(obj.edgeMatrix(pos(1), pos(2)) == 0)
                fprintf('<strong>Error</strong> while quering edge indices [%d, %d]\n', pos(1), pos(2));
                error('EdgeIndexer: wrong edge position queried!');
            else
                index = obj.edgeMatrix(pos(1), pos(2));
            end
        end
    end
end

