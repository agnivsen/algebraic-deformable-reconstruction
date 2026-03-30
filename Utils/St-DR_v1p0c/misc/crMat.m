function M = crMat(v,n)
% M = crMat(v,n)
%
% Returns the (3x3) matrix M such that Mz = v x z for all z \in \R^3.
% The second optional argument makes the function return only the first two
% rows of the crMat, ie a (2x3) matrix.

if nargin == 1    
    M = [
        0 -v(3) v(2)
        v(3) 0 -v(1)
        -v(2) v(1) 0
        ];
else
    M = [
        0 -v(3) v(2)
        v(3) 0 -v(1)
        ];
end;