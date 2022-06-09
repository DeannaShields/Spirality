function [ROTATED] = RotateTheta(INMATRIX,THETA,R0,C0)

% Doug Shields
% February 2014

% This function rotates the matrix INMATRIX by THETA radians 
% counterclockwise about the point at the R0th row and the C0th column.  
% The output matrix, ROTATED, has the same dimensions as the input, 
% INMATRIX.  The unused elements are left as zeros.

SZ = size(INMATRIX);
ROWS = SZ(1);
COLS = SZ(2);
ROTATED = zeros(ROWS,COLS);        

R = [cos(THETA) -sin(THETA);sin(THETA) cos(THETA)];  %Rotation matrix

for ROW = 1:ROWS                         % For every row & col in the 
    for COL = 1:COLS                     %    input matrix
        
        ROW_CTR = ROW-R0;        % Counting rows & cols from the
        COL_CTR = COL-C0;        %    pivot
        NEW_CTR = R*[ROW_CTR;COL_CTR];   % New location of the element
        NEWROW = NEW_CTR(1)+R0;  % Back to counting rows & cols
        NEWCOL = NEW_CTR(2)+C0;  %   from the upper left
%        NEWROW = round(NEWROW);
%        NEWCOL = round(NEWCOL);
        if NEWROW<ROWS+1
            if NEWROW>0                  % If any part of the new pixel
                if NEWCOL<COLS+1         %     lies on the new matrix grid
                    if NEWCOL>0
                        
                        ELEMENT = INMATRIX(ROW,COL);
                        NEWROWMOD = mod(NEWROW,1);  % Just decimal part 
                        NEWCOLMOD = mod(NEWCOL,1);  %   of new row & col
                        NEWROWFLR = floor(NEWROW);  % Just the whole 
                        NEWCOLFLR = floor(NEWCOL);  %   number part
                        
                        % Now we divvy the element value among all its
                        % adjacent pixels, depending on how much of the old
                        % pixel falls onto each new pixel.
                        
                        if NEWROW>=1
                            if NEWCOL>=1
                                ROTATED(NEWROWFLR,NEWCOLFLR) = ...
                                    ROTATED(NEWROWFLR,NEWCOLFLR) + ...
                                    ELEMENT * (1 - NEWROWMOD) * ...
                                    (1 - NEWCOLMOD);
                            end
                        end
                        
                        if NEWROW>=1
                            if NEWCOL<COLS
                                ROTATED(NEWROWFLR,NEWCOLFLR+1) = ...
                                    ROTATED(NEWROWFLR,NEWCOLFLR+1) + ...
                                    ELEMENT * (1 - NEWROWMOD) * ...
                                    NEWCOLMOD;
                            end
                        end
                        
                        if NEWROW<ROWS
                            if NEWCOL>=1
                                ROTATED(NEWROWFLR+1,NEWCOLFLR) = ...
                                    ROTATED(NEWROWFLR+1,NEWCOLFLR) + ...
                                    ELEMENT * (NEWROWMOD) * ...
                                    (1 - NEWCOLMOD);
                            end
                        end
                        
                        if NEWROW<ROWS
                            if NEWCOL<COLS
                                ROTATED(NEWROWFLR+1,NEWCOLFLR+1) = ...
                                    ROTATED(NEWROWFLR+1,NEWCOLFLR+1) + ...
                                    ELEMENT * NEWROWMOD * ...
                                    NEWCOLMOD;
                            end
                        end
                                               
%                        ROTATED(NEWROW,NEWCOL) = INMATRIX(ROW,COL);
                    end
                end
            end
        end
    end
end

return
        


