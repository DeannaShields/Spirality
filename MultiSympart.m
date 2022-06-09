function [MULTI_SYM,MULTI_RESID] = MultiSympart(IMAGE,HIGH_M,C0,R0)

% Doug Shields
% February 2014

% This function subtracts the asymmetric parts (or, equivalently, the high-
% mode symmetric parts) from the IMAGE matrix.  Specifically, the code 
% keeps the 2-fold, 3-fold, ... , and HIGH_M-fold symmetric parts in 
% the MULTI_SYM matrix and puts the remainder in the MULTI_RESID matrix.
%
% The center of symmetry is the point at row R0 and column C0.
%
% HIGH_M mist be an integer greater than 1.

% Other functions needed:
%
% - Sympart.m
% - RotateTheta.m
% fitsread.m

% QUIRKS
%
% When using FITSREAD to import a file, the center (X0, Y0) maps to (C0,
% R0).


SZ = size(IMAGE);
ROWS = SZ(1);
COLS = SZ(2);
MULTI_SYM = zeros(ROWS,COLS);
RESID = IMAGE;                               % Initialize outputs

for M = HIGH_M:-1:2                          % For every symmetry mode  
                                             %    you want to keep 
                                             
    IMAGE = RESID;                           % Look for symmetry in the
                                             %    previous mode's residual
                                             
    [SYM,RESID] = Sympart(IMAGE,M,C0,R0);
    MULTI_SYM = MULTI_SYM + SYM;             % Add the symmetry to the 
end                                          %    output MULTI_SYM

MULTI_RESID = RESID;

return
    
    
 