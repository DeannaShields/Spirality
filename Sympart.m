function [SYM,RESID] = Sympart(IMAGE,M,C0,R0)

% Doug Shields
% February 2014

% This function expresses an image (matrix) as the sum of its 
% M-fold symmetric part and its residual, where M = 360 degrees divided by
% the symmetry angle.
%
% The function compares every matrix element to its M-1 conjugate elements
% symmetric about the point at the R0th row and the C0th column.  The 
% function the reduces the higher values to the lowest value, resulting 
% in the M-fold symmetric matrix SYM. 
%
% Mathematically, SYM = (M-1)*IMAGE - XXXX
%   where XXXX is sum from j=1 to j=m-1 of (IMAGE - IMAGE rotated by 
%   2*pi*j/M)_truncated.
%
% In the above formula, "truncated" refers to negative values being changed
% to zero.  The method is described in more detail in Elmegreen et al, 1992 
% (ApJS, 79, 37).
%
% The difference between the original matrix IMAGE and the symmetric matrix
% SYM is the residual matrix RESID.


% OTHER FUNCTIONS NEEDED
%
% RotateTheta.m

% QUIRKS
%
% When using FITSREAD to import a file, the center (X0, Y0) maps to (C0,
% R0).

SZ = size(IMAGE);
ROWS = SZ(1);
COLS = SZ(2);
SYM = zeros(ROWS,COLS);    
RESID = zeros(ROWS,COLS);                 % Initialize outputs

SYM = (M-1)*IMAGE;

for j=1:1:M-1                             % For every image rotation
    THETA=2*pi*j/M;
    IMROT=RotateTheta(IMAGE,THETA,R0,C0); % Rotate the image
    SUB = IMAGE - IMROT;                  % Subtract rotated from original
    SUB = 0.5*(SUB+abs(SUB));             % Truncate
    SYM = SYM - SUB;
end

RESID = IMAGE - SYM;
return
    