

% **** ALTER CODE BEGINNING HERE ****

% FSTRING = 'Synthetic28';


% Identify the computer and navigate to the appropriate folder.

% [~,COMPUTER] = system('hostname');
% if COMPUTER(1:6) == 'scooby'
%    CDSTRING = ['/Users/Doug/Documents/Pitch_m/doublepeak/' FSTRING];
% else if COMPUTER(1:3) == 'Taz' 
%    CDSTRING = ['/Users/Doug/Documents/AGES/Side_Projects/Pitch_m/'...
%    'DoublePeak/' FSTRING];
%    else if COMPUTER(1:7) == 'Foghorn'
%            CDSTRING = ['Ben Ben Ben Ben Ben Ben Ben Ben' FSTRING]
%        end
%    end
% end
% cd(CDSTRING) 

% cd /Users/Doug/Documents/AGES/Side_Projects/Pitch_m/TestTest


% Input pitch.m parameters
FSTRING = 'Synthetic28';
FILE = 'Spiral_5_28_0_50_2_0_NoGrad.fits';
X0 = 53;                % x-coordinate of the spiral center
Y0 = 53;                % y-coordinate of the spiral center
VIS_INNER = 6;          % Visible inner radius of spiral
VIS_OUTER = 32;         % Visible outer radius of spiral
MSMT_INNER1 = 6;        % Innermost inner radius of measurement annuli 
MSMT_INNER2 = 20;       % Outermost inner radius of measurement annuli
InnerRadiusSpacing = 2; % Spacing between consecutive inner radii
MSMT_OUTER = 45;        % Outer radius of measurement annuli
NAXIS = 100;            % Number of spiral axes to compute for each pitch angle
MINP = 20;              % Minimum of pitch angle domain
MAXP = 35;              % Maximum of pitch angle domain
PSTEP = 1;             % Spacing between points on pitch angle domain
AxisPointSpacing = 0.1; % Spacing, in pixels, between computation points on each axis 
SMOOTH=1;               % Smooth the Variance of Means (fitness) function?  
                        %   (1 = yes)
Save2D=0;               % Save a 2-D graph for each inner radius (1 = yes)?  
Save3D=1;               % Save a single 3-D graph for all inner radii
                        %   (1 = yes)?


% Symmetry to be computed

SYM = 4;   



% **** NO NEED TO ALTER CODE BEYOND THIS POINT ****


% Import image

IMAGE = fitsread(FILE,'full',0);
C0 = X0;
R0 = Y0;

% Generate filenames

SYMSTRING = num2str(SYM);
SYMFILE = [FSTRING 'SingleSym' SYMSTRING '.fits'];
RESIDFILE = [FSTRING 'SingleResid' SYMSTRING '.fits'];

% Compute the symmetrical component

[SYMMETRIC,RESIDUAL] = Sympart(IMAGE,SYM,C0,R0);
fitswrite(SYMMETRIC', SYMFILE);
fitswrite(RESIDUAL',RESIDFILE);

% Measure the pitch angle of the symmetric component

FILE = SYMFILE;
                 
                    
[PITCHvsINNER, BESTFITPITCH, ERROR] = Spirality(FILE,X0,Y0,...
    VIS_INNER,VIS_OUTER,MSMT_INNER1,MSMT_INNER2,InnerRadiusSpacing,...
    MSMT_OUTER,NAXIS,MINP,MAXP,PSTEP,AxisPointSpacing,SMOOTH,Save2D,Save3D)

close all

