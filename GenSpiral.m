function ...
   IMAGE = GenSpiral(M,PCONST,PSLOPE,RADIUS,THICK,INVSNR,GRADIENT,...
   FILESAVE,BULGERADIUS)
 
% This function generates an artificial spiral.  It outputs the square 
% bitmap matrix IMAGE of the spiral on empty space, where the spiral has 
% pixel value MAXPIX (for now, MAXPIX = 255) and empty space has pixel 
% value 0.  The function also creates and saves a .fits image and a .jpg 
% image of the spiral.


% ***** INPUTS
%
% M is the number of spiral arms.
% 
% The user inputs the pitch angle as a function of radius.  For now,
% the function must be linear: P(RR) = PSLOPE*RR + PCONST, where PSLOPE is 
% the slope in degrees per outer radius, RR is the radial coordinate in 
% in units of outer radii, and and PCONST is the intercept in degrees.  
% Note that RR is not an input.  Put simply, the pitch angle varies 
% linearly from PCONST at the center to PCONST + PSLOPE at the outer 
% radius.
%
% The output spiral has outer radius RADIUS in pixels and spiral arm 
% thickness THICK in pixels.  
%
% INVSNR determines how noisy the image will be.  Specifically it is the
% inverse of the signal-to-noise ratio.  For a noiseless image, INVSNR = 0.
% For an image with SNR = 10, INVSNR = 0.1.
%
% GRADIENT determines whether the spiral will have a light profile like a
% real galaxy or whether the spiral will have the same intensity
% throughout:
%
%   If GRADIENT = 0, then all pixel on the spiral have the same
%   pixel value, subject to the noise.  Empty space has the same amount of 
%   noise.  
%
%   If GRADIENT = 1 or 2, then the spiral's luminosity is modeled after UGC 
%   463, as computed via bulge-disk decomposition by Matt Hartley:
%
%   I = I_bulge + I_disk, where R is the radial coordinate, assuming that 
%   the visible edge of the galaxy has radius 70 pixels, and:
%
%       I_bulge = 2300*exp(-(R/2.7))^(1/1.3) 
%       I_disk = 454*exp(-R/127)
%
%   GRADIENT = 1 is more like a nearby galaxy, while GRADIENT = 2 is more 
%   like a high-redshift galaxy. 
%
%   For GRADIENT = 1, the SNR is constant throughout the galaxy, as though 
%   the galaxy itself were supplying the noise. 
%
%   For GRADIENT = 2, the standard deviation of the pixels (not the SNR) 
%   remains constant throughout the image, as though the detector were
%   supplying the noise.
%
% FILESAVE determines whether to save the IMAGE matrix to image files.  If 
% FILESAVE is set to 0, then the function only outputs a bitmap matrix.   
% Otherwise the bitmap matrix will be saved as both a FITS and a JPEG.


% Limitations on inputs include:
%     M must be a positive integer.
%     P should be between -90 and 90, and must not be 0.
%     RADIUS must be an integer.
%     THICK must be an integer.
%     GRADIENT must be 0, 1 or 2.


% ***** OUTSIDE FUNCTIONS USED
%
% PeriodToDash.m
% fitswrite.m


% TROUBLESHOOTING
%
% This code relies on fitswrite.m to write the .fits files.  If fitswrite.m
% writes the wrong pixel values, you might need to tell fitswrite.m to 
% bitswap.  
%
% Fitswrite.m contains a series of if/else statements that looks for your
% architecture to determine whether bitswapping is necessary. The code
% isn't familiar with all architectures.  For my Imac64, I had to insert:
% 
%       elseif strmatch(friend,'MACI64')
%         bswap = 'b';
%
% into the series of if/else statements.

% ***** METHOD
%
% The code computes a mathematical spiral with its origin at the center of
% the image.  It then establshes points spaced D pixels apart along each 
% spiral arm.  The code adds thickness to the spiral by recognizing 
% points perpendicular to the spiral, also D pixels apart, along the 
% spiral arm.  If one or more of these recognized points lies within a
% given pixel, then the pixel's value is 1.  Otherwise the pixel value
% remains zero.
%
% For now we let the point spacing be 0.25 pixels.



% If the bulge radius is bigger than the spiral radius, throw a warning.

if BULGERADIUS >= RADIUS
    warning(['Bulge radius equals or exceeds spiral radius. '...
        'Output will be a solid disk.'])
    BULGERADIUS = RADIUS;
end


% If the bulge is bigger than the bar, throw a warning.
%
% if BULGERADIUS >= BARHALFLENGTH
%    if BARHALFLENGTH > 0
%        warning(['Bulge radius equals or exceeds bar half length. '...
%        'Output will not include a bar.'])
%    end
% end
    
% ***** COMPUTATIONS

D = 0.25;                 % Point spacing in pixels
% P = P*pi/180;      
PCONST = PCONST*pi/180;   % Convert pitch angles to radians
PSLOPE = PSLOPE*pi/180;

% B = -tan(P);       % Compute the spiral constant
MAXPIX = 1;          % Value of the brightest pixels, excluding noise

% Set up the image matrix

N = 2*RADIUS + 2*THICK + 1;   % Length of the side of the image, in pixels
IMAGE = zeros(N,N);           % Initialize matrix
X0 = RADIUS + THICK + 1;      % Establish the origin
Y0 = RADIUS + THICK + 1;


% Compute the arclength S, in pixels, of each spiral arm

dR = 0.25;                          % Radial derivative element
S = 0;                              % Initialize arclength
for RR = dR:dR:RADIUS               % For each radial element
    PP = PCONST + PSLOPE*RR/RADIUS;    % Compute local pitch
    B = -tan(PP);                      % Comput local spiral constant
    dSdR = sqrt(1 + B^2)/abs(B);       % Compute arclength derivative dS/dR
    S = S + dR*dSdR;                   % New arclength
end


% S = abs(RADIUS/B)*sqrt(1+B*B);


NT = S/D;   % Number of points on each spiral arm, at D pixels per point.

% Compute the spiral.  Temporary variables include:
%
%   SS = the arclength at a given point
%
%   RR = the radius at a given point on the spiral proper
%   DR = change in radius (radial derivative element) of the spiral
%   RRTHICK = radius at a given point on the thickness line
%
%   THTH = the azimuth at a given point on the spiral proper
%   RDTH = radius times the change in azimuth (azimuthal derivative 
%            element) of the spiral
%   THTHTHICK = azimuth at a given point on the thickness line
%
%   NORM = normalization coefficient of the derivative elements
%
%   XX = the x-value of the pixel encompassing a given point
%   YY = the y-value of the pixel encompassing a given point



for TH0 = 2*pi/M:2*pi/M:2*pi      % For each spiral arm
%   SREAL = 0;                        % Arclentgh test quantity
    SS = D;                           % Initialize arclength  
    RR = D;                           % Initialize radial coordinate
    THTH = -TH0;                      % Initialize azumithal coordinate
    XX = X0;                          % Initialize Cartesian coordinates
    YY = Y0;
    for TSTEP = 1:NT              % For each point on the spiral arm
        
%        SS = TSTEP*D;
%        B = abs(tan(PCONST + PSLOPE*RR/RADIUS)); % Abs val of spiral const
%        DR = SS*abs(B)/sqrt(1+B*B) - RR;
%        RR = RR + DR;
%        N;
%        RDTH = RR*(log(RR)/B - TH0 - THTH);
%        THTH = THTH + RDTH/RR;
%        NORM = sqrt(RDTH*RDTH+DR*DR);
      
        PP = PCONST + PSLOPE*RR/RADIUS;    % Local pitch
%        PPdeg = PP*180/pi                
 
        B = tan(PP);                       % Local spiral constant
        dS = D;                            % Arclength element
        dSdR = sqrt(1 + B^2)/abs(B);
        dR = dS / dSdR;                    % Radial element
        THTH;
%        dTHdR = 1 / (B * exp(B*THTH))
        dTHdR = 1 / (B * RR);
        dTH = dR * dTHdR;                  % Azimuthal element
        
        SS = SS + dS;                % New arclength coordinate
        RR = RR + dR;                % New radial coordinate 
        THTH = THTH + dTH;           % New azimuthal coordinate
        
        RdTH = RR * dTH;             % Azimuthal distance element
        NORM = sqrt(RdTH^2 + dR^2);  % Normalization constant
        
        for THICKSTEP = (-THICK/2+0.5):D:(THICK/2-0.5)   
                                                 % For each point on the 
                                                 % thickness line, perpen-
                                                 % dicular to the spiral
            RRTHICK = RR + THICKSTEP*RdTH/NORM;         % Establish points
            THTHTHICK = THTH - THICKSTEP*(dR/RR)/NORM;     
            XX = round(X0 + RRTHICK*cos(THTHTHICK));    % Find the pixel
            YY = round(Y0 + RRTHICK*sin(THTHTHICK));   
            IMAGE(XX,YY)=MAXPIX;                        % Color the pixel
%            IMAGE(XX,YY)=IMAGE(XX,YY)*...               % Apply circular
%                cos(pi*THICKSTEP/(THICK-1));            %    cross section
        end
    end
end

% Apply bar

% For all pixels
% Compute radial coordinate
% If pixel is inside bar half length
% Erase pixel
% How the fuck do i do that?



% Apply bulge

for XX = 1:N
    for YY = 1:N                               % For all pixels
        RR = sqrt((XX - X0)^2 + (YY - Y0)^2);  % Compute radial coordinate
        if RR <= BULGERADIUS                   % If pixel is inside bulge
            IMAGE(XX,YY) = MAXPIX;             % Draw pixel
        end
    end
end


% Apply noise and gradient

if GRADIENT==0                                        % If no gradient
    IMAGE = IMAGE + MAXPIX*INVSNR*randn(N);              % just add noise
else
    for XX = 1:N
        for YY = 1:N
            RR = sqrt((XX - X0)^2 + (YY - Y0)^2);     % Radial coordinate
            RR = RR*70/RADIUS;                        % Scale RR to 70-
                                                      %   pixel radius
            I_bulge = 2300*exp(-(RR/2.7)^(1/1.3));    % Bulge profile
            I_disk = 454*exp(-RR/127);                % Disk profile
            I = (I_disk + I_bulge)/(2300+454);        % Scale the gradient            
            IMAGE(XX,YY) = IMAGE(XX,YY)*I;            % Apply the gradient

            if GRADIENT==1                             % If low-redshift
                IMAGE(XX,YY) = ...                     %   then apply noise
                    IMAGE(XX,YY) + I*INVSNR*randn(1);  %   here
            else                                       % If high-redshift
                IMAGE(XX,YY) = IMAGE(XX,YY) + ...      %   then apply noise
                     MAXPIX*INVSNR*randn(1);           %   here
            end
            
                
        end
    end
end


if FILESAVE ~= 0    % If we want to generate image files

    % Establish the filename of the output .fits anf .jpg images.
    PCONST = PCONST*180/pi;          % Convert pitches back to degrees
    PSLOPE = PSLOPE*180/pi;          

    MSTRING = num2str(M);             % Convert inputs to strings

    PCONSTSTRING = num2str(PCONST);                       
    PCONSTSTRING = PeriodToDash(PCONSTSTRING);

    PSLOPESTRING = num2str(PSLOPE);                       
    PSLOPESTRING = PeriodToDash(PSLOPESTRING);

    RADIUSSTRING = num2str(RADIUS);

    THICKSTRING = num2str(THICK);

    INVSNRSTRING = num2str(INVSNR);
    INVSNRSTRING = PeriodToDash(INVSNRSTRING);

    if GRADIENT == 0
        GRADIENTSTRING = 'NoGrad';
    else
        if GRADIENT == 1
            GRADIENTSTRING = 'LowZ';
        else
            if GRADIENT == 2
                GRADIENTSTRING = 'HighZ';
            end
        end
    end

    BULGERADIUSSTRING = num2str(BULGERADIUS);
    BULGERADIUSSTRING = PeriodToDash(BULGERADIUSSTRING);
        
    FNAME = ['Spiral_' MSTRING '_' PCONSTSTRING '_' PSLOPESTRING '_' ...
        RADIUSSTRING '_' THICKSTRING '_' INVSNRSTRING '_' ...
        GRADIENTSTRING '_' BULGERADIUSSTRING];
        FITSNAME = [FNAME '.fits'];
    JPGNAME = [FNAME '.jpg'];

    % Write the .fits filename
    fitswrite(IMAGE',FITSNAME);   % Generate the .fits file.
                                  % The .jpg writer flips the array, so we 
    IMAGE2 = fliplr(IMAGE);       %    flip it back. 
    imwrite(IMAGE2,JPGNAME);      % Generate the .jpg.

end
return






    
