function ...
    [PITCHvsINNER, BESTFITPITCH, ERROR] = Spirality(FILE,X0,Y0,...
    VIS_INNER,VIS_OUTER,MSMT_INNER1,MSMT_INNER2,InnerRadiusSpacing,...
    MSMT_OUTER,NAXIS,MINP,MAXP,PSTEP,AxisPointSpacing,SMOOTH,Save2D,Save3D)




% INTRODUCTION
%
% This function measures a spiral galaxy's best-fit pitch angle.  The code 
% computes spiral coordinate systems, or templates, for every pitch angle 
% from MINP to MAXP degrees in steps of PSTEP.  It then attempts to fit the 
% galaxy to each template, and outputs a fitting function called Variance 
% of Means.  The galaxy's best-fit pitch corresponds to a peak in the 
% fitting function.
%
% COMPUTING THE FITTING FUNCTION: For every pitch angle template, the code 
% computes the mean pixel value along each spiral axis.  The variance of 
% these means is recorded as a point in the fitting function, Variance of 
% Means vs. Pitch.  For a noiseless, logarithmic spiral, the fitting 
% function shows an absolute maximum at the spiral's true pitch angle.  For 
% a noisy galaxy image, the fitting function shows a local maximum at the 
% galaxy's best-fit pitch angle.
% 
% ERROR BARS: In order to find the one-sigma confidence interval, the inner
% radius of the measurement annulus is varied from MSMT_INNER1 to 
% MSMT_INNER2 equidistant steps of InnerRadiusSpacing pixels.  The mean of 
% these measurements is considered the best-fit pitch BESTFITPITCH.  
%
% The random error is the standard deviation of these measurements.  The 
% random error is scaled in order to penalize small stable measurement 
% regions and reward large ones.  The scale factor is the difference 
% between the visible inner and outer radii of the galaxy, VIS_OUTER - 
% VIS_INNER, divided by the difference between the two inner radii, 
% MSMT_INNER2 - MSMT_INNER1.  The scaled result, or the "stretched error" 
% is added in quadrature with the pitch angle domain spacing PSTEP, 
% resulting in the total error ERR.



% OUTPUTS: PASSED VARIABLES
%
% PITCHvsINNER - A two-column array showing pitch angle in degrees as a 
% function of inner measurement radius.
%
% BESTFITPITCH - The mean pitch angle in the PITCHvsINNER array, or the 
% best-fit pitch angle of the galaxy.
%
% ERR - The total error in the BESTFITPITCH measurement.  It is the
% standard deviation of the pitch angles in the PITCHvsINNER array, scaled 
% by the radial range of the visible spirals divided by the range of inner 
% measurement radii, then added in quadrature with the spacing between
% successive pitch angle templates.



% INPUTS
%
% FILE - The filename of the galaxy image.  The file must be in *.fits 
% format, and the galaxy must be face-on or deprojected to circular.
%
% X0, Y0 - The center of the galaxy in Cartesian pixel coordinates.  Be 
% sure to use "image" coordinates, not "physical" coordinates.  Must be 
% positive reals.
%
% VIS_INNER, VIS_OUTER - Visually estimated inner and outer radii, in 
% pixels, of the galaxy's spirals.  These inputs are used to compute the 
% error bar, not to compute pitch angle itself.  Must be positive reals 
% such that VIS_INNER < VIS_OUTER.
%
% MSMT_INNER1, MSMT_INNER2, MSMT_OUTER - The code first measures the galaxy
% on an annulus with inner radius MSMT_INNER1 and outer radius MSMT_OUTER.  
% It then repeats the process, increasing the inner radius incrementally.  
% The final measurement is on an annulus from MSMT_INNER2 to MSMT_OUTER.  
% The best-fit pitch is the mean of pitch angles measured at all inner 
% radii from MSMT_INNER1 to MSMT_INNER2.   Must be positive reals, such 
% that MSMT_INNER1 < MSMT_INNER2.  All radii are in pixels.
%
% InnerRadiusSpacing - Spacing, in pixels, between inner radii.  If the
% user wishes to measure only one inner radius, there are two ways to
% accomplish this:
%
%    A) Set InnerRadiusSpacing=0.  In that case, only the inner radius
%    MSMT_INNER1 will be measured, while MSMT_INNER2 will be ignoresd.  
% 
%    or
%
%    B) Set MSMT_INNER2 = MSMT_INNER1.  In that case, only the inner radius
%    defined by MSMT_INNER2 and MSMT_INNER1 will be measured.  The value of
%    InnerRadiusSpacing will be ignored.
%
% MSMT_OUTER - The outer radius, in pixels, of all measurement annuli.  
% Ideally, it should be set just outside the visible outer radius of the 
% spirals, in empty space.  However, it should not include pixel 
% contamination from foreground stars.  Must be a real number greater than 
% MSMT_INNER2.
%
% NAXIS - Number of spiral axes in each spiral template.  1000 is often
% sufficient.  Insufficient values of NAXIS will result in high-frequency,
% periodic variations in the fitting function, partcularly in the loose
% end of pitch angle domain (that is, as |P| approaches 90.  In order to 
% ensure that each pixel is measured at least once, NAXIS should be at 
% least 2*pi*MSMT_OUTER.  Must be a positive integer.
%
% MINP, MAXP - Minimum and maximum pitch angles, in degrees, of the spiral 
% templates created for fitting.  Must be real numbers such that
% -90 <= MINP <= MAXP <= 90.
%
% PSTEP - Spacing, in degrees, between pitch angles of successive
% templates.  
%
% AxisPointSpacing - The spacing, in pixels, between computation points on
% each axis.  Must be a positive integer.  The code works very well with a 
% value of 0.1.  Lower values will asymptote the pitch angle result toward 
% the correct value.  Computation time varies inversely with this quantity. 
%
% SMOOTH - A toggle for applying a 5-point moving average to the Variance
% of Means vs. pitch angle fitting function.  If SMOOTH is 1, the moving
% average is applied; otherwise it is not.  This feature be useful in 
% smoothing high-frequency variations caused by an insufficient value of 
% NAXIS.  However, also affect the location of the peaks, so use with 
% caution.
%
% Save2D and Save3D - Toggles for saving the output files.  If Save2D==1,
% a 2-D graph of the fitting function vs. pitch angle will be generated for
% each inner radius.  If Save3D==1, a 3-D graph of the fitting function vs.
% pitch angle and inner radius will be generated.  If either variable is
% set to 1, then a text file summarizing the results will be generated.
%



% OUTPUTS: FILES GENERATED
%
% If Save2D==1 or if Save3D=1, the following file is generated:
%
%   - *.txt: A summary text file showing the passed 
%     variables PITCHvsINNER, BESTFITPITCH, and ERR.
%
% If Save2D==1, the following files is generated:
%
%   - 2D*.fig: A MATLAB figure showing the fitting function, Variance
%     of Means vs. pitch angle.  One such file is created for each inner
%     radius.
%
%   - 2D*.eps: Same as *.fig, but an .eps image.
%
% If Save3D==1, the following files is generated:
%
%   - 3D*.fig: A MATLAB figure showing the fitting function, Variance
%     of Means vs. pitch angle.  Only one such file is created.
%
%   - 3D*.eps: Same as *.fig, but an .eps image.
%



% OUTSIDE FUNCTIONS NEEDED
%
%  Extract_Filename
%
%  Extract_Extention 
%       [Only necessary for older versions of pitch.m]
%
%  fitsread
%
%  fitsheader
%
%  PeriodToDash



% TROUBLESHOOTING
%
% This code relies on fitsread.m to read the .fits files.  If fitsread.m
% reads the wrong pixel values, you might need to tell fitsread.m to 
% bitswap.  
%
% Fitsread.m contains a series of if/else statements that looks for the
% endianness of your architecture - that is, to determine whether 
% bitswapping is necessary.  The code isn't familiar with all 
% architectures.  For my Imac64, I had to add the following line to
% fitsread.m:
%
%       elseif strmatch(friend,'MACI64')
%         bswap = 'b';


% METHOD
%
% The code will record the pixel values along the spiral axes of several 
% spiral coordinate systems.  Then the code will measure the mean of all 
% pixel values along each coordinate axis.  The variance of these means 
% will be recorded for each coordinate system.  The true pitch angle will 
% yield the maximum variance of means.  The true pitch is passed as 
% PITCHvsINNER.


% MATH REVIEW
%
% Let a spiral be given by 
%
%       R = exp(B*THETA).  
%
% Then, using the astronomer's sign convention, the pitch angle is 
%
%       P = arctan(-B), and
%
% the arclength S from the origin to radius MSMT_OUTER is 
%
%       S = abs((MSMT_OUTER/B)*sqrt(1+B^2)).
%
% Source: Wolfram (http://mathworld.wolfram.com/LogarithmicSpiral.html)



%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% INITIAL HOUSEKEEPING
%
% Sometimes a research-weary user will swap the MINP and MAXP values.
% Instead of crashing, the code will correct the mistake.
%
 if MINP > MAXP     % If MINP is bigger than MAXP, swap the two
     TEMPP = MAXP;
     MAXP = MINP;
     MINP = TEMPP;
 end
%
%
% The code will tnterpret "InnerRadiusSpacing = 0" to mean the user only 
% wants to measure a single inner radius (specifically, MSMT_INNER1).
%
%
 if InnerRadiusSpacing==0       
     MSMT_INNER2=MSMT_INNER1;  
     InnerRadiusSpacing=1;     
 end
%
%
% If instructed by the user, we convert inputs into strings for use in
% output filenames.
% 
 if Save2D==1 || Save3D==1

    % Convert inputs to strings 
     
    [~,FILESTRING] = Extract_Filename(FILE);    % Chop path from filename
    FILESTRING2 = PeriodToDash(FILESTRING);     % Turn periods into dashes

    X0STRING = num2str(X0);                     % Convert inputs to strings
    X0STRING2 = PeriodToDash(X0STRING);

    Y0STRING = num2str(Y0);
    Y0STRING2 = PeriodToDash(Y0STRING);
    
    VIS_INNERSTRING = num2str(VIS_INNER);
    VIS_INNERSTRING2 = PeriodToDash(VIS_INNERSTRING);

    VIS_OUTERSTRING = num2str(VIS_OUTER);
    VIS_OUTERSTRING2 = PeriodToDash(VIS_OUTERSTRING);

    MSMT_INNER1STRING = num2str(MSMT_INNER1);
    MSMT_INNER1STRING2 = PeriodToDash(MSMT_INNER1STRING);
    
    MSMT_INNER2STRING = num2str(MSMT_INNER2);
    MSMT_INNER2STRING2 = PeriodToDash(MSMT_INNER2STRING);    

    InnerRadiusSpacingSTRING = num2str(InnerRadiusSpacing);
    InnerRadiusSpacingSTRING2 = PeriodToDash(InnerRadiusSpacingSTRING);
    
    MSMT_OUTERSTRING = num2str(MSMT_OUTER);
    MSMT_OUTERSTRING2 = PeriodToDash(MSMT_OUTERSTRING);

    NAXISSTRING = num2str(NAXIS);
    NAXISSTRING2 = PeriodToDash(NAXISSTRING);
    
    MINPSTRING = num2str(MINP);
    MINPSTRING2 = PeriodToDash(MINPSTRING);

    MAXPSTRING = num2str(MAXP);
    MAXPSTRING2 = PeriodToDash(MAXPSTRING);
    
    PSTEPSTRING = num2str(PSTEP);
    PSTEPSTRING2 = PeriodToDash(PSTEPSTRING);
           
    AxisPointSpacingSTRING = num2str(AxisPointSpacing);
    AxisPointSpacingSTRING2 = PeriodToDash(AxisPointSpacingSTRING);
           
    if SMOOTH==1 
        SMOOTHSTRING = '1';
    else SMOOTHSTRING = '0';
    end
    
    if Save2D==1
        Save2DString = 'Yes';  
    else
        Save2DString = 'No'; 
    end
    
    if Save3D==1
        Save3DString = 'Yes'; 
    else
        Save3DString = 'No'; 
    end
    
    % Create a general output string.  This string will be part of all
    % output filenames.
    
    GeneralOutString = [FILESTRING2 '_' X0STRING2 '_' Y0STRING2 '_' ...
        VIS_INNERSTRING2 '_' VIS_OUTERSTRING2 '_' MSMT_INNER1STRING2 ...
        '_' MSMT_INNER2STRING2 '_' InnerRadiusSpacingSTRING2 ...
        '_' MSMT_OUTERSTRING2 '_' NAXISSTRING2 '_' MINPSTRING2 '_' ...
        MAXPSTRING2 '_' PSTEPSTRING2 '_' AxisPointSpacingSTRING2 ...
        '_' SMOOTHSTRING];
    
    % Create the output filename for the textfile that summarizes the
    % results.
    
    TextFileName = ['Spirality_' GeneralOutString '.txt']; 
        
 end

 
%% Import data
%
% Read data from .fits file into local array
%
  FILEDATA = fitsread(FILE,'full',1); 

  
%% INITIALIZE SHIT & DO THE ONE-TIME COMPUTATIONS


% FITTING FUNCTION & DOMAINS  
% Our output will be the fitting function (Variance of Means) as a function 
% of pitch angle and inner radius.  So we first establish the two domains,
% then create an empty array for the fitting function.

% All inner measurement radii, in pixels.  The array InnerRadiusDomain has 
% dimension NumberOfInnerRadii.
  
  InnerRadiusDomain = MSMT_INNER1:InnerRadiusSpacing:MSMT_INNER2;           
  NumberOfInnerRadii = numel(InnerRadiusDomain);                            

% Since the user determines the inner radius spacing, the outermost inner 
% computation radius may be less than the outermost inneer radius 
% MSMT_INNER2 supplied by the user.  For example:
%
% Let MSMT_INNER1 = 0
% Let MSMT_INNER2 = 12
% Let InnerRadiusSpacing = 5
% Then InnerRadiusDomain = [0 5 10]
% 
% In this case, the outermost inner computation radius (10) is less than
% the outermost inner radius MSMT_INNER2 (12).  So we capture the outermost
% inner computation radius OutermostInnerRadius.  This array has dimension
% 1, and its value is in pixels.

  OutermostInnerRadius = InnerRadiusDomain(NumberOfInnerRadii);
  
% Pitch angles of all templates, in degrees and radians.  The arrays
% PitchAngleDomain and PitchAngleDomainInRadians have dimension
% NumberOfPitchAngleTemplates.

  PitchAngleDomain = MINP:PSTEP:MAXP;                     % Pitches in deg. 
  PitchAngleDomainInRadians = PitchAngleDomain*pi/180;    % Pitches in rad. 
  NumberOfPitchAngleTemplates = numel(PitchAngleDomain);  % Array size
 
% Compute the spiral constant for each pitch angle.  Spiral constants are 
% unitless.  These arrays have dimension NumberOfPitchAngleTemplates.

  NegPitches = -1*PitchAngleDomainInRadians;                                
  SpiralConstants = tan(NegPitches);                      % Spiral Constant 
  AbsSpiralConstants = abs(SpiralConstants);                                
  
% Computation time diverges as the template pitch angle approaches zero, 
% and a divide-by-zero error occurs when the pitch equals zero. 
%
% We avoid these problems by flagging zero-value or near-zero spiral 
% constants. 
% 
% For now, we set the tightest spiral we're willing to measure to be 0.01 
% degrees, or 0.00017 radians.  Due to the small angle approximation 
% (tanTheta ~= Theta), the spiral constant is nearly equal to the pitch 
% angle, with a negative sign, for very small angles.

  TightestPitchWeAreWillingToMeasureInDegrees = .01;  
  TightestSpiralWeAreWillingToMeasure = ...
      TightestPitchWeAreWillingToMeasureInDegrees*pi/180;                   

% Flag the unreasonably tight spirals.  This is an array of dimension 
% NumberOfPitchAngleTemplates.  Its values are 0 for unflagged elements,
% and 1 for flagged elements.

  TinySpiralConstants = AbsSpiralConstants < ...
      TightestSpiralWeAreWillingToMeasure;                                  

% For now, we just set the unreasonably tight spiral constants to some
% arbitrary value.  This will avoid the divide-by-zero errors in the next
% few computations.  Farther down the code, when the pitch angle templates
% are actually being computed, these unreasonably tight spirals will be
% ignored.
%
% This is an array of dimension NumberOfPitchAngleTemplates.  Its values
% are unitless.
  
  SpiralConstants(TinySpiralConstants) = 1;                                 
  AbsSpiralConstants(TinySpiralConstants) = 1;                              

% Compute the radius-to-arclength ratio for each pitch angle.  These arrays
% have dimension NumberOfPitchAngleTemplates.  Their values are unitless.
  
  Denominators = sqrt(1 + SpiralConstants.^2);                              
  RadiusToArclengthRatio = AbsSpiralConstants./Denominators;                

% Subtract the innermost inner radius from all inner radii.  This gives the
% distance from the innermost inner radius MSMT_INNER1 to all inner radii.
% This array has dimension NumberOfInnerRadii; its values are in pixels.

  DistancesFromInnermostRadiusToAllInnerRadii = ...
      InnerRadiusDomain - MSMT_INNER1;

% Arclengths, in pixels, from the innermost radius MSMT_INNER1 to each 
% inner radius.  This 2D array has dimensions (NumberOfInnerRadii, 
% NumberOfPitchAngleTemplates). 

  ArclengthsFromInnermostRadiusToAllInnerRadii = ...
      DistancesFromInnermostRadiusToAllInnerRadii' * ...
      RadiusToArclengthRatio.^-1;

% Transpose the array ArclengthsFromInnermostRadiusToAllInnerRadii.  Its
% new dimensions are (NumberOfPitchAngleTemplates, NumberOfInnerRadii).
  
  ArclengthsFromInnermostRadiusToAllInnerRadii = ...
    ArclengthsFromInnermostRadiusToAllInnerRadii';                          

% Index of each inner radius on the computation point array for each pitch 
% angle template.  The index depends on the pitch angle (the template) but 
% not on the phase angle (the individual axis).  This unitless 2D array has
% dimensions (NumberOfPitchAngleTemplates, NumberOfInnerRadii).

  IndicesOfInnerRadii = floor(...
      ArclengthsFromInnermostRadiusToAllInnerRadii / AxisPointSpacing + 1);

% Compute the arclength, in pixels, from the innermost measurement radius
% to the outer measurement radius for each pitch angle template.  These 
% arrays have dimension NumberOfPitchAngleTemplates.

  ArclengthsOriginToInner = MSMT_INNER1./RadiusToArclengthRatio;                       
  ArclengthsOriginToOuter = MSMT_OUTER./RadiusToArclengthRatio;             
  ArclengthsInnerToOuter = ArclengthsOriginToOuter - ...                    
      ArclengthsOriginToInner;               

% Compute the number of computation points per axis.  Equivalently, the 
% indices of the outer measurement radius MSMT_OUTER on the computation 
% point array for each pitch angle template.  This depends only on pitch 
% (the template), not on phase angle (the individual axis).  This array has 
% dimension NumberOfPitchAngleTemplates.

  NumbersOfPointsPerAxis = floor(ArclengthsInnerToOuter/ ...
      AxisPointSpacing +1);                                                 
  
% Maximum number of computation points for any axis.  This array has
% dimension 1.
  
  MaxNumberOfPointsPerAxis = max(NumbersOfPointsPerAxis);                     
 
% Initialize the fitting function.  This is a 2-D array with dimensions
% NumberOfInnerRadii by NumberOfPitchAngleTemplates. 

  FittingFunction = zeros(NumberOfInnerRadii,NumberOfPitchAngleTemplates);
  
% Phase angle, in radians, of each axis.  This array applies 
% to every pitch angle template.

  AxisPhaseAngles = 2*pi/NAXIS:2*pi/NAXIS:2*pi; 

% Initialize the radial coordinate, in pixels, at each computation point 
% for each pitch angle template.  This 2-D array has dimensions 
% (NumberOfPitchAngleTemplates, MaxNumberOfPointsPerAxis).

   RadiiOfComputationPoints = ...
      zeros(NumberOfPitchAngleTemplates, MaxNumberOfPointsPerAxis);
  
% OUTPUT ARRAY
% This array shows the peak pitch angle vs. inner radius.  It has dimension 
% (NumberOfInnerRadii,2) with the first column being the inner radius
% domain and the second column being the pitch peak angle.

  PITCHvsINNER = zeros(NumberOfInnerRadii,2);

% Populare the first column with all the inner measurement radii.
  
  PITCHvsINNER(:,1) = InnerRadiusDomain;    
  
  
  
  
  
%% DISPLAY PROGRESS
%
% Blow smoke up the user's ass.

 disp('Computing pitch angle templates.')


%% INITIALIZE PITCH ANGLE TEMPLATE ARRAY

% Now we can finally initialize the pitch angle template array.  This is a
% 3-D array with dimensions NumberOfPitchAngleTemplates by NAXIS by
% NumberOfPointsPerAxis.  The array elements are pixel values.

  PitchAngleTemplates = zeros(NumberOfPitchAngleTemplates,NAXIS,...
      MaxNumberOfPointsPerAxis);

%% COMPUTE PITCH ANGLE TEMPLATES
%
% We compute all the pitch angle templates first, then run computations on 
% them later.  Since computing the templates is the bulk of the 
% computation, this method significanly reduces computation time compared 
% to the old method, in which the templates were recomputed for each inner 
% radius.  
%
% The current method is as follows:
%   For each template
%       For each spiral axis
%           For each computation point
%               Locate the right pixel, record its value in the
%               PitchAngleTemplates array
%           end
%       end
%   end
%
% Shall we begin?


% For each pitch angle template

  for PitchCount = 1:NumberOfPitchAngleTemplates
      TemplatePitchAngle = PitchAngleDomainInRadians(PitchCount);
          
      % If the template pitch angle is too tight, emit a warning and ignore
      % the template.
      
      if TinySpiralConstants(PitchCount)==1
          
          
          % generate strings for the warning message.
              
            StringTightestPitchWeAreWillingToMeasureInDegrees = ... 
               num2str(TightestPitchWeAreWillingToMeasureInDegrees);
              
            % Convert template pitch to degrees so the humans will
            % understand.  Then generate a string.
              
            TemplatePitchAngleInDegrees = TemplatePitchAngle*180/pi;
              
            StringTemplatePitchAngleInDegrees = ...
                num2str(TemplatePitchAngleInDegrees);
                        
            % Spit out the warning.
            
            warning(['Templates tighter than ' ...
                StringTightestPitchWeAreWillingToMeasureInDegrees ...
                ' degrees are not computed.  Template of pitch '...
                StringTemplatePitchAngleInDegrees ' degrees has ' ....
                'been ignored.']) 
          

      % Else, if the template pitch is reasonable, then proceed with 
      % computation.  

      else
          
          % Ratio of radius to arclength for this pitch angle template.
          % This is a unitless array of dimension 1.
          
          RadiusToArclengthForThisPitch = ...
              RadiusToArclengthRatio(PitchCount);
          
          % Spiral constant for this pitch angle template.
          % This is a unitless array of dimension 1.
          
          SpiralConstantForThisTemplate = SpiralConstants(PitchCount);        
          
          % Number of computation points per axis for this template
          % This is a unitless array of dimension 1.
          
          NumberOfComputationPointsPerAxisForThisTemplate = ...
              NumbersOfPointsPerAxis(PitchCount);  
          
          % Arclength, in pixels, from innermost radius to outer radius for 
          % this template.  This is an array of dimension 1.
          
          ArclengthInnerToOuterForThisPitch = ...
              ArclengthsInnerToOuter(PitchCount);
          
          % Arclengths, in pixels, from the innermost radius to each 
          % computation point.  This is an array of dimension
          % NumberOfComputationPointsPerAxisForThisTemplate.
         
          ArclengthsFromInnerToComputationPoint = 0:AxisPointSpacing: ...
              ArclengthInnerToOuterForThisPitch;
          
          % Radial distance, in pixels, from innermost radius to each 
          % computation point.  This is an array of dimension
          % NumberOfComputationPointsPerAxisForThisTemplate.
          
          RadialDistances = ArclengthsFromInnerToComputationPoint * ... 
              RadiusToArclengthForThisPitch;                               
          
          % Radial coordinate, in pixels, from spiral center to each 
          % computation point.  This is an array of dimension
          % NumberOfComputationPointsPerAxisForThisTemplate.
          
          RadiiOfComputationPointsForThisTemplate = ...
              RadialDistances + MSMT_INNER1;                               

          % The azimuthal angle of the dead center of a spiral is -inf. The
          % code is happy to tell you that.  However, later on, if it tried
          % to use -inf to compute things, it would get cranky.  So, in
          % order to avoid such errors, we check to see if the inner
          % measurement radius is zero (which would cause the first azimuth
          % to be -Inf).  If so, we offset the innermost measurement point
          % by exp(-10) pixels.
          
          if MSMT_INNER1==0
              RadiiOfComputationPointsForThisTemplate(1) = ...
                  exp(-10);
          end

          
          % Store RadiiOfComputationPointsForThisTemplate in the 2-D array 
          % RadiiOfComputationPoints.  The 2-D array has dimensions 
          % (NumberOfPitchAngleTemplates, MaxNumberOfPointsPerAxis).  Its
          % values' units are pixels.
                    
          RadiiOfComputationPoints(PitchCount,1: ...
              NumberOfComputationPointsPerAxisForThisTemplate) = ...
              RadiiOfComputationPointsForThisTemplate;                     
          
                   
          % Azimuthal angle (Polar coordinate theta), in radians, of each
          % computation point, notwithstanding the phase angle of each 
          % axis.  This is an array of dimension
          % NumberOfComputationPointsPerAxisForThisTemplate.
          
          AzimuthsWithoutPhaseAngles = ...
              log(RadiiOfComputationPointsForThisTemplate)/...
              SpiralConstantForThisTemplate;                               
          
          % For each spiral axis
              
          for AxisCount = 1:NAXIS
              
              % Phase angle, in radians, of this particular axis.  This is
              % an array of dimension 1.
              
              PhaseAngle = AxisPhaseAngles(AxisCount);                     
              
              % Azimuthal angle (Polar coordinate theta), in radians, of
              % each computation point, now taking into account the phase
              % angle of this particular axis.  This is an array of
              % dimension NumberOfComputationPointsPerAxisForThisTemplate.
              
              Azimuths = AzimuthsWithoutPhaseAngles - PhaseAngle;          
              
              % X-coordinate, rounded to the nearest pixel, of each 
              % computation point.  This is an array of dimension
              % NumberOfComputationPointsPerAxisForThisTemplate.
              
              Xpixels = round(X0 + ...
                  RadiiOfComputationPointsForThisTemplate.*cos(Azimuths)); 
                  

              % Y-coordinate, rounded to the nearest pixel, of each 
              % computation point.  This is an array of dimension
              % NumberOfComputationPointsPerAxisForThisTemplate.
              
              Ypixels = round(Y0 + ...
                  RadiiOfComputationPointsForThisTemplate.*sin(Azimuths)); 
                  
              
              % For each computation point on this axis
              
              for ComputationPointCount = 1:...
                      NumberOfComputationPointsPerAxisForThisTemplate
                  
                  % X-value of the pixel at this particular computation
                  % point.  This is an array of dimension 1.
                  
                  XX = Xpixels(ComputationPointCount);                     
                  
                  % Y-value of the pixel at this particular computation
                  % point.  This is an array of dimension 1.

                  YY = Ypixels(ComputationPointCount);                     
                  
                 % Record the pixel value at the computation point in 
                 % the array PitchAngleTemplates.
                  
                  PitchAngleTemplates(PitchCount,AxisCount,...
                      ComputationPointCount) = FILEDATA(YY,XX);               
                                               
              end % End this computation point
          end     % End this spiral axis
      end         % End the condition that the pitch angle be reasonable
  end             % End this pitch angle template
      

 disp('Templates computed.')
 disp('')              
 disp(['NumberOfInnerRadii = ' num2str(NumberOfInnerRadii)])
 disp('')
 disp('  #     INNER  PITCH')  


 
 
 
 
%%%% COMPUTE VARIANCE OF MEANS (i.e., the fitting function).
%
% PSEUDO-CODE
%
% For each template
%
%   Establish an array (MeanPixelValuesForEachAxis) with the mean pixel
%   value of each axis for each inner radius.  This array will have
%   dimension NAXIS.
%      
%
%   For each inner radius (from OutermostInnerRadius backward to
%   MSMT_INNER1)
%
%       Compute the mean pixel value for each axis.
%       Take the variance of all the means; record it in the fitting
%       function
%
%   End this inner radius
%
% End this template       
%



% REAL CODE
%
% For each pitch angle template

  for PitchCount2 = 1:NumberOfPitchAngleTemplates

      % If the pitch is loose enough to measure
      
      if TinySpiralConstants(PitchCount)==0
          
      % Number of computation points per axis for this template
      
      NumberOfComputationPointsPerAxisForThisTemplate = ...                 
              NumbersOfPointsPerAxis(PitchCount2);
                           
      % Establish a 2-D array with the sum of pixel value of each axis for 
      % each inner radius.  This array has dimension NAXIS.
            
      SumOfPixelValuesForEachAxis = zeros(1,NAXIS);
      
      % We save computation time by being clever.  Start by finding the sum
      % of the spiral segment from the outermost inner radius MSMT_INNER2
      % to the outer radius MSMT_OUTER.  Compute the mean.  Then we move
      % inward, to the next inner radius.  Instead of summing all those
      % numbers again, we simply add the next bit of spiral segment to the
      % previously computed sum.  Compute the new mean.  Etc.
      %
      % We start by establishing boundary markers for the spiral segment
      % we're about to sum.  The boundary markers each have dimension 1.
      %
      % Initially, inner boundary is the computation point element 
      % corresponding to the outermost inner radius OutermostInnerRadius.  
      % The outer boundary corresponds to the outer radius MSMT_OUTER.  
           
      
      SpiralSegmentOuterBoundary = ...
          NumberOfComputationPointsPerAxisForThisTemplate;

      SpiralSegmentInnerBoundary = ...                                      
          IndicesOfInnerRadii(PitchCount2,NumberOfInnerRadii);                  
      
      % Compute the number of computation points per axis for this inner
      % radius (initially, for MSMT_INNER2).
      
      NumberOfComputationPointsForThisInnerRadius = ...
          SpiralSegmentOuterBoundary - SpiralSegmentInnerBoundary + 1;
      
      % For each inner radius
      
      for InnerRadiusCount = NumberOfInnerRadii:-1:1
          
          % Compute the sum of pixel values along the spiral segment; add
          % to the previous sum.  
          % 
          % Syntax note: sum(blah,3) means we're summing over the 3rd 
          % dimension of blah.  In this case, we're summing over 
          % the 3rd dimension - the computation points - of the 
          % PitchAngleTemplates array.  Note that we're holding the first
          % dimension - the pitch angle template - constant.
          %
          % These arrays have dimension NAXIS.
          
          
          ThingsToAddToSumOfPixelValuesForEachAxis = ...
              sum(PitchAngleTemplates(PitchCount2,:,...
              SpiralSegmentInnerBoundary:SpiralSegmentOuterBoundary),3);

          SumOfPixelValuesForEachAxis = ...
              SumOfPixelValuesForEachAxis + ...
              ThingsToAddToSumOfPixelValuesForEachAxis;
          
          % Compute the mean pixel value for each axis.  This array has
          % dimension NAXIS.
          
          MeanPixelValuesForEachAxis = SumOfPixelValuesForEachAxis / ...
              NumberOfComputationPointsForThisInnerRadius;
          
          % Compute the variance of the mean pixel values.  This array has
          % dimension 1.
          
          VarianceOfMeans = var(MeanPixelValuesForEachAxis,1);
          
          % Record the result in the fitting function.
          
          FittingFunction(InnerRadiusCount,PitchCount2) = ...
              VarianceOfMeans;
          
          % Move outer boundary marker to just inside the inner boundary.
          
          SpiralSegmentOuterBoundary = SpiralSegmentInnerBoundary - 1;
          
          % If we're not already at the innermost inner radius, move the
          % inner boundary marker inward to the next inner radius.
          
          NextInnerRadiusCount = InnerRadiusCount - 1;
          
          if NextInnerRadiusCount > 0
              SpiralSegmentInnerBoundary = ...
                    IndicesOfInnerRadii(PitchCount2,NextInnerRadiusCount);
          end
          
          
      end     % End this inner radius
      end     % End the condition that the pitch angle be reasonable
  end         % End this pitch angle template
      
      
  

% Smooth the fitting function along the pitch angle axis, if the user 
% requested it.  Otherewise, let the "smoothed" fitting function equal the
% original fitting functino.  This is a 2-D array with dimensions
% NumberOfInnerRadii by NumberOfPitchAngleTemplates. 

SmoothFittingFunction = FittingFunction;

if SMOOTH==1
    for InnerRadiusCount2 = 1:NumberOfInnerRadii
          SmoothFittingFunction(InnerRadiusCount2,:) = ...
              smooth(FittingFunction(InnerRadiusCount2,:));
    end
end


%                             -----


%%%%%%%%%%%%%%%%%%%%%%% Output: PITCHvsINNER %%%%%%%%%%%%%%%%%%%%%%% 

% Now that the fitting function is complete, we analyze the results.  We
% look at each row (each inner radius) of the fitting function, and find
% the peak value.

% For each inner radius

for InnerRadiusCount3 = 1:NumberOfInnerRadii
    
    % Import the data from the fitting function, both original and smooth.  
    % Specifically, import the slice of each fitting function that 
    % corresponds to this inner radius.  Each of these arrays has dimension 
    % NumberOfPitchAngleTemplates.
    
    FitForThisInnerRadius = FittingFunction(InnerRadiusCount3,:);
    SmoothFitForThisInnerRadius = ...
        SmoothFittingFunction(InnerRadiusCount3,:);

    % For each inner radius, we index the peak of the fitting function.  
    % The index will correspond to the pitch angle template that 
    % corresponds to the strongest pitch angle for this inner radius.  The
    % index array has dimension 1, unless multiple templates show equally 
    % maximum fits.
    
    WhereThePeakIs=find...
       (SmoothFitForThisInnerRadius==max(SmoothFitForThisInnerRadius));   
      
    % If multiple templates show equally equally maximum fits, the above
    % index array will have dimension >1.  We then take the mean of the 
    % indices and round.  The array will then have dimension 1.
      
    if numel(WhereThePeakIs) > 1
        WhereThePeakIs=round(mean(WhereThePeakIs));
    end
      
    % Capture the pitch angle of the template corresponding to the peak 
    % for this inner radius. 
    
    PeakPitchForThisInnerRadius = ...
         PitchAngleDomain(WhereThePeakIs);
    
    % Record the measured pitch angle in the output array PITCHvsINNER.
    % The output array PITCHvsINNER has dimensions (NumberOfInnerRadii,2).
      
    PITCHvsINNER(InnerRadiusCount3,2)=PeakPitchForThisInnerRadius;
        
    % Establish the inner radius, in pixels.  This array has dimension 1.
    
    ThisInnerRadius = InnerRadiusDomain(InnerRadiusCount3);
     
    % If the user requested, generate a 2-D graph of the fitting function
    % vs. pitch angle for this inner radius.
    
    if Save2D == 1
        
        % Convert this inner radius to a string
                
        ThisInnerRadiusString = num2str(ThisInnerRadius);

        ThisInnerRadiusString2 = PeriodToDash(ThisInnerRadiusString);
        
        % Generate the output filename for this graph.
        
        FileName2D = ['Spirality2D_' GeneralOutString '_' ...
            ThisInnerRadiusString2];
        
        % Generate the graph, save as both a MATLAB figure (.fig) and an 
        % encapsulated post scrips (.eps), 
        
        h2 = figure('NumberTitle','off',...
        'Visible','off'); 
        plot(PitchAngleDomain,FitForThisInnerRadius,...
            PitchAngleDomain,SmoothFitForThisInnerRadius,'Visible','on')
        xlabel('Pitch (degrees)')
        ylabel('Fitting Function (arbitrary units)')
        saveas(h2,FileName2D,'fig')
        saveas(h2,FileName2D,'eps')
        close(h2)
        

    end % End the 2-D graph for this inner radius
         
    %% Display the inner radius number & Inner Radius vs. Pitch on screen

    % PeakPitchString = num2str(PeakPitchForThisInnerRadius);    
    format =  '%3d   %6.2f   %6.2f';
    [s,errmsg] = sprintf(format,InnerRadiusCount3,...
        ThisInnerRadius,PeakPitchForThisInnerRadius); %#ok<NASGU>
    disp(s)
    
end % End this inner radius

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%                             -----

%%%%%%%%%%%%%%%%%%%%%%% Output: BESTFITPITCH %%%%%%%%%%%%%%%%%%%%%%% 
%
% This is the best-fit pitch angle for the spiral.  It is the mean of the 
% strongest pitch angles measured for each inner radius.  in the second 
% column of PITCHvsINNER.

BESTFITPITCH = mean(PITCHvsINNER(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                             -----

%%%%%%%%%%%%%%%%%%%%%%%%%%% Output: ERROR %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The error for the best fit pitch is computed in three steps.  The first
% step is the random error, which is the standard deviation of the pitch 
% angles given in the second column of PITCHvsINNER.

RANDERR = std(PITCHvsINNER(:,2));  % Random error        

% The error is then penalized, or "stretched", if the measurement is based
% on a small inner radius annulus.  It is stretched by the ratio the radius 
% segments of two annuli: the visible spiral annulus to the inner radius 
% annulus.  
%
% On the rare occasion that the inner radius annulus is larger than the
% visible spiral annulus, the code will not compress the error bars.
% Rather, the code sets the stretching coefficient to 1.

STRETCHERR = RANDERR*max(1,(VIS_OUTER - VIS_INNER)/ ...  % Stretched error
    (OutermostInnerRadius - MSMT_INNER1));                                                        

% Finally, the stretched error is added in quadrature with PSTEP, which the
% distance in degrees between successive pitch angle templates.  PSTEP is
% the precision of the measurement, and is therefore the smallest possible
% error bar.

ERROR = sqrt(STRETCHERR^2 + PSTEP^2);                      % Total error

% Check to see whether the precision error PSTEP conributed significantly
% to the total error.  If so, throw a warning.
%
% For now, we define the maximum acceptible error contribution by PSTEP to 
% be 1%.  In future versions, we may ask the user to input the maximum
% acceptable error contribution, or else include the contribution as an
% output.

MaximumAcceptibleErrorContributionByPSTEP = .01;
ErrorContributionByPSTEP = (ERROR - STRETCHERR)/STRETCHERR;

if ErrorContributionByPSTEP > MaximumAcceptibleErrorContributionByPSTEP
    
    % Convert error contribution to percent
    
    ErrorContributionByPSTEP = 100*ErrorContributionByPSTEP;
    
    % Convert error contribution to string
    
    StringErrorContributionByPSTEP = ...
        num2str(ErrorContributionByPSTEP);
    
    % Throw a warning
    
    warning(['The precision error PSTEP added ' ...
        StringErrorContributionByPSTEP ...
        ' percent to your total error.']) %#ok<*WNTAG>
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                             -----


% If the user asked for it, save the results in a text file.

if Save2D == 1 || Save3D == 1
    
    % Open a new text file.
    
    fid = fopen(TextFileName,'wt');
    
    % Write the header text
    
    fprintf(fid, 'PITCH ANGLE MEASUREMENT \n');
    fprintf(fid, 'Angles in degrees, radii in pixels \n \n');
    fprintf(fid, 'INPUTS: \n');

    % Write the inputs
    
    formatSpec = '    File                    %s \n';
    fprintf(fid,formatSpec,FILESTRING);

    formatSpec = '    X0                      %s \n';
    fprintf(fid,formatSpec,X0STRING);

    formatSpec = '    Y0                      %s \n';
    fprintf(fid,formatSpec,Y0STRING);

    formatSpec = '    VIS_INNER               %s \n';
    fprintf(fid,formatSpec,VIS_INNERSTRING);
    
    formatSpec = '    VIS_OUTER               %s \n';
    fprintf(fid,formatSpec,VIS_OUTERSTRING);

    formatSpec = '    MSMT_INNER1             %s \n';
    fprintf(fid,formatSpec,MSMT_INNER1STRING);

    formatSpec = '    MSMT_INNER2             %s \n';
    fprintf(fid,formatSpec,MSMT_INNER2STRING);

    formatSpec = '    InnerRadiusSpacing      %s \n';
    fprintf(fid,formatSpec,InnerRadiusSpacingSTRING);

    formatSpec = '    MSMT_OUTER              %s \n';
    fprintf(fid,formatSpec,MSMT_OUTERSTRING);
    
    formatSpec = '    NAXIS                   %s \n';
    fprintf(fid,formatSpec,NAXISSTRING);

    formatSpec = '    MINP                    %s \n';
    fprintf(fid,formatSpec,MINPSTRING);

    formatSpec = '    MAXP                    %s \n';
    fprintf(fid,formatSpec,MAXPSTRING);

    formatSpec = '    PSTEP                   %s \n';
    fprintf(fid,formatSpec,PSTEPSTRING);

    formatSpec = '    AxisPointSpacing        %s \n';
    fprintf(fid,formatSpec,AxisPointSpacingSTRING);

    formatSpec = '    SMOOTH                  %s \n';
    fprintf(fid,formatSpec,SMOOTHSTRING);

    formatSpec = '    2-D Graphs Saved        %s \n';
    fprintf(fid,formatSpec,Save2DString);

    formatSpec = '    3-D Graph Saved         %s \n';
    fprintf(fid,formatSpec,Save3DString);

    fprintf(fid, '\n');
    
    % Write the individual results, pitch vs. inner radius
    
    fprintf(fid, 'RESULTS \n');
    fprintf(fid, '    INNER       PITCH \n');
    formatSpec = '    %9.5f   %9.5f \n';
    for i = 1:NumberOfInnerRadii
        fprintf(fid,formatSpec,PITCHvsINNER(i,1),PITCHvsINNER(i,2));
    end
    fprintf(fid, ' \n');

    % Write the combined results for all inner radii
    
    formatSpec = 'Best Fit Pitch      %9.5f \n';
    fprintf(fid,formatSpec,BESTFITPITCH);
    formatSpec = 'Random Error        %9.5f \n';
    fprintf(fid,formatSpec,RANDERR);
    formatSpec = 'Stretched Error     %9.5f \n';
    fprintf(fid,formatSpec,STRETCHERR);
    formatSpec = 'Total Error         %9.5f \n';
    fprintf(fid,formatSpec,ERROR);

    % Close the text file
    
    fclose(fid);

end
    

% 3-D graph
%

if Save3D==1
    
    % Generate the 3D file string
    
    FileName3D = ['Spirality3D_' GeneralOutString];
    
    % Generate the 3D graph
    
        h4 = figure('NumberTitle','off',...
        'Visible','off'); 
        contourf(PitchAngleDomain,InnerRadiusDomain,...
            SmoothFittingFunction,'Visible','on')
        xlabel('Pitch (degrees)')
        ylabel('Inner Radius (pixels)')
        saveas(h4,FileName3D,'fig')
        saveas(h4,FileName3D,'eps')
        close(h4)
        
end


return
end

