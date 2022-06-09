function ...
    [OutputArray] = SpiralArmCount(FILE,X0,Y0,INNER,OUTER,PITCH)

% Written 2015 
% Doug Shields & Ben Davis
% University of Arkansas

%%%%% INTRO
% This function counts the spiral arms in a spiral image.  The input file
% muse be in .fits format.  The pixel location (X0, Y0) of the spiral
% center, the inner (INNER) and outer (OUTER) radii of the measurement
% annuli, and the spiral's pitch angle must be specified.

% The output are two plots: median pixel value vs. phase angle and that
% function's FFT.  Each of those plots is output in both .fit and .eps
% format, for a total of 4 output files.  

% By counting the bumps on the median pixel value plot, the user can
% visually determine the number of spiral arms.  Alternately, the FFT plot
% will count the spiral arms directly.

% The algorithm is the following: First the function computes a spiral 
% coordinate system (template) with the given pitch angle.  It computes the 
% median pixel value along each axis.  It then outputs the mean pixel value 
% vs. Phase angle, resulting in one peak for each spiral arm.  Finally, it
% computes the FFT (power vs. mode), the peak of which occurs at the true 
% number of spiral arms.

%%%%% INPUTS

% - FILE (string) is the name of the image file.  It must be a .fits image.

% - X0 and Y0 (doubles) are the x- and y-pixels of the galaxy's center.

% - INNER (double) is the inner radius of the measurement annulus, in 
%   pixels.  It should be just outside the bar or bulge.

% - OUTER (double) is the outer radius of the measurement annulus, in 
%   pixels.  It should coincide with the outer radius of the galaxy.

% - PITCH (double) is the pitch angle, in degrees, of the galaxy.  For  
%   S-windings, the pitch is positive.  For z-windings, the pitch is 
%   negative.  PITCH should not be zero (or 2pi times any integer), lest
%   computation time diverge.

%%%%% ARRAY OUTPUT

% - MedianVsPhaseAngle (double, 2 x NumberOfAxes) is a two-row array with 
%   as many columns as axes in the spiral coordinate system.
%   The first row shows the sum of pixel values along each spiral axis.
%   The second row shows the phase angle, in radians, of each spiral axis.

%%%%% FILE OUTPUTS

% - Median pixel value vs. phase angle (.fig)
% - Median pixel value vs. phase angle (.eps)
% - FFT (Power vs. Mode) (.fig)
% - FFT (Power vs. Mode) (.eps)

%%%%% OUTSIDE FUNCTIONS NEEDED

% - fitsread
% - PeriodToDash
% - num2str




%%%%% PRELIMINARY COMPUTATIONS

UseMedians=0;   % Toggle for using mean pixel value (UseMedian=0) or
UseMedians=1;   % median pixel value (UseMedian=1).      

FILEDATA = fitsread(FILE,'full',1);   % Read data from .fits file

NumberOfAxes = round(4*pi*OUTER)+1;   % (integer) Number of axes in the 
                                      % spiral coordinate system
                                  
AxisPointSpacing = 0.1;           % (double) Number of pixels between 
                                  % successive computation points on each 
                                  % spiral axis.
                                  
PitchInRadians = PITCH*pi/180;      % (double) Convert pitch to radians.

SpiralConstant = ...         % (double) Spiral Constant for this pitch.     
    -1*tan(PitchInRadians);  % A logarithmic spiral is defined by
                             %    R = exp(SpiralConstant*Theta)
                             % and the pitch angle is
                             %    PITCH = arctan(-1*SpiralConstant)

RadiusToArclengthRatio = ...      % (double) Ratio of the radial coordinate 
    (abs(SpiralConstant))/ ...    % to the arclength of the spiral.  A 
    sqrt(1 + SpiralConstant^2);   % useful constant for later computations.

ArclengthOriginToInner = ...      % (double) Spiral axis arc length, in 
    INNER/RadiusToArclengthRatio; % pixels, from the origin to the inner  
                                  % radius of the measurement annulus.

ArclengthOriginToOuter = ...      % (double) Spiral axis arc length, in 
    OUTER/RadiusToArclengthRatio; % pixels, from the origin to the outer  
                                  % radius of the measurement annulus.

ArclengthInnerToOuter = ...       % (double) Spiral axis arc length, in 
    ArclengthOriginToOuter - ...  % pixels, from the inner radius to the   
    ArclengthOriginToInner;        % outer radius.

NumberOfPointsPerAxis = ...           % Number of computation points
    floor(ArclengthInnerToOuter/ ...   % on each spiral axis
    AxisPointSpacing +1);                                                 



%%%%% COMPUTATIONS COMMON TO ALL SPIRAL AXIS

ArclengthsInnerToComputationPoint = ...  % (double, 
     0:AxisPointSpacing: ...             % 1 x NumberOfPointsPerAxis)
     ArclengthInnerToOuter;              % Spiral arc length, in pixels,
                                         % from the inner radius of the
                                         % measurement annulus to each
                                         % computation point on the spiral
                                         % axis.
                                         
 RadialDistances = ...
     ArclengthsInnerToComputationPoint * ... 
     RadiusToArclengthRatio;             % (double, 
                                         % 1 x NumberOfPointsPerAxis)
                                         % Radial distance, in pixels,
                                         % from the inner radius of the
                                         % measurement annulus to each
                                         % computation point on the spiral
                                         % axis.
          
RadiiOfComputationPoints = ...           % (double,
    RadialDistances + INNER;       % 1 x NumberOfPointsPerAxis) 
                                         % Radial coordinate, in pixels, of
                                         % each computation point.
                                         
if INNER==0                                  % Avoid divide-by-zero errors
    RadiiOfComputationPoints(1) = exp(-10);  % by moving the inner radius
end                                          % of the measurement annulus                                            
                                             % slightly outward if it is
                                             % zero.
        
AzimuthsWithoutPhaseAngles = ...         % (double,
    log(RadiiOfComputationPoints)/...    % 1 x NumberOfPointsPerAxis)
    SpiralConstant;                      % Azimuthal angle (Polar 
                                         % coordinate theta), in radians, 
                                         % of each computation point, 
                                         % notwithstanding the phase angle 
                                         % of the axis. 
          
          
                                         
 %%%%% COMPUTATIONS FOR INDIVIDUAL SPIRAL AXES

PitchAngleTemplate = ...            % (double, NumberOfAxes x 
    zeros(NumberOfAxes,...          % NumberOfPointsPerAxis)
    NumberOfPointsPerAxis);         % Pixel value at every computation 
                                    % point along every spiral axis
                                                                                                       
Means = zeros(1,NumberOfAxes);      % (double, 1 x NumberOfAxes) Mean pixel
                                    % value for each spiral axis
                        
Medians = zeros(1,NumberOfAxes);    % (double, 1 x NumberOfAxes) Median 
                                    % pixel value for each spiral axis
                                    
                                    
PhaseAngles = ...                   % (double, 1 x NumberOfAxes) Phase  
    linspace(0,2*pi,NumberOfAxes);  % angle, in radians, of each spiral
                                    % axis.

for AxisCount = 1:NumberOfAxes      % For each spiral axis
    
    PhaseAngle = PhaseAngles(AxisCount);  % (double) Phase angle, in 
                                          % radians, of this 
                                          % particular axis.
    
    Azimuths = ...                        % (double, 
        AzimuthsWithoutPhaseAngles - ...  % 1 x NumberOfPointsPerAxis)
        PhaseAngle;                       % Azimuthal angle (Polar 
                                          % coordinate theta), in radians, 
                                          % of  each computation point, now 
                                          % taking into account the phase
                                          % angle of this particular axis.
                                          
     Xpixels = round(X0 + ...             % (integer,
         RadiiOfComputationPoints.* ...   % 1 x NumberOfPointsPerAxis)
         cos(Azimuths));                  % X-coordinate, rounded to the 
                                          % nearest pixel, of each
                                          % computation point.
                                          
     Ypixels = round(Y0 + ...             % (integer,
         RadiiOfComputationPoints.* ...   % 1 x NumberOfPointsPerAxis)
         sin(Azimuths));                  % Y-coordinate, rounded to the 
                                          % nearest pixel, of each
                                          % computation point.
                                          
      for ComputationPointCount = 1:...   % For each computation point on
              NumberOfPointsPerAxis       % this axis   
          
          XX = ...                             % (integer) X-coordinate
              Xpixels(ComputationPointCount);  % of the pixel at this
                                               % computation point.
                                               
          YY = ...                             % (integer) Y-coordinate
              Ypixels(ComputationPointCount);  % of the pixel at this
                                               % computation point.
                                               
          PitchAngleTemplate(AxisCount,...     % (double) Pixel value at 
              ComputationPointCount) = ...     % this computation point
              FILEDATA(YY,XX);                 % on this spiral axis

      end                                      % End this computation point
      
      AxisSum =  ...                      % (double) Sum of pixel values at 
          sum(PitchAngleTemplate ...      % all computation points along
          (AxisCount,:));                 % this spiral axis 
                                                                                              
      Means(1,AxisCount) = ...            % (double) Mean pixel value along
          AxisSum*AxisPointSpacing;       % this spiral axis
      
      Medians(1,AxisCount) = ...          % (double) Mean pixel value at
          median(PitchAngleTemplate...    % all computation points along
          (AxisCount,:));                 % % this spiral axis
      
        
end                                       % End this spiral axis

Means = Means/max(Means);                 % Normalize the means

Medians = Medians/max(Medians);           % Normalize the medians
%Medians = Medians*-1;                     % Flip the function about y=0
%Medians = Medians + 1 - max(Medians);     % Increase so that the max is 1


%%%%% OUTPUT FILENAMES
 
[~,FILESTRING] = Extract_Filename(FILE);    % Chop path from filename
FILESTRING2 = PeriodToDash(FILESTRING);     % Turn periods into dashes

X0STRING = num2str(X0);                     % Convert inputs to strings
X0STRING2 = PeriodToDash(X0STRING);

Y0STRING = num2str(Y0);
Y0STRING2 = PeriodToDash(Y0STRING);
    
INNERSTRING = num2str(INNER);
INNERSTRING2 = PeriodToDash(INNERSTRING);

OUTERSTRING = num2str(OUTER);
OUTERSTRING2 = PeriodToDash(OUTERSTRING);

PITCHSTRING = num2str(PITCH);
PITCHSTRING2 = PeriodToDash(PITCHSTRING);
    
GeneralOutString =  ...                     % (string) Create a general
    [FILESTRING2 '_'...                     % string to be used in the 
    X0STRING2 '_' Y0STRING2 '_' ...         % output filenames for both the
    INNERSTRING2 '_' OUTERSTRING2 '_' ...   % .fig file and the .eps file.
    PITCHSTRING2];                          

MedianPixelOutString = ...                  % Filename of median pixel             
    ['SpiralArmCount_MedianPixel_' ...      %   value vs. phase angle graph
    GeneralOutString];

FFTOutString = ...                           % Filename of FFT graph
    ['SpiralArmCount_FFT_' ...
    GeneralOutString];


TITLESTRING= ...                            % (string) Title of graphs
    ['Pitch = ' PITCHSTRING ' deg.'];



%%%%% OUTPUT ARRAY

MeanVsPhaseAngle = [PhaseAngles;Means];   % (double, 2 x NumberOfAxes)
                                          % Mean pixel value (Row 1) vs. 
                                          % phase angle (Row 2) for each 
                                          % spiral axis
                                        
MedianVsPhaseAngle = [PhaseAngles;Medians];   % (double, 2 x NumberOfAxes)
                                          % Mean pixel value (Row 1) vs. 
                                          % phase angle (Row 2) for each 
                                          % spiral axis
                                        
if UseMedians==1
    OutputArray=MedianVsPhaseAngle;
else
    OutputArray=MeanVsPhaseAngle;
end
                                          
%%%%% OUTPUT FIGURE: MEAN/MEDIAN PIXEL VALUE VS. PHASE ANGLE

h2 = figure...                            % Generate a figure
    ('NumberTitle','off',...              % Don't number the figure
    'Visible','off');                     % Don't display on screen

axes1 = axes('Parent',h2,...              % Set axis properties
    'XTick',[0 pi/2 pi 3*pi/2 2*pi],...   % x-axis tick mark locations
    'XTickLabel',...                      % x-axisick mark labels.  We set 
    {'0','pi/2','pi','3pi/2','2pi'});
xlim(axes1,[0 2*pi]);                     % Set domain


set(gca,'FontSize',20)                    % Font size for axis tick marks

box(axes1,'on');                          % Draw axes
hold(axes1,'all');                        % Keep axis properties when 
                                          %   plotting
if UseMedians==1
  plot(PhaseAngles,Medians,...            % Plot Mean pixel value vs. 
    ...                                   %   Phase Angle
    'Visible','on',...                    % Draw plot in figure.
    'LineWidth',3, ...                    % Plot line is thick.
    'Color',[0 0 0]);                     % Plot line is black.
  ylabel...                                 
    ('Median relative pixel value (arbitrary units)',...  % y-labels 
    'FontSize',15,'fontname','helvetica')                 % and title
else
  plot(PhaseAngles,Means,...              % Plot Mean pixel value vs. 
    ...                                   %   Phase Angle
    'Visible','on',...                    % Draw plot in figure.
    'LineWidth',3, ...                    % Plot line is thick.
    'Color',[0 0 0]);                     % Plot line is black.
  ylabel...                                 
    ('Mean relative pixel value (arbitrary units)',...
    'FontSize',15,'fontname','helvetica') 
end



xlabel('Phase angle','FontSize',15,...    % x-labels and title. 
    'fontname','helvetica')              
          
title(TITLESTRING,'FontSize',15,...       
    'fontname','helvetica');

saveas(h2,MedianPixelOutString,'fig')         % Save a MATLAB figure
saveas(h2,MedianPixelOutString,'eps')         % Save as .eps 

close(h2)                                 % Close figure



% TAKE FFT of result.  In other words, count the number of spiral arms.

if UseMedians==1
    Y = fft(Medians);
else
    Y = fft(Means);
end

Y(1) = [];           % Chop first element from fft result (why?)
n = length(Y);       % Count the elements in the fft array
N = n/(2*pi);        % Count the number of whole-period cycles 
power = abs(Y(1:floor(n/2))).^2;     % Strength of each mode (Xi* Xi)
power = power/max(power);            % Scale down to unity 
nyquist = 1/2;                       % Nyquist frequency
freq = (1:n/2)/(n/2)*nyquist*N;      % Frequency of each fft data point
period = 1./freq;                    % Period of each fft data point
m = (2*pi)./period;                  % Mode (number of arms)

%%%%% OUTPUT FIGURE: STRENGTH VS. MODE


h3 = figure...                            % Generate a figure
    ('NumberTitle','off',...              % Don't number the figure
    'Visible','off');                     % Don't display on screen
plot(m,power,...                          % Plot strength vs. mode 
    'Visible','on',...                    % Draw plot in figure.
    'LineWidth',3, ...                    % Plot line is thick.
    'Color',[0 0 0]);                     % Plot line is black.
axis([0 12 0 3000]);
axis 'auto y';
ylabel('Power (Relative units)');
xlabel('Harmonic Mode (m)');

hold on;
index = find(power == max(power));
mainPeriodStr = num2str(m(index));
plot(m(index),power(index),'r.', 'MarkerSize',25);
% text(m(index)+0.5,power(index)-0.05,['m = ',mainPeriodStr]);
text(m(index)-0.3,power(index)+0.035,['m = ',mainPeriodStr]);
hold off;
                
saveas(h3,FFTOutString,'fig')         % Save a MATLAB figure
saveas(h3,FFTOutString,'eps')         % Save as .eps 

close(h3)

end

