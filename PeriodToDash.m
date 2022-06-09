function DASHED = PeriodToDash(STRINGNAME)

% This function inputs the string STRINGNAME, converts all the periods to
% dashes, and outputs the result as DASHED.  The function was created to
% fix filenames, since some systems have difficulty with filenames
% containing periods.



% EXAMPLE
%
% >> STRINGNAME = 'abc.de.f';
% >> DASHED = PeriodToDash(STRINGNAME)
%
% DASHED =
%
% abc-de-f



% COMPUTATION

N = numel(STRINGNAME);           % Number of characters in the input
for j = 1:N
    if STRINGNAME(j)=='.'        % Find the periods
        STRINGNAME(j)='-';       % Change the periods to dashes
    end
end

DASHED = STRINGNAME;

return


