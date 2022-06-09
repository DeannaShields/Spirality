function [PATH,FILE] = Extract_Filename(PATHFILE)

% This function inputs a path and file name PATHFILE, which must be a
% string.  The function separates PATHFILE into a pathname (PATH) and a 
% file name (FILE).  
%
% The function searches for all forward slashes ('/') and separates in
% input string at the rightmost forward slash.


% EXAMPLE
%
% >> PATHFILE = '/Users/test.jpg';
% >> [PATH,FILE] = Extract_Filename(PATHFILE)
% PATH =
% /Users/
% FILE =
% test.jpg


% COMPUTATION

% Find the rightmost forward slash

MARK = 0;                 % Location marker for the rightmost forward slash
N = numel(PATHFILE);      % Number of characters in the input
for j = 1:N
    if PATHFILE(j)=='/'   % Find the rightmost forward slash
        MARK = j;
    end
end

% What if there are no forward slashes?

if MARK==0
    PATH = '';
    FILE = PATHFILE;          
    return
end

% What if the final character is a forward slash?

if MARK==N
    PATH = PATHFILE;          
    FILE = '';
    return
end 

% Otherwise:

PATH = PATHFILE(1:MARK);
FILE = PATHFILE(MARK+1:N);

return




    

