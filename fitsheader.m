function out=fitsheader(filename,verbose)
%FITSHEADER returns header information from a FITS file.
%
%  FITSHEADER determines the value of a keyword or the number 
%  of 36 card header blocks in a FITS file. As a side effect 
%  the function can display all of the header information in the file 
%  onto the screen.
%
%  n_hdu=fitsheader(filename)
%  value=fitsheader(filname,'keyword')
%  n_hdu=fitsheader(filename,'silent_mode')
% 
%  The first form displays all header information on the
%  screen and returns the number of 36 card header blocks in the file.
%
%  The second form displays only the value corresponding to
%  a specified keyword.
% 
%  The third form does not display any header information, but
%  returns the number of 36 card header blocks in the file.
%
%  Useful FITS Definitions:
%
%  "card" = 80 byte line of ASCII data in a file. This 
%           contains keyword/value pairs and/or comments.
%           The FITS standard requires that cards occur
%           36 at a time. Blank lines are used to fill
%           in space if required.

% Version 1.0 
% R. Abraham, Institute of Astronomy, Cambridge University
%
% For FITS Info, check out the following document: 
% 
%      NASA/Science Office of Standards and Technology
%      Definition of the Flexible Image Transport System"
%      NOST 100-1.0
%
%      This and other FITS documents are available on-line at:
%      http://www.gsfc.nasa.gov/astro/fits/basics_info.html

if nargin<2
	verbose='FULL';
else
	verbose=upper(verbose);
end;

%Open the file
fid=-1;
if ~isstr(filename)
	filename=setstr(filename);
end;
if (isempty(findstr(filename,'.'))==1)
	filename=[filename,'.fits'];
end
[file,message] = fopen(filename,'r','l');
if file == -1
	error(message);
end


%Read the header information
n_card=0;
keyword=[];
while(~strcmp(upper(keyword),'END'))
    n_card=n_card+1;
    card=setstr(fread(file,80,'uchar')');
    [keyword,value,comment]=parse_card(card);
	if strcmp(verbose,'FULL')
		disp(card)
	elseif strcmp(keyword,verbose)
		out=value;
		return
	end
end;

%Clean up and output data 
fclose(file);

if (strcmp(verbose,'FULL') | strcmp(verbose,'SILENT_MODE'))
	out=ceil(n_card/36); %Return number of HDUs
else
	disp('Keyword not found.');
	out=[];
end


function [keyword,value,comment] = parse_card(s)
%Parses a FITS header card. 
%Reference:
%                NASA/Science Office of Standards and Technology
%           "Definition of the Flexible Image Transport System (FITS)"
%                        NOST 100-1.0    June 19, 1993
      
%Set defaults
keyword=[];value=[];comment=[];

%Get keyword in bytes 1 - 8
keyword=s(1:8);
if nargout==1
   return;
end

%If keyword is blank then the line is a comment
if strmatch(keyword,'       ')
    keyword=[];
	value=[];
	comment=deblank(s(11:80));
	return;
end;


%Keyword is non-blank. Check if there is a corresponding value/comment.
%If not then the only possibilities are that bytes 11:80 are a comment
%or that they are blank
if ~strmatch(s(9:10),'= ')
    keyword=deblank(keyword);
	value=[];
	comment=deblank(s(11:80));
	return;
end;

%Card is a standard keyword/value/comment structure. Break the value/comment
%string (bytes 11 - 80) into separate strings by tokenizing on "/" character.
%Remove the leading and trailing blanks on the value and the trailing blanks
%on the comment.

keyword=deblank(keyword);
[value,comment]=strtok(s(11:80),'/');
comment=deblank(comment);
value=fliplr(deblank(fliplr(deblank(value))));

%Now figure out whether to output the value as a string or as a number.
%The FITS standard requires all values that are strings to be in single
%quotes like this: 'foo bar', so I can simply look for occurences of a
%single quote to flag a string. However, logical variables can take the
%values T or F without having any single quotes, so I'll have to look
%out for those also.

%Test for logical. Return logical as a string.
if strmatch(upper(value),'T') | strmatch(upper(value),'F')
	return;
end;

%Test for string. Return string unconverted.
if length(findstr('''',value)) ~= 0
	return;
end;

%Only thing left is a number. Convert string to number.
value=str2num(value);
