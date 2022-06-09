function output_image=fitsread(filename,subset,scaleQ)

% FITSREAD reads a FITS image into a Matlab variable.
%
%  data = fitsread(filename); 
%  data = fitsread(filename,range);
%  data = fitsread(filename,range,scaleQ);
%
%  The first form indicates that the full image should be read in. The
%  data is scaled by the BZERO and BSCALE keywords, which are assumed to
%  exist in the file.
%
%  The second form indicates that a subimage should be read in. The
%  subimage is specified by a string of the form 'xlow:xhigh,ylow:yhigh'.
%  The string 'full' can also be used if the full image is to be read in.
%  The data will be scaled using the BZERO and BSCALE keywords, which
%  are assumed to exist in the file.
%
%  The third form is like the second form, except a 1 or 0 is used to
%  indicate whether the BZERO and BSCALE keywords should be searched
%  for and used to scale the data (1=yes, 0=no). This should be
%  unnecessary unless you're reading in non-standard FITS files where
%  the BZERO and BSCALE keywords do not exist.
%
%  Notes:
%    1. If a filename with no extension is supplied, this routine looks
%    for a file with an extension of '.fits'. 
%
%    2. If a sub-image is read in only the memory required to hold the
%    sub-image is used by the function (ie. the whole image is not read
%    in, which would be wasteful and possibly very slow).
%
%    3. To compare the image with one viewed with a standard image viewer
%    such as xv or SAOimage be sure to display the image with the 
%    following settings: 
%
%              axis xy; 
%              axis image;
%
%    so that the origin is at the lower left and the correct aspect ratio
%    is shown.
%
%  Examples:
%    im=fitsread('n2715_J');                        % Assumes .fits extension
%    im=fitsread('n2715_lJ.fits','10:313,15:175');  % Reads a subimage 
%    im=fitsread('n2715_lJ.fits','full',0);         % Full image, no scaling 
%
%  Known deficiences:
%  1. Does not read anything more complicated than a simple 2D image
%     ie. no binary tables or data cubes.
%
%  2. Only does a very cursory inspection of possible architectures
%     in order to figure out if the data needs to be byte swapped.
%     This is because I only test for architectures I am familiar with.
%     (I believe this routine will work with PC-compatibles, Macs, Sun
%     workstations, and possibly DEC Alphas). 
%
% Version 2.5
% R. Abraham, Institute of Astronomy, Cambridge University


% Changes
% Jan. 2014 -- Bitswapped for MACI64 architecture (Doug Shields)
% Mar. 99 -- added the capability of extracting a subimage from the 
%            full image without having to read the whole image into 
%            memory.
% Dec. 98 -- removed option to specify number of 36 card header units.
%            An option to scale the data has been substituted instead.
% Oct. 98 -- added code to detect 80X86 and Alpha machines and to read
%            data into these with the big endian switch.


%Useful FITS Definitions:
%
%  "card" = 80 byte line of ASCII data in a file. This
%           contains keyword/value pairs and/or comments.
%
%  "Header area" = a set of 36 "cards". Note that
%          FITS files must have an integer number of header
%          areas (typically between 1 and 6, depending on how much
%          header information is stored in the file.)
%
%
%Reading FITS files:
%
%The first few cards of the header area must give information
%in a pre-defined order (eg. the data format, number of axes,
%size etc) but after that the header keywords can come in any
%order. The end of the last card giving information is flagged by
%the keyword END, after which blank lines are added to pad the
%last header area so it also contains 36 cards. After the last card
%the data begins, in a format specified by the BITPIX keyword.
%The dimensions of the data are specified by the NAXIS1 and
%NAXIS2 keywords.
%
%Reference:   NASA/Science Office of Standards and Technology
%           "Definition of the Flexible Image Transport System"
%                            NOST 100-1.0
%
%           This and other FITS documents are available on-line at:
%           http://www.gsfc.nasa.gov/astro/fits/basics_info.html



% First try to figure out if we need to swap bytes. This is
% imperfect as I don't know the endian-ness of each
% architecture, so I'm only looking for ones I know for
% sure.
friend = computer;
if strmatch(friend,'PCWIN')
   bswap = 'b';
elseif strmatch(friend,'LNX86')
   bswap = 'b';
elseif strmatch(friend,'ALPHA')
   bswap = 'b';
elseif strmatch(friend,'MACI64')
   bswap = 'b';
else
   bswap = 'l';
end


%Figure out how to scale the data by looking for the BSCALE and
%BZERO keywords. I think these can be anywhere in the header (ie.
%they don't have mandatory locations). We look for these using the
%"fitsheader" routine which scans through the whole header block
%looking for specific keywords (a slow process). We want to
%do this here before opening the file, because the fitsheader function 
%also opens the file internally and we want to avoid mixing calls to 
%"fitsheader" with explicit reads from the data file within this
%M-file, to avoid screwing up the  internal file indexing.
if nargin<2 
   subset = 'full';   
   scaleQ = 1;
elseif nargin<3
   scaleQ = 1;
end

if scaleQ==1
   bscale=fitsheader(filename,'BSCALE');
   bzero=fitsheader(filename,'BZERO');
   if isempty(bscale) | isempty(bzero)
      warning('BSCALE or BZERO keyword missing from the FITS file.')
      disp('Assuming BSCALE=1, BZERO=0');
      bscale=1.;
      bzero=0.;
   end
else
   bscale=1;
   bzero=0;
end


%Now we open the file and grab the required keywords from the first
%five cards, which MUST contain specific information according to
%the FITS standard. We take advantage of the mandatory locations 
%to try to speed things up by reading all the required keywords in
%a single pass.
fid=-1;
if ~isstr(filename)
   filename=setstr(filename);
end;
if (isempty(findstr(filename,'.'))==1)
   filename=[filename,'.fits'];
end
[file,message] = fopen(filename,'r',bswap);
if file == -1
   error(message);
end
[d,simple,d]=parse_card(setstr(fread(file,80,'uchar')'));
[d,bitpix,d]=parse_card(setstr(fread(file,80,'uchar')'));
[d,naxis,d]=parse_card(setstr(fread(file,80,'uchar')'));
[d,naxis1,d]=parse_card(setstr(fread(file,80,'uchar')'));
[d,naxis2,d]=parse_card(setstr(fread(file,80,'uchar')'));
n_card=5;

%Keep reading cards until one turns up with the keyword 'END'.
keyword='NIL';
while(~strcmp(deblank(upper(keyword)),'END'))
   n_card=n_card+1;
   card=setstr(fread(file,80,'uchar')');
   [keyword]=parse_card(card);
end;
%Go past the blank lines of padding before the start of the data
if (rem(n_card,36) ~= 0)
   n_blanks = 36 - rem(n_card,36);
   dummy=fread(file,n_blanks*80,'uchar');
end;

%work out the range of rows and columns to read
if strcmp('FULL',upper(deblank(subset)))
   lx=1;
   hx=naxis1;
   ly=1;
   hy=naxis2;
   nr=naxis2;
   nc=naxis1;
else
   [lx,hx,ly,hy] = lowhigh(subset);   %note order here
   nc = (hx-lx)+1;
   nr = (hy-ly)+1;
end


%
%   Geometry of image reading: Matlab's fread begins reading
%   at the lower left of the image and moves to the right, then
%   up. The x-axis is naxis1, and the y-axis is naxis2. 
%
%                ----------------------------------------------
%                |                                            |
%                |                                            |
%                |                                            |
%                |                                            |
%                |                                            |
%                |                                            |
%       naxis2   |                                            |
%                |                                            |
%                | NB. Matlab starts reading at X (lower left |
%                |     of frame) and moves to the right...    |
%                |                                            |
%                | X------->                                  |
%                ----------------------------------------------
%                                  naxis1
%


%work out how many pixels to skip to get to first byte of data
pixel_start_skip = (ly - 1)*naxis1 + (lx - 1);

%work out how many pixels to skip from the start of a row to the beginning
%of the next row
pixel_row_skip = naxis1 - nc; 

% move to first byte 
if bitpix==-64
   byte_start_skip = pixel_start_skip*8;
   byte_column_skip = pixel_row_skip*8;
   type = 'float32';   
   fseek(file,byte_start_skip,'cof');
elseif bitpix==-32
   byte_start_skip = pixel_start_skip*4;
   byte_column_skip = pixel_row_skip*4;   
   type = 'float';
   fseek(file,byte_start_skip,'cof');   
elseif bitpix==8
   byte_start_skip = pixel_start_skip*1;
   byte_column_skip = pixel_row_skip*1;  
   type = 'uint8';   
   fseek(file,byte_start_skip,'cof');
elseif bitpix==16
   byte_start_skip = pixel_start_skip*2;
   byte_column_skip = pixel_row_skip*2;
   type = 'short';   
   fseek(file,byte_start_skip,'cof');
elseif bitpix==32
   byte_start_skip = pixel_start_skip*4;
   byte_column_skip = pixel_row_skip*4;
   type = 'long';   
   fseek(file,byte_start_skip,'cof');
else
   error('data type specified by BITPIX keyword is not -64, -32, 8, 16, or 32');
end;



% now start reading in the data. We read a column at a time and then
% skip to the next row.
X=zeros(nr*nc,1);
startIndex = 1;
finIndex = nc;
for i=1:nr  
   if i>1
      startIndex=finIndex+1;
      finIndex=startIndex+nc-1;
   end
   X(startIndex:finIndex)=fread(file,nc,type);
   fseek(file,byte_column_skip,'cof');
end

%Scale the data
X = bscale*X + bzero;

%Clean up and output data
fclose(file);

%Reshape to a 2D image
output_image=reshape(X,nc,nr)';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lx,hx,ly,hy]=lowhigh(str)
% LOWHIGH extracts range elements from a range string 
% of the form 'XL:XH,YL:YH'. The elements are returned
% as numbers not as strings.

[dx,dy]=strtok(str,',');
[lx,hx]=strtok(dx,':');
[ly,hy]=strtok(dy,':');

% remove bogus characters
hx = strrep(hx,':','');
ly = strrep(ly,',','');
hy = strrep(hy,':','');

% return as a number instead of as a string
lx = str2double(lx);
hx = str2double(hx);
ly = str2double(ly);
hy = str2double(hy);
