function f=convolveCircle(f, diameter)
% Purpose: 
% convole an image f(y,x) containing beads position with a disk filter
% to generate finite size beads on the image 
% note: the output image needs to be convolve with the PSF again to get the
% real image
% with a *scalar* diameter.  Periodic boundary conditions are
% used. 
% 

% Single pixel psf has been useful for a number of tests and does not
% require convolution at all. This speeds up code execution considerably.
% So we treate the case when convolution is unnecessary at all separately.

if diameter <= 1
    % do nothing
else
    
    % Next, define the filter size.
    % the size of the disk filter is always odd
   
    filter=fspecial('disk', diameter);
    % convert the input image to double precision to avoid integer arithmetic
    if ~isfloat(f)
        f = double(f);
    end
    % this typically takes most of the time
    f=imfilter(f,filter,'circular','conv');
end