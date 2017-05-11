function f=convolveGaussian(f, standard_deviation)
% function f=convolveGaussian(f, standard_deviation)
% Purpose: 
% convole an image f(y,x) with a Gaussian filter (PSF)
% with a *scalar* standard_deviation.  Periodic boundary conditions are
% used. 
% 
% Kirill Titievsky
% kir at mit dot edu
% 24-Aug-2009
%
%

% Single pixel psf has been useful for a number of tests and does not
% require convolution at all. This speeds up code execution considerably.
% So we treate the case when convolution is unnecessary at all separately.

if standard_deviation <= 0.1
    % do nothing
else
    % first, fspecial uses a two dimensional gaussian of exp(-x^2/(2*s^2)),
    % Therefore, I use
    s = max(eps, standard_deviation);   % modified by Syuan-Ming 2012
    % which produces a 2-D gaussian with the the standard deviation of the PSF
    % standard_deviation.
    % The miniminum used here is the smallest representable double.  This
    
    % Next, define the filter size.
    % Keep the filter size odd.  Otherwise, minimal std produces a 2x2 filter.
%     N = max(3, 6*ceil(standard_deviation) + 1);
    N = max(3, 2*ceil(4*standard_deviation) + 1);
%     filter=fspecial('gaussian', N, s);
    filter=gaussian_filter(N, s);
    % scale the filter to make the maximum 1
%     filter = filter/filter(ceil(N/2),ceil(N/2)) ;
    % convert the input image to double precision to avoid integer arithmetic
    if ~isfloat(f)
        f = double(f);
    end
    % this typically takes most of the time
    f=imfilter(f,filter,'circular','conv');
end