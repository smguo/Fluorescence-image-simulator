% Creat a 2D Gaussian filter where each pixel is the integral of a gaussian
% within the pixel area. This is more realistic than the MATLAB built-in
% gaussian filter
% 
function filter = gaussian_filter(N, s)
     siz   = (N-1)/2*[1 1];     
     [x,y] = meshgrid(-siz(2):siz(2),-siz(1):siz(1));
     x1 = (x + 0.5)./sqrt(2)./s ; % scale sigma to sqrt(2)
     y1 = (y + 0.5)./sqrt(2)./s ;
     x2 = (x - 0.5)./sqrt(2)./s ;
     y2 = (y - 0.5)./sqrt(2)./s ;
     filter = (erf2(x1)-erf2(x2)).*(erf2(y1)-erf2(y2)) ;                   
end
% error funciton with integration range [-Inf x]
function y = erf2(x)
   y = zeros(size(x)) ;
   y(x>=0) = (erf(x(x>=0))+1)/2 ;
   y(x<0) = (1-erf(-x(x<0)))/2 ;
end