function [f o] = points_from_nodes(statex, o) % statex is the particle positions with dimensions(n,3)
% We assume that the image is cropped from the center of the 
% simulation box, assuming periodic boundary conditions.
% Specifically,  assume that the image is at [0.5*Lbpx(1)-0.5*Lpx(1), 0.5*Lbpx(1)+0.5*Lpx(1))
% x [0.5*Lbpx(2)-0.5*Lpx(2), 0.5*Lbpx(2)+0.5*Lpx(2)), while the
% simulation system is at [0, Lbpx(1)) x [0, Lbpx(2))
% set up some convenient shortcuts
% 
% shortcut to box size dimension in pixels
% statex = state.x ;
if size(statex,1)== 0
    f = 0 ;
else    
Lpx = o.box_size_px;
Lbpx = o.sim_box_size_um/o.um_per_px ;
% shortcut to box size in physical units
L = Lpx*o.um_per_px;

% calculate the pixel positions corresponding to the particles 
iy = ceil(statex/o.um_per_px);

% find the particles in the field of view
ivalid = ones(size(statex,1),1);
 for dim = 1:2
     % the pixel coordinate must be inside the closed interval [1, Lpx]
    ivalid = ivalid & (iy(:,dim) <= 0.5*Lbpx(dim)+0.5*Lpx(dim));
    ivalid = ivalid & (iy(:,dim) > 0.5*Lbpx(dim)-0.5*Lpx(dim)); 
 end 

iy = iy(ivalid,:);
% shift the particle coordnates
for dim = 1:2
iy(:,dim) = iy(:,dim) - (0.5*Lbpx(dim)-0.5*Lpx(dim));
end
% add to the intensity of corresponding pixels
%   allocate image storage
f = zeros(Lpx(2), Lpx(1));
%   calculate linear index to pixels
linear_index = iy(:,1)  + (iy(:,2)-1)*size(f,1);

%   add the contribution of each particle to the corresponding pixel
if o.n_dims == 2
    f(1:max(linear_index)) = accumarray(linear_index, 1.0);% reshape(histc(linear_index, 1:numel(f)), size(f));
elseif  o.n_dims ==3        
    % calculate the reduction in intensity due to out-of-plane z-position
    w =  1.0*exp(-(statex(ivalid,3)-L(3)/2).^2/(2*o.psf_sigma_um(3).^2)) ; 
    f(1:max(linear_index)) = accumarray(linear_index, w);
end
end
end
