function [f o] = points_from_nodes_bug(state_cld, o) % statex is the particle positions with dimensions(n,3)
% A simulation cylinder is placed at the center of the image, assuming periodic boundary conditions.
% The image size is smaller than the length but larger than the width of the cyliner
% input is in cylindrical coordinate: state_cld = [z rho*theta]. rho is
% o.sim_box_size_um/2/pi
% x [0.5*Lbpx(2)-0.5*Lpx(2), 0.5*Lbpx(2)+0.5*Lpx(2)), 
% set up some convenient shortcuts
% 
% shortcut to box size dimension in pixels
% statex = state.x ;
if size(state_cld,1)== 0
    f = 0 ;
else    
Lpx = o.box_size_px;
Lb = o.sim_box_size_um ;
Lbpx = o.sim_box_size_um/o.um_per_px ;
statex = zeros(size(state_cld,1), 3) ;
statex(:,1) = state_cld(:,1) ;
% shortcut to box size in physical units
L = Lpx*o.um_per_px;
rho = Lb(2)/2/pi ;
[statex(:,2), statex(:,3)] = pol2cart(state_cld(:,2)/rho,rho) ;% convert from cylindrical to Cartesian
statex(:,3) = statex(:,3)+ rho ; % shift z from [-rho, rho] to [0, 2*rho]
% calculate the pixel positions corresponding to the particles 
iy = ceil(statex/o.um_per_px);
iy(:,2) = iy(:,2)+ 0.5*Lpx(2) ; % shift the center of the cylinder to the center of the image

% find the particles in the field of view
ivalid = ones(size(statex,1),1);
 for dim = 1:1  % only the dimension along the long axis of the cylinder needs to be cropped
     % the pixel coordinate must be inside the closed interval [1, Lpx]
    ivalid = ivalid & (iy(:,dim) <= 0.5*Lbpx(dim)+0.5*Lpx(dim));
    ivalid = ivalid & (iy(:,dim) > 0.5*Lbpx(dim)-0.5*Lpx(dim)); 
 end 

iy = iy(ivalid,:);
% shift the particle coordnates
for dim = 1:1
iy(:,dim) = iy(:,dim) - (0.5*Lbpx(dim)-0.5*Lpx(dim));
end
% add to the intensity of corresponding pixels
%   allocate image storage
f = zeros(Lpx(1), Lpx(2));
%   calculate linear index to pixels
linear_index = iy(:,1)  + (iy(:,2)-1)*size(f,1);

% Scale the intensity according to the z-position
%     w =  1.0*exp(-(statex(ivalid,3)-L(3)/2).^2/(2*o.psf_sigma_um(3).^2)) ; % Confocal
    w =  1.0*exp(-statex(ivalid,3)/o.psf_sigma_um(3)) ; % TIRF
    f(1:max(linear_index)) = accumarray(linear_index, w);
% end
end
end
