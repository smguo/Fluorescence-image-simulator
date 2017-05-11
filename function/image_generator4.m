function [imageFinal, o] = image_generator4 (logs,varargin)
% varargin=o
% interpret the input
o_base = struct( ...
    ...%%%%%%%% parameters for Brownian dynamics simulation %%%%%%%%%%%%%%%%
    'um_per_px', 0.0645 ...			    % resolution, microns per pixel
    ... 
    ...%%%%%%%% parameters for image generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    , 'box_size_px', [2 2 2].^8 ...  	% size of the images.
    , 'psf_type', 'g' ...            	% only gaussian psf ('g') is currently supported
    , 'psf_sigma_um', [0.4 0.4 0.7] ...  	% standard deviation of the psf in microns, indepedent for all three directions
    , 'renderer', 'points_from_nodes' ...  	%   'lines_from_bonds' or 'points_from_nodes'.
    ...    
    , 'brightness',        100 ...          %   observed intensity per node (added 09/27/2010)   
    ...
    , 'finer_grid' , 3 ... %%% to more accurately simulate bead position    
);


o = merge_ops(varargin, o_base);


imageFinalSize = [o.box_size_px(1)*2,o.box_size_px(2)*2]; % define output image size
imageFinal = zeros(imageFinalSize(1),imageFinalSize(2),o.num_frames); % define output image
% B_sf = mvnpdf([0 0],[0 0],(o.psf_sigma_um(1)/o.um_per_px)^2*[1 1]) ;

o.box_size_px = o.box_size_px*o.finer_grid*2; % expand image to have finer grid
o.um_per_px = o.um_per_px/o.finer_grid; % recale calibration factor

B = o.brightness*o.time_step ; % rescale the brightness so that it is always intensity/bead/frame time
if isfield(o, 'brightness1')
    B1 = o.brightness1*o.time_step ;
    B2 = o.brightness2*o.time_step ;
else
end

%% adding PSF and noise to create the images

% Preallocate intermediate images (resized, and need to be resized back)
image = zeros(o.box_size_px(1),o.box_size_px(2),o.num_frames);

dt = o.time_step ; %(added 09/27/2010)
exr = o.exposure ;
n_obs = ceil(exr/dt);% number of steps observed
n_ptl = size(logs.x,1); % number of particles in logs
% n_steps = ceil(o.sec_per_frame/dt);

%%%%%%%%%%%% add the particles to the images, and the same time add signal background to images %%%%%%%%%%%
% check if input contains bead position
if isfield(logs, 'bead')      
% convert unit of the radius of bead to #pixel
radius = 0.5*(o.bead_diameter/o.um_per_px) ;
B_scale = 0.5*o.bead_diameter/max([0.3*o.psf_sigma_um(1) o.um_per_px]) ;
% B_scale = 0.5*o.bead_diameter/o.um_per_px ;
    for t = 1:o.num_frames;    
        image(:,:,t) = image(:,:,t) + B.* 0.5*n_obs*B_scale^2 ... % scale the brightness of beads according area ratio
            * feval(o.renderer, logs.bead(:,:,t), o); % same as points_from_nodes(state, o)
        image(:,:,t) = convolveCircle(image(:,:,t), radius) ;
    end
else
end

if exr>dt
    for t = 1:o.num_frames;        
        avr = zeros(n_ptl*n_obs, o.n_dims);
        for j= 1:n_obs
        avr(n_ptl*(j-1)+(1:n_ptl),:) = logs.tra(:,:,j,t) ;
        end

        image(:,:,t) = image(:,:,t) + B.* feval(o.renderer, avr, o); % same as points_from_nodes(state, o)
    end
elseif exr == dt
     tic
     if isfield(o, 'brightness1')
         dom_pos_grid =  o.dom_pos_grid ;
         for t = 1:o.num_frames;  
             ind_dom = get_dom_ind2(logs.x(:,:,t), dom_pos_grid, o.dom_radius) ; % check states of particles (inside or outside domains)             
             ind_indom = ind_dom>0 ; % get indices
             ind_outdom = ind_dom == 0 ;
             image(:,:,t) = image(:,:,t) + B1.* feval(o.renderer, logs.x(ind_outdom,:,t), o); % same as points_from_nodes(state, o)
             image(:,:,t) = image(:,:,t) + B2.* feval(o.renderer, logs.x(ind_indom,:,t), o); % same as points_from_nodes(state, o)
         end
     else
         for t = 1:o.num_frames;     
           image(:,:,t) = image(:,:,t) + B.* feval(o.renderer, logs.x(:,:,t), o); % same as points_from_nodes(state, o)
         end
     end

     toc
end
%%%%%%%%%%%% pick a psf %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(o.psf_type,'g')
    psf = @(f)convolveGaussian(f,o.psf_sigma_um(1)/o.um_per_px);
elseif strcmp(PSFType,'l')
    psf = @(f) convolveGaussianLine(f,filter_size(1:2),o.psf_sigma_um(1)/o.um_per_px);
else
    psf = @(f) convolveAiry(f,filter_size(1:2),o.psf_sigma_um(1)/o.um_per_px);
end

%%%%%%%%%%%% Convolve new positions, resize back the image %%%%%%%%%%
tic
for t = 1:o.num_frames
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fprintf(1,['No.',num2str(t),' frame completed !\n']) % show progress       
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    image(:,:,t)  = psf(image(:,:,t)); % psf is a picked psf function


 end
 toc
 tic

imageFinal = bin_image_3(image,o.finer_grid) ;

 toc
 o.box_size_px = o.box_size_px/o.finer_grid/2; % rescale back image size
 imageFinal = imageFinal(0.5*o.box_size_px(1)+1:1.5*o.box_size_px(1),...
     0.5*o.box_size_px(2)+1:1.5*o.box_size_px(2), :) ;

o.um_per_px = o.um_per_px*o.finer_grid; % recale back calibrationfactor

imageFinal = imageFinal + o.signal_background;
imageFinal = poissrnd(imageFinal) ; % add photon counting noise 
imageFinal = gamrnd(imageFinal ,o.EMgain); % EM noise
imageFinal = imageFinal./o.ADCgain ;
imageFinal = imageFinal + o.offset + o.readout_noise*randn(size(imageFinal)); % add offset and readout noise

end

function ind_dom = get_dom_ind2(state_now, dom_pos_grid, dom_radius)
ind_dom = zeros(size(state_now,1), 1) ;

    for k = 1:size(dom_pos_grid,1)
         ind_dom_bi = (state_now(:,1)- dom_pos_grid(k,1)).^2+(state_now(:,2)- dom_pos_grid(k,2)).^2 < dom_radius.^2 ;         
%          ind_dom_s = find(ind_dom_bi>0, 1) ;         
%          ind_dom_s = find(ind_dom_bi>0, 1);       
            ind_dom(ind_dom_bi) = k ;                           
    end
end