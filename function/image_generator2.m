function [imageFinal, o] = image_generator2 (logs,varargin)
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
    , 'signal_level', 100 ...               %signal level above background (disabled 09/27/2010)
    , 'brightness',        100 ...          %   observed intensity per node (added 09/27/2010)
    , 'signal_background', 200 ...
    , 'counting_noise_factor', 1.327  ...   % counting noise  factor (noise = sqrt(o.counting_noise_factor*imageFinal).*randn(im_dims) 
    , 'dk_background',  189.462 ...   % dark background average
    , 'dk_noise', 7.265 ...           % dark noise rms(std)
    ...
    , 'finer_grid' , 3 ... %%% to more accurately simulate bead position
    , 'store_x', 1 ...
);


o = merge_ops(varargin, o_base);


imageFinalSize = [o.box_size_px(1),o.box_size_px(2)]; % define output image size
imageFinal = zeros(imageFinalSize(1),imageFinalSize(2),o.num_frames); % define output image
o.box_size_px = o.box_size_px*o.finer_grid; % expand image to have finer grid
o.um_per_px = o.um_per_px/o.finer_grid; % recale calibration factor


%% adding PSF and noise to create the images

% Preallocate intermediate images (resized, and need to be resized back)
image = zeros(o.box_size_px(2),o.box_size_px(1),o.num_frames);


%%%%%%%%%%%% add the particles to the images, and the same time add signal background to images %%%%%%%%%%%
% check if input contains bead position
if isfield(logs, 'bead')      
% convert unit of the radius of bead to #pixel
radius = 0.5*(o.bead_diameter/o.um_per_px) ;
B_scale = 0.5*o.bead_diameter/max([0.3*o.psf_sigma_um(1) o.um_per_px]) ;
% B_scale = 0.5*o.bead_diameter/o.um_per_px ;
    for t = 1:o.num_frames;    
        image(:,:,t) = image(:,:,t) + o.brightness.* 0.5*n_obs*B_scale^2 ... % scale the brightness of beads according area ratio
            * feval(o.renderer, logs.bead(:,:,t), o); % same as points_from_nodes(state, o)
        image(:,:,t) = convolveCircle(image(:,:,t), radius) ;
    end
else
end

     tic
     parfor t = 1:o.num_frames;     
       image(:,:,t) = image(:,:,t) + o.brightness.* feval(o.renderer, logs.x(:,:,t), o); % same as points_from_nodes(state, o)
     end 
     toc

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
parfor t = 1:o.num_frames
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% fprintf(1,['No.',num2str(t),' frame completed !\n']) % show progress       
% %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    image(:,:,t)  = psf(image(:,:,t)); % psf is a picked psf function


 end
 toc
 tic

imageFinal = bin_image_3(image,o.finer_grid) ;

 toc
 o.box_size_px = o.box_size_px/o.finer_grid; % rescale back image size
 

o.um_per_px = o.um_per_px*o.finer_grid; % recale back calibrationfactor

imageFinal = imageFinal + o.signal_background;

imageFinal = imageFinal + sqrt(o.counting_noise_factor*imageFinal).*randn(size(imageFinal)); % add photon counting noise (signal noise). 

imageFinal = imageFinal + o.dk_background + abs( o.dk_noise*randn(size(imageFinal))); % add dark backround and noise

end