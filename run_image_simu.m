% dir_name ='D:\Manuscript FCS application\simulated_images\' ;
dir_name ='C:\data\simulated_images\' ;
B_v = [5 10 20 30].*1e3 ;
for q = 1:1
% for q = 1:numel(B_v)
%%
o = struct( ...
    ...%%%%%%%% parameters for Brownian dynamics simulation %%%%%%%%%%%%%%%%
     'BD_type', 'free' ...  % BD simulation type. 'free', 'mdomain', or 'corral'
    , 'n_dims', 2 ...      %  number of dimensions in which the simulation is conducted  
    , 'sim_box_size_um',  NaN ... % the simulated box is of this size. See code below for actual default value (should be box_size_px*um_per_px)
    , 'num_frames',     5e4 ...			% number of frames in the stack produced by the code
    , 'density',        20 ...              % density of the particles, in um^-o.n_dims
    , 'um_per_px',      0.24 ...			% resolution, microns per pixel (Laura: 0.088)
    , 'sec_per_frame',  1e-3 ...		% seconds per frame (interval)
    , 'exposure',       1e-3 ...    % exposure time (sec)
    , 'time_step',      1e-3 ...   %  time step for BD simulation (sec)
    ...
    , 'diff_coeff',     3.6 ...              	% diffusion coefficient in micron^2/sec    
    , 'u_convection',   [0 0] ...        	% convection velocity in microns/sec
    , 'make_gradient' , 0 ...  
    , 'num_particle',   1e3 ...                 % number of particles
    , 'bound_condi',    'periodic' ...    % boundary condition, 'periodic' or 'random'
    ...
    ...%%%%%%%% parameters for image generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    , 'box_size_px',    [2 2 2].^4 ...  	% size of the images.
    , 'psf_type', 'g' ...            	% only gaussian psf ('g') is currently supported
    , 'psf_sigma_um',   [0.16 0.16 0.3] ...  	% standard deviation of the psf in microns, indepedent for all three directions
    , 'bead_diameter',  1 ...               %
    , 'renderer', 'points_from_nodes' ...  	%   'lines_from_bonds' or 'points_from_nodes'.
    ...   
    , 'brightness',        2e3 ...          %   count rate per particle at the center of detection volume
    , 'signal_background',  0 ...  % background light     
    , 'offset',      96 ...   % camera offset
    , 'readout_noise',         3.19  ...         % dark noise rms(std) 6
    , 'EMgain',            164 ...
    , 'ADCgain',           16.92 ... 
    ...
    , 'finer_grid' , 2 ... %%% to more accurately simulate bead position   
);
o.sim_box_size_um =  o.box_size_px* o.um_per_px*8 ;
% o.Pe = o.psf_sigma_um(1) * sqrt(o.u_convection*o.u_convection')/o.diff_coeff;
%%
% w_eff = 2*o.psf_sigma_um(1)*((2.311*2*o.psf_sigma_um(1)/o.um_per_px + 0.06765)^-2+1) ; % correction factor for finite pixel size. see Sisan BPJ 2006.
w_eff = 2*o.psf_sigma_um(1) ;

% o.charac_time_1 = w_eff.^2/4/o.diff_coeff_1;
% o.charac_time_2 = w_eff.^2/4/o.diff_coeff_2;

o % show the parameters

%% Add immobile beads
%   state1 = state_rand_nodes(o) ;
% %   o.diff_coeff = 25e-10 ;
% [logs1 state1 o] = BD_simul8tr(state1, o);
% 
% logs.bead = logs1.x ;
% [im o_im] = image_generator2(logs,o);
% figure(16)
% imsequence_play(im,.1);
%% two components
% o.num_particle_1 = floor(0.5* o.num_particle) ;
% o.num_particle_2 = floor(0.5* o.num_particle) ;
% o.num_particle = o.num_particle_1 + o.num_particle_2 ;
% % o.diff_coeff_1 = 10.^-0.5* o.diff_coeff ;
% % o.diff_coeff_2 = 10.^0.5* o.diff_coeff ;
% % o.diff_coeff = 10.^0.5* o.diff_coeff ;
% o.diff_coeff_1 = 10.^-1* o.diff_coeff ;
% o.diff_coeff_2 = 10.^1* o.diff_coeff ;
% o.diff_coeff = cat(1,repmat(o.diff_coeff_1, [o.num_particle_1 1]),repmat(o.diff_coeff_2, [o.num_particle_2 1])) ;


o.charac_time = w_eff^2./4./o.diff_coeff;
imf = zeros(o.box_size_px(1),o.box_size_px(2), o.num_frames, 'uint16') ;
%% mobile beads
state = state_rand_nodes(o) ;
o.num_frames = o.num_frames/10 ;
for j = 1:10
tic
[logs state o] = BD_simul8tr(state, o);
toc
%%
% [im o_im] = image_generator2(logs,o);

%% vary brightness
% br_v = [2:4] ;
% for j = 1:numel(br_v)
% o.brightness = 10.^br_v(j) ;
[im o_im] = image_generator4(logs,o);
im = uint16(im) ;
imf(:,:,(j-1)*o.num_frames+1:j*o.num_frames) = im ;
clear logs
end
o.num_frames = o.num_frames*10 ;
im = imf ;
%     if max(im(:))>2^16-1 ;
%         disp('Intensity values are over-scale. Images will be saved in the double-precistion format.')
%     else
%         im = uint16(im) ;
%     end
file_name = ['im_D=' num2str(o.diff_coeff) ...
    '_B=' num2str(o.brightness) '_N=' num2str(o.num_particle)...
    '_EM=' num2str(o.EMgain)...
    '_sb=' num2str(o.signal_background)...
    '_' o.bound_condi...
    '_box=' num2str(o.sim_box_size_um(1)/o.box_size_px(1)/o.um_per_px)...
    '_fm' num2str(size(im,3)/1000) 'k'] ;
save([dir_name file_name '.mat'],'o','im')  ;

% end
%%
figure(35)
imm = mean(im, 3) ;
imshow(imm, 'InitialMagnification', 'fit')
% imshow(im(:,:,1), 'InitialMagnification', 'fit')
    %axis([1 size(X,1) 1 size(X,2)]); 
    axis on
    colormap(gray)
    h = colorbar ;
    caxis auto
    title('mean intensity','FontSize',15)
    format_fig2(2)
%%
print(gcf ,'-dpng','-r300', [fig_path,num2str(fig_num, '%03d'),'.png'])
fig_num = fig_num +1 ;
% clear logs
%%
% figure(16)
% imsequence_play(im,.5);
% imsequence_play(im);
% hold on
%%
% TACF_curve_cov_calculation_simulation
end