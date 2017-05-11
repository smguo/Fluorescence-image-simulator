function [logs state o] = BD_simul8tr(state, varargin)
% interpret the input
o_base = struct( ...
    ...%%%%%%%% parameters for Brownian dynamics simulation %%%%%%%%%%%%%%%%
    'BD_type', 'free' ...  % BD simulation type. 'free', 'mdomain', or 'corral'
    , 'n_dims', 2 ...      %  number of dimensions in which the simulation is conducted  
    , 'sim_box_size_um',  [2 2 2].^10 * 0.0645 ...  % the simulated box is of this size. See code below for actual default value (should be box_size_px*um_per_px)
    , 'num_frames', 10 ...			% number of frames in the stack produced by the code
    , 'density', .1 ...              % density of the particles, in um^-o.n_dims
    , 'sec_per_frame', 5 ...		% seconds per frame
    ...
    , 'diff_coeff', 0.01 ...              	% diffusion coefficient in micron^2/sec
    , 'u_convection', [0 0] ...        	% convection velocity in microns/sec
    , 'bonds_per_atom', 0 ...            	% # of bonds per atom
    , 'spring_const', 1 ...              	% spring constant 1/micron^2
    , 'topo', 'nn' ...               	% ('nn', 'random'). network topology nearest neighbor or random connections
    , 'make_gradient' , 0 ...  
    , 'ic', 'rand' ...               % initial conditions.  Leave at 'rand'.    ...
    ...%%%%%%%% parameters for image generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ...  
    , 'finer_grid' , 3 ...         %%% to more accurately simulate bead position
);


o = merge_ops(varargin, o_base);

% set up the log structure
logs.x = [];

%% simulate the system 
fprintf(1,'Running Brownian dynamics simulation....  ')
switch o.BD_type
    case 'free'
        [state,o, logs] = simul8trMovement(state,o, logs);
    case 'mdomain'
        [state,o, logs] = mdomain_simul8tr(state,o, logs);
end
        
fprintf(1,'Simulation done, generating images....')

end
