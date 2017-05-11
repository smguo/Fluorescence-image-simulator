function [state,o, logs] = simul8trMovement(state,o,logs)
% Simulates particle diffusion and flow, given a 2D matrix, inputObjects
% (containing only 0, 1, 2, etc) and diffCoeff, state.u_convection(1),state.u_convection(2)
% Does not support multiple states -- call a new simul8trMovement separately,
% and pass it the inputObjects matrix of a different state

% August 30, 2004
% By DK

% Do different kinds of particle movement here
%
% All calculations are performed assuming units (length, time, energy) = (micrometers, seconds,
% k_BT)
%

% ------------------------
% set up some short cuts
% ------------------------

D = o.diff_coeff;
L = o.sim_box_size_um;
exr = o.exposure ;
v = o.u_convection ;

if exr>o.time_step % simulate finite exposure time effect
    dt = o.time_step ; %(added 09/27/2010)
    n_obs = ceil(exr/dt); % number of steps observed
    n_steps = ceil(o.sec_per_frame/dt);
    logs.x = zeros([o.num_particle, o.n_dims, o.num_frames]);
    logs.tra = zeros([o.num_particle, o.n_dims, n_steps, o.num_frames]);

    for t = 1:o.num_frames; 
        sigma = sqrt(2*D*dt);   % Gaussian step size
        % record the position after each step, like logs, for the exposure time
        % feature
        traject= zeros([size(state.x),n_steps]) ; 
        %% new diffusion simulation
        dx = bsxfun(@times, randn(size(traject)), sigma) ; % all the dx due to diffusion over history are generated
        dx = bsxfun(@plus, dx, v*dt) ; % flow
        switch o.bound_condi
        case 'periodic'
            x = cumsum(dx,3) ;
            traject = bsxfun(@plus, x, state.x) ;
            state.x = traject(:,:,end) ;
            state.tra = traject ; 
            % -----------------------------
            % Periodic boundary conditions
            % -----------------------------
            for dim = 1:o.n_dims 
            state.x(:,dim) = mod(state.x(:,dim), L(dim));
                for j = 1:n_steps
                    state.tra(:,dim,j) = mod(state.tra(:,dim,j), L(dim));
                end
            end
        case 'random'
            for j = 1:n_steps
                 state.x = state.x + dx(:,:,j);
                 for dim = 1:o.n_dims 
                    ind_out = find(state.x(:,dim)<0|| state.x(:,dim)>=L(dim));            
                    state.x(:,dim)= mod(state.x(:,dim), L(dim));
                    dim2 = mod(dim+1,o.n_dims) ;
                    state.x(ind_out,dim2) = rand(numel(ind_out),1)*L(dim2);
                 end
                 traject(:,:,j) = state.x;
            end                       
        end
        state.tra = traject ; 
        logs.x(:,:,t) = state.x ; % Store the coordinates in matrix logs.x (particle, dim, time)
        logs.tra(:,:,:,t) = state.tra ;        
    end
elseif exr == o.time_step % fast simulation where no intermidiate time step within a frame is simulated
    dt = o.sec_per_frame ; %(added 09/27/2010)
    n_steps = o.num_frames;
    sigma = sqrt(2*D*dt);   % Gaussian step size
    % record the position after each step, like logs, for the exposure time
    % feature
    
    
    %% new diffusion simulation    
    switch o.bound_condi
    case 'periodic'       
        traject= zeros([size(state.x),n_steps]) ; 
        dx = zeros([size(state.x),n_steps]) ; 
         if n_steps> 1000
        for j = 1: floor(n_steps/1e3)
            dx(:,:,(j-1)*1000+1:j*1000) = randn([size(state.x) 1000]) ;
        end

        if ~mod(n_steps,1000)==0
            dx(:,:,j*1000+1:end) = randn([size(state.x) mod(n_step,1000)]) ;
        end
        else
                dx = randn([size(state.x) n_steps]) ;
    end 
    dx = bsxfun(@times, dx, sigma) ; % all the dx due to diffusion over history are generated
    dx = bsxfun(@plus, dx, v*dt) ; % flow
        x = cumsum(dx,3) ;
        clear dx
        traject = bsxfun(@plus, x, state.x) ;
        state.x = traject(:,:,end) ;
        logs.x = traject ;
        % -----------------------------
        % Periodic boundary conditions
        % -----------------------------
        for dim = 1:o.n_dims         
            for j = 1:n_steps
                  logs.x(:,dim,j) = mod(logs.x(:,dim,j), L(dim));
            end
        end   
     case 'random'    
         traject= zeros([size(state.x),n_steps]) ; 
         dx = zeros([size(state.x),n_steps]) ; 
          if n_steps> 1000
        for j = 1: floor(n_steps/1e3)
            dx(:,:,(j-1)*1000+1:j*1000) = randn([size(state.x) 1000]) ;
        end

        if ~mod(n_steps,1000)==0
            dx(:,:,j*1000+1:end) = randn([size(state.x) mod(n_step,1000)]) ;
        end
        else
                dx = randn([size(state.x) n_steps]) ;
    end 
    dx = bsxfun(@times, dx, sigma) ; % all the dx due to diffusion over history are generated
    dx = bsxfun(@plus, dx, v*dt) ; % flow
        for j = 1:n_steps
             state.x = state.x + dx(:,:,j);
             for dim = 1:o.n_dims 
                ind_out = find(state.x(:,dim)<0 | state.x(:,dim)>=L(dim));
%                 marg = 0.001 ;
%                 bound = [marg marg L(1)-marg L(2)-marg] ;
%                 state.x(ind_out,:) = bound(ceil(4*rand(numel(ind_out), o.n_dims))) ;
                state.x(:,dim)= mod(state.x(:,dim), L(dim));
                dim2 = dim+1 ;
                if dim2 > o.n_dims
                   dim2 = 1 ;
                end
                state.x(ind_out,dim2) = rand(numel(ind_out),1)*L(dim2);
             end
             traject(:,:,j) = state.x;   
        end                    
        logs.x = traject ;
      case '9boxes'    
          traject= zeros([size(state.x),n_steps]) ; 
          dx = zeros([size(state.x),n_steps]) ; 
        if n_steps> 1000
            for j = 1: floor(n_steps/1e3)
                dx(:,:,(j-1)*1000+1:j*1000) = randn([size(state.x) 1000]) ;
            end

            if ~mod(n_steps,1000)==0
                dx(:,:,j*1000+1:end) = randn([size(state.x) mod(n_step,1000)]) ;
            end
            else
                dx = randn([size(state.x) n_steps]) ;
        end 
    dx = bsxfun(@times, dx, sigma) ; % all the dx due to diffusion over history are generated
    dx = bsxfun(@plus, dx, v*dt) ; % flow
        for j = 1:n_steps
             state.x = state.x + dx(:,:,j);             
             for dim = 1:o.n_dims 
                ind_out = find(state.x(:,dim)<0 | state.x(:,dim)>=L(dim));
                new_x = rand(10*numel(ind_out),1)*L(dim) ;
                new_x(new_x>=L(dim)/3 & new_x<2*L(dim)/3) = [];
                state.x(ind_out,dim) = new_x(1:numel(ind_out));
             end    
             
                traject(:,:,j) = state.x;   
        end                    
        logs.x = traject ;
       case 'flow'        
           v_dim = find(v);           
           if ~numel(v_dim)==1
               error('flow boundary condition can only apply to simulation with flow in one of dimensions')
           end 
           step_v = floor(L(v_dim)/v(v_dim)/dt) ;   % time steps needed for particles to flow through the box
           state_up = state_rand_nodes(o) ; % generate random particle positions
           state_down = state_rand_nodes(o) ; % generate random particle positions
           state_up.x(:,v_dim) = state_up.x(:,v_dim) - L(v_dim) ;
           state_down.x(:,v_dim) = state_down.x(:,v_dim) + L(v_dim) ;
           state.x = [state.x; state_up.x ; state_down.x] ;
           n_ptl = size(state.x,1) ;
           dx = zeros([size(state.x),n_steps]) ; 
           traject= zeros([size(state.x),n_steps]) ; 
            if n_steps> 1000
                for j = 1: floor(n_steps/1e3)
                    dx(:,:,(j-1)*1000+1:j*1000) = randn([size(state.x) 1000]) ;
                end

                if ~mod(n_steps,1000)==0
                    dx(:,:,j*1000+1:end) = randn([size(state.x) mod(n_step,1000)]) ;
                end
            else
                dx = randn([size(state.x) n_steps]) ;
            end 
            dx = bsxfun(@times, dx, sigma) ; % all the dx due to diffusion over history are generated
            dx = bsxfun(@plus, dx, v*dt) ; % flow
           for j = 1:n_steps                                                              
                state.x = state.x + dx(:,:,j);  
               if mod(j, step_v)==0
                   ind_out = find(state.x(:,v_dim)>=L(v_dim));
                   ptl_out = state.x(ind_out,:); 
                   state.x(ind_out,:) = []; 
                   state_up = state_rand_nodes(o) ; % generate random particle positions
                   state_up.x(:,v_dim) = state_up.x(:,v_dim) -L(v_dim) ;
                   state.x = [state.x; state_up.x] ;
                   state.x = [state.x; ptl_out] ;
                   state.x = state.x(1:n_ptl,:) ;                   
               end 
               for dim = 1:o.n_dims 
                   if dim == v_dim
                     ind_out = find(state.x(:,dim))<0 ;
                     state.x(ind_out,dim)= mod(state.x(ind_out,dim), L(dim));
                   else
                     state.x(:,dim)= mod(state.x(:,dim), L(dim));
                   end
               end                
             traject(:,:,j) = state.x;   
            end                    
        logs.x = traject ;                 
            
    end         
end
end

