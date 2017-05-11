function [state] = state_rand_nodes(o)
% generate initial state of randomly distributed nodes in the simulation box
% c is the number of nodes
% o is the structure of simulation parameters

c = o.num_particle ;
L = o.sim_box_size_um;
L = L(1:o.n_dims);
xc = [rand(c,1)*L(1) rand(c,1)*L(2)];
state.x = xc;
end