function [F_wall] = wall_flux(W, nu_ii)

V = W_to_V_2D(W);

p = V(4);

F_wall = [0; p*nu_ii(1); p*nu_ii(2); 0];

end

