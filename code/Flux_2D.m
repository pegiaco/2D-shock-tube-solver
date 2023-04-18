function [F] = Flux_2D(W)

global g;

rho = W(1);
u = W(2)/W(1);
v = W(3)/W(1);
E = W(4);
p = (g-1)*(E - 0.5*rho*(u^2 + v^2));

F_x = [rho*u; rho*(u^2) + p; rho*u*v; (E + p)*u];
F_y = [rho*v; rho*u*v; rho*(v^2) + p; (E + p)*v];

F = [F_x, F_y];

end

