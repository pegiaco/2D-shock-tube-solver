function [A,B] = jacobians_2D(W)

global g;

V = W_to_V_2D(W);

rho = V(1);
u = V(2);
v = V(3);
p = V(4);
E = W(4);

A = [0, 1, 0, 0;
    -u^2 + 0.5*(g-1)*(u^2 + v^2), 2*u - (g-1)*u, -(g-1)*v, (g-1);
    -u*v, v, u, 0;
    -u*E/rho - (g-1)*u/rho*(E - 0.5*rho*(u^2 + v^2)) + (g-1)*u*0.5*(u^2 + v^2), E/rho + (g-1)/rho*(E - 0.5*rho*(u^2 + v^2)) - (g-1)*u^2, -(g-1)*u*v, g*u];

B = [0, 0, 1, 0;
    -u*v, v, u, 0;
    -v^2 + 0.5*(g-1)*(u^2 + v^2), -(g-1)*u, 2*v - (g-1)*v, (g-1);
    -v*E/rho - (g-1)*v/rho*(E - 0.5*rho*(u^2 + v^2)) + (g-1)*v*0.5*(u^2 + v^2), -(g-1)*u*v, E/rho + (g-1)/rho*(E - 0.5*rho*(u^2 + v^2)) - (g-1)*v^2, g*v];

end

