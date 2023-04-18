function [V] = W_to_V_2D(W)

global g;

rho = W(1);
u = W(2)/W(1);
v = W(3)/W(1);
E = W(4);
p = (g-1)*(E - 0.5*rho*(u^2 + v^2));

V = [rho; u; v; p];

end

