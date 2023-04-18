function [F_Riemann] = Roe_2D_flux(W_i, W_j, nu_ij)

global g;

%left
V_L = W_to_V_2D(W_i);

rho_L = V_L(1);
u_L = V_L(2);
v_L = V_L(3);
p_L = V_L(4);
E_L = W_i(4);
H_L = E_L + p_L;

%right
V_R = W_to_V_2D(W_j);

rho_R = V_R(1);
u_R = V_R(2);
v_R = V_R(3);
p_R = V_R(4);
E_R = W_j(4);
H_R = E_R + p_R;

% Roe averaged
a = sqrt(rho_L);
b = sqrt(rho_R);

u_RL = (a*u_L + b*u_R)/(a+b);
v_RL = (a*v_L + b*v_R)/(a+b);
h_RL = (H_L/a + H_R/b)/(a+b);
rho_RL = a*b;
E_RL = h_RL*rho_RL/g + ((g-1)/(2*g))*rho_RL*(u_RL^2 + v_RL^2);

W_RL = [rho_RL; rho_RL*u_RL; rho_RL*v_RL; E_RL];

% Jacobians
[A_RL, B_RL] = jacobians_2D(W_RL);

% matrix
A_ij = nu_ij(1)*A_RL + nu_ij(2)*B_RL;

% matrix abs
[Q_inv, L] = eig(A_ij);
A_ij_abs = Q_inv*abs(L)/Q_inv;

% left and right fluxes
F_L = Flux_2D(W_i)*transpose(nu_ij);
F_R = Flux_2D(W_j)*transpose(nu_ij);

% Final flux
F_Riemann = 0.5*(F_L + F_R) - 0.5*A_ij_abs*(W_j - W_i);

end

