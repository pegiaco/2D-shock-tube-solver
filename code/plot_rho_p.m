function [] = plot_rho_p(time, t_max, rho, p, name)

global nx;
global ny;
global fontsize;

set(gca,'FontSize', fontsize);
subplot(2,1,1);
plot(time, real(rho), 'k', 'LineWidth', 1);
title("Density of probe "+name+"  |  t = "+t_max+" s  |  n_{x} = "+nx+"  |  n_{y} = "+ny);
legend("\rho (kg.m^{-3})");
subplot(2,1,2);
plot(time, real(p), 'r', 'LineWidth', 1);
title("Pressure of probe "+name+"  |  t = "+t_max+" s  |  n_{x} = "+nx+"  |  n_{y} = "+ny);
legend("p (N.m^{-2})");

end

