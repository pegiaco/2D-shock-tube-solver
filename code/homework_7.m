clear all; close all; clc;

% Parameters

global nx;
global ny;
global x_min;
global x_max;
global y_min;
global y_max;
global g;
global W_L;
global W_R;
global delta_x;
global delta_y;
global fontsize;

g = 1.4;

x_min = -10;
x_max = 10;
y_min = -1;
y_max = 1;

p_L = 1e5; % N.m^{-2}
rho_L = 1; % kg.m^{-3}
u_L = 100; % m.s^{-1}

p_R = 1e4; % N.m^{-2}
rho_R = 0.125; % kg.m^{-3}
u_R = 50; % m.s^{-1}

W_L = [rho_L; rho_L*u_L; 0; p_L/(g-1) + 0.5*rho_L*(u_L^2)];
W_R = [rho_R; rho_R*u_R; 0; p_R/(g-1) + 0.5*rho_R*(u_R^2)];

delta_t = 1e-4;

fontsize = 14;

%% Grid generator with normals plotted (nx = 12, ny = 6)

nx = 12;
ny = 6;
delta_x = (x_max-x_min)/(nx-1);
delta_y = (y_max-y_min)/(ny-1);

XY = xy();

AREA = area(XY);

EDGE = edge();

EDGE_area_normal = edge_area_normal(XY);

BC_type = bc_type();

BC_area_normal = bc_area_normal(bc_type);

plot_normals = true;
figure(1);
plot_grid(XY, EDGE, EDGE_area_normal, BC_area_normal, BC_type, plot_normals);

%% Grid generator without normals plotted (nx = 100, ny = 20) (otherwise, long time to plot)

nx = 100;
ny = 20;
delta_x = (x_max-x_min)/(nx-1);
delta_y = (y_max-y_min)/(ny-1);

XY = xy();

AREA = area(XY);

EDGE = edge();

EDGE_area_normal = edge_area_normal(XY);

BC_type = bc_type();

BC_area_normal = bc_area_normal(bc_type);
plot_normals = false;
figure(2);
plot_grid(XY, EDGE, EDGE_area_normal, BC_area_normal, BC_type, plot_normals);


%% Probes time-histories (still nx = 100, ny = 20)

t_max = 0.01;
N_time = ceil(t_max/delta_t);

W = [];

% initialize states W
for i = 1:nx
    for j = 1:ny
        if XY((j-1)*nx+i,1) < 0
            W(:,(j-1)*nx+i) = W_L;
        else
            W(:,(j-1)*nx+i) = W_R;
        end
    end
end

W_current = W;

% initialize probes time-histories
A_ind = [nx/2, ny/2];
B_ind = [nx/2 + 1, ny/2];
C_ind = [nx/2 + 1, ny/2 + 1];
D_ind = [nx/2 + 1, ny];

V_A = W_to_V_2D(W(:,(A_ind(2)-1)*nx+A_ind(1)));
rho_A = [V_A(1)];
p_A = [V_A(4)];
V_B = W_to_V_2D(W(:,(B_ind(2)-1)*nx+B_ind(1)));
rho_B = [V_B(1)];
p_B = [V_B(2)];
V_C = W_to_V_2D(W(:,(C_ind(2)-1)*nx+C_ind(1)));
rho_C = [V_C(1)];
p_C = [V_C(2)];
V_D = W_to_V_2D(W(:,(D_ind(2)-1)*nx+D_ind(1)));
rho_D = [V_D(1)];
p_D = [V_D(2)];

% solve
for n = 1:N_time
    n
    for j = 1:ny
        for i = 1:nx
            F = zeros(4,1);
            
            % get the fluxes
            k = (j-1)*nx+i;
            for m = 1:length(EDGE)
                if k == EDGE(m,1)
                    W_i = W_current(:,k);
                    W_j = W_current(:,EDGE(m,2));
                    F(:,end+1) = Roe_2D_flux(W_i, W_j, EDGE_area_normal(m,:));
                elseif k == EDGE(m,2)
                    W_i = W_current(:,k);
                    W_j = W_current(:,EDGE(m,1));
                    F(:,end+1) = Roe_2D_flux(W_i, W_j, -EDGE_area_normal(m,:));
                end
            end

            if ismember(k, BC_type(:,1))
                index = find(k == BC_type(:,1));
                if BC_type(index,2) == 1
                    F(:,end+1) = wall_flux(W_current(:,k), BC_area_normal(index,:));
                elseif BC_type(index,2) == 2
                    F(:,end+1) = far_field_flux(W_current(:,k), BC_area_normal(index,:));
                end
            end
            
            % update W
            W(:,k) = real(W_current(:,k) - (delta_t/AREA(k))*sum(F,2));
        
            % update V
            V(:,k) = W_to_V_2D(W(:,k));

            % construct probes time-histories
            if isequal([i,j], A_ind)
                rho_A(end+1) = V(1,k);
                p_A(end+1) = V(4,k);
            elseif isequal([i,j], B_ind)
                rho_B(end+1) = V(1,k);
                p_B(end+1) = V(4,k);
            elseif isequal([i,j], C_ind)
                rho_C(end+1) = V(1,k);
                p_C(end+1) = V(4,k);
            elseif isequal([i,j], D_ind)
                rho_D(end+1) = V(1,k);
                p_D(end+1) = V(4,k);
            end

        end
    end

    W_current = W;
end

% plot probes time-histories

time = 0:delta_t:t_max;

figure(3);
plot_rho_p(time, t_max, rho_A, p_A, "A");

figure(4);
plot_rho_p(time, t_max, rho_B, p_B, "B");

figure(5);
plot_rho_p(time, t_max, rho_C, p_C, "C");

figure(6);
plot_rho_p(time, t_max, rho_D, p_D, "D");

%% Plot density, pressure, velocity at t = [0.002, 0.004, 0.006, 0.008, 0.01] s (still nx = 100, ny = 20)

% initialize the data used in the plots
names = {'Density', 'X velocity', 'Y velocity', 'Pressure'};
legends = ["\rho (kg.m^{-3})", "v_x (m.s^{-1})", "v_y (m.s^{-1})", "p (N.m^{-2})"];
X = x_min:delta_x:x_max;
Y = y_min:delta_y:y_max;

% 1st figure number (not to overlap with previous sections)
q = 7;


t_max = 0.01;
N_time = ceil(t_max/delta_t);

W = [];

% initialize states W
for i = 1:nx
    for j = 1:ny
        if XY((j-1)*nx+i,1) < 0
            W(:,(j-1)*nx+i) = W_L;
        else
            W(:,(j-1)*nx+i) = W_R;
        end
    end
end

W_current = W;

% solve
for n = 1:N_time
    n
    for j = 1:ny
        for i = 1:nx
            F = zeros(4,1);
            
            % get the fluxes
            k = (j-1)*nx+i;
            for m = 1:length(EDGE)
                if k == EDGE(m,1)
                    W_i = W_current(:,k);
                    W_j = W_current(:,EDGE(m,2));
                    F(:,end+1) = Roe_2D_flux(W_i, W_j, EDGE_area_normal(m,:));
                elseif k == EDGE(m,2)
                    W_i = W_current(:,k);
                    W_j = W_current(:,EDGE(m,1));
                    F(:,end+1) = Roe_2D_flux(W_i, W_j, -EDGE_area_normal(m,:));
                end
            end

            if ismember(k, BC_type(:,1))
                index = find(k == BC_type(:,1));
                if BC_type(index,2) == 1
                    F(:,end+1) = wall_flux(W_current(:,k), BC_area_normal(index,:));
                elseif BC_type(index,2) == 2
                    F(:,end+1) = far_field_flux(W_current(:,k), BC_area_normal(index,:));
                end
            end
            
            % update W
            W(:,k) = real(W_current(:,k) - (delta_t/AREA(k))*sum(F,2));
        
            % update V
            V(:,k) = W_to_V_2D(W(:,k));

        end
    end

    W_current = W;

    % plot
    if ismember(n, [round(N_time/5), round(2*N_time/5), round(3*N_time/5), round(4*N_time/5), N_time])
        figure(q);
        for i = 1:4
            subplot(2,2,i);
            surf(X, Y, reshape_array(real(V(i,:))));
            instant = n*delta_t;
            set(gca,'FontSize', fontsize);
            title(names(i)+"  |  t = "+instant+" s  |  n_{x} = "+nx+"  |  n_{y} = "+ny);
            c = colorbar;
            c.Label.String = legends(i);
        end
        q = q+1;
    end
end
