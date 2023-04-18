function [outputArg1,outputArg2] = plot_grid(xy, edge, edge_area_normal, bc_area_normal, bc_type, plot_normals)

global nx;
global ny;
global fontsize;

% initialize data used in the plots
x = xy(:,1);
y = xy(:,2);

colors = [0, 0, 1;
          0, 1, 0;
          1, 0, 0;
          0, 0.4, 0.6;
          0.4, 0.6, 0;
          0.6, 0, 0.4;
          0, 0, 0];

% plot grid
set(gca,'FontSize', fontsize);
hold on;
scatter(x,y,'k','filled');
for k = 1:length(edge)
    edge_temp = edge(k,:);
    ind_1 = edge_temp(1);
    ind_2 = edge_temp(2);
    A = xy(ind_1,:);
    B = xy(ind_2,:);
    X = [A(1), B(1)];
    Y = [A(2), B(2)];
    plot(X, Y, 'k', 'LineWidth', 1);
end

% plot normals if desired
if plot_normals
    for j = 1:ny
        for i = 1:nx
            % get all the normals of the dual cell of the (i,j) node
            normals = normals_dual_cell(i,j,edge_area_normal, bc_area_normal, bc_type);
            
            % define the positions of the normals in the plot
            P = zeros(7,2);
            if i <= nx-1 && j <= ny-1 && i >= 2 && j >= 2
                P(1,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i,:));
                P(2,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i+1,:));
                P(3,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i+1,:));
                P(4,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i,:));
                P(5,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i-1,:));
                P(6,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i-1,:));
            elseif j == ny && ~ismember(i, [1,nx])
                P(3,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i+1,:));
                P(4,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i,:));
                P(5,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i-1,:));
                P(6,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i-1,:));
            elseif j == 1 && ~ismember(i, [1,nx])
                P(1,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i,:));
                P(2,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i+1,:));
                P(3,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i+1,:));
                P(6,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i-1,:));
            elseif i == 1 && ~ismember(j, [1,ny])
                P(1,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i,:));
                P(2,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i+1,:));
                P(3,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i+1,:));
                P(4,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i,:));
            elseif i == nx && ~ismember(j, [1,ny])
                P(1,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i,:));
                P(4,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i,:));
                P(5,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i-1,:));
                P(6,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i-1,:));
            elseif i == 1 && j == 1
                P(1,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i,:));
                P(2,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i+1,:));
                P(3,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i+1,:));
            elseif i == 1 && j == ny
                P(3,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i+1,:));
                P(4,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i,:));
            elseif i == nx && j == 1
                P(1,:) = 0.5*(xy((j-1)*nx+i,:) + xy(j*nx+i,:));
                P(6,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i-1,:));
            elseif i == nx && j == ny
                P(4,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i,:));
                P(5,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-2)*nx+i-1,:));
                P(6,:) = 0.5*(xy((j-1)*nx+i,:) + xy((j-1)*nx+i-1,:));
            end
            P(7,:) = xy((j-1)*nx+i,:);
    
            % plot the normals of the dual cell of the (i,j) node
            for n = 1:7
                if ~isequal(normals(n,:),[0,0])
                    quiver(P(n,1), P(n,2), normals(n,1), normals(n,2), 'color', colors(n,:));
                end
            end
            
        end
    end
end
hold off;
% title adapts if normals are desired or not
if plot_normals
    title("Grid  |  n_x = "+nx+"  |  n_y = "+ny+"  |  with normals");
else
    title("Grid  |  n_x = "+nx+"  |  n_y = "+ny+"  |  without normals");
end
grid on;
axis equal;

end

