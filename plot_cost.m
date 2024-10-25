%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_cost.m
% By William Clark, 23 October 2024
%
% This script runs "find_trajectories.m" many times to draw the value
% function along with the initial control values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The domain
x = linspace(0.05,1,100);
y = linspace(0.05,1,100);
[X,Y] = meshgrid(x,y);
J = X; U = X;

% Loop over initial conditions
for i = 1:length(x)
    parfor j = 1:length(y)
        [J(j,i),U(j,i)] = find_trajectories([x(i);y(j)]);
        disp([i,j]);
    end
end

% Plot
figure; grid;
surf(X,Y,J,'EdgeColor','none');
xlabel('$x$','Interpreter','Latex','FontSize',14);
ylabel('$y$','Interpreter','Latex','FontSize',14);
zlabel('$J$','Interpreter','Latex','FontSize',14);
title('Value Function');

figure; grid;
surf(X,Y,U,'EdgeColor','none');
xlabel('$x$','Interpreter','Latex','FontSize',14);
ylabel('$y$','Interpreter','Latex','FontSize',14);
zlabel('$u$','Interpreter','Latex','FontSize',14);
title('Initial Control');