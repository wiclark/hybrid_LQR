%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamic_programming.m
% By William Clark, 04 November 2024
%
% This script determines the optimal control for the bouncing ball via
% dynamic programming. In particular, this allows for Zeno solutions.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize
clc; clear; tic;

%% Initial parameters

% Cost weights
alpha = 10; beta = 10; kappa = 1;
% Set Costs
h = @(x)alpha*(x(1)-1)^2+beta*x(2)^2;
g = @(x,u,t)kappa*1/2*u^2;

% The time duration and discretization
T = 10;
times = linspace(0,T,150)';
dt = times(2)-times(1);
N = length(times);

% The controls discretization
U = linspace(-1,3,150);

% The state-space discretization
x1 = linspace(0,2,150);
x2 = linspace(-2,2,150);
m1 = (length(x1)-1)/(x1(end)-x1(1));
m2 = (length(x2)-1)/(x2(end)-x2(1));
[X1,X2] = meshgrid(x1,x2);

%% Discretization of the dynamics and setting trajectory information

% Dynamics and running cost
fd = @(x,u)next_step(x, dt, u);
gd = @(x,u,t)(dt*g(x,u,t));

% All the trajectory information
u_star = zeros(length(x1),length(x2),length(times));
J_star = u_star; % All the cost information

%% Beginning the loop by setting the terminal cost
% The terminal condition
J_old = zeros(length(x1),length(x2));
for i = 1:length(x1)
    for j = 1:length(x2)
        J_old(i,j) = h([x1(i);x2(j)]);
    end
end

%% The main loop
% Looping through all the times
figure;
for k = 1:length(times)-1
    index = N - k;
    
    % Loop through all the controls
    J_new = J_old;
    F_J = griddedInterpolant(J_old, 'linear', 'none');
    parfor i = 1:length(x1)
        dummy = zeros(length(x2),1);
        u_dum = dummy; reset_dum = dummy;
        for j = 1:length(x2)
            COST = U; % The cost to the next step
            RESET = U;
            for c = 1:length(U)
                % The cost is made up of two parts, the current
                c1 = gd([x1(i);x2(j)],U(c),times(index));
                % and the tail
                pt_new = fd([x1(i);x2(j)],U(c));
                % This point is not a point on our grid so we interpolate
                o1 = m1*(pt_new(1)-x1(1))+1;
                o2 = m2*(pt_new(2)-x2(1))+1;
                c2 = F_J(o1,o2);
                COST(c) = c1 + c2;
            end
            % Determine the best choice
            [M,I] = min(COST);
            dummy(j) = M;
            u_dum(j) = U(I);
        end
        % Record the information
        J_new(i,:) = dummy;
        u_star(i,:,index) = u_dum;
        J_star(i,:,index) = dummy;
    end
    J_old = J_new;
    % Display the output every 5 steps
    if mod(k,5)==0
        disp(k);
        surf(X1,X2,J_star(:,:,index)','edgecolor','none'); view([0,90]);
        pause(0.01);
    end
end

%% Display the value function
figure; L = linspace(0.8,1.599,8);
contourf(X1,X2,J_star(:,:,index)',L ,'ShowText',true,'FaceAlpha',0.5);
axis([0,1.25,-1.25,1.25]);
xlabel('$x_0$','Interpreter','Latex','FontSize',14);
ylabel('$v_0$','Interpreter','Latex','FontSize',14);
title('Value Function','Interpreter','Latex','FontSize',14);

%% Collecting and presenting the results
% Integrate the trajectories
u_fcn = griddedInterpolant(u_star);
u_f   = @(t,x)(u_fcn(m1*(x(1)-x1(1))+1,m2*(x(2)-x2(1))+1,t/dt));

% The flow
T = times; X = T; P = T; X(1) = 1; P(1) = 0;
for k = 2:length(T)
    x0 = X(k-1); p0 = P(k-1);
    x_out = next_step([x0;p0], dt, u_f(T(k-1),[x0;p0]));
    X(k) = x_out(1); P(k) = x_out(2);
end

% Plot
figure; hold on;
% Plot the optimal trajectory
subplot(211); hold on;
axis([0 times(length(times)) 0 1.25]); grid;
plot(T,X,'k','LineWidth',2);
xlabel('$t$','Interpreter','Latex','FontSize',14);
ylabel('$x$','Interpreter','Latex','FontSize',14);
% Plot the optimal control
subplot(212);
controls = T;
for j = 1:length(T)
    controls(j) = u_f(T(j),[X(j),P(j)]);
end
plot(T,controls,'k','LineWidth',2);
grid; axis([0,times(length(times)),-1,3]);
xlabel('$t$','Interpreter','Latex','FontSize',14);
ylabel('$u$','Interpreter','Latex','FontSize',14);

% Repeat with different initial condition
T = times; X = T; P = T; X(1) = 0.5; P(1) = 0;
for k = 2:length(T)
    x0 = X(k-1); p0 = P(k-1);
    x_out = next_step([x0;p0], dt, u_f(T(k-1),[x0;p0]));
    X(k) = x_out(1); P(k) = x_out(2);
end

% Plot
% Plot the optimal trajectory
subplot(211);
plot(T,X,'b--','LineWidth',2);
% Plot the optimal control
subplot(212); hold on;
controls = T;
for j = 1:length(T)
    controls(j) = u_f(T(j),[X(j),P(j)]);
end
plot(T,controls,'b--','LineWidth',2);

sgtitle('Optimal Trajectory and Control','Interpreter','Latex','FontSize',14);
%% End
toc;