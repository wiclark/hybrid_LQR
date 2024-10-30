%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% double_spring.m
% By William Clark, 30 October 2024
%
% This script solves numerically finds optimal trajectories for the double
% mass spring system in Section VI-C.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Script parameters
k1 = 2; k2 = 1;
d1 = 1; d2 = 1;
delta = 2; tf = 5;

% Building the matrices
A = [0,1,0,0;...
    -k1,0,0,0;...
    0,0,0,1;...
    0,0,-k2,0];
b = [0;k1*d1;0;k2*d2];
C = [1,0,0,0;...
    0, 0, 0, -1;...
    0,0,1,0;...
    0,-1,0,0];
lambda = [1;0;1;0];
Q = zeros(4,4); F = 50*eye(4); R = 1;
B = [0;1;0;0]; y = [0;0;delta;0];

F = diag([1,1,20,1]);

% Set parameters
param.A = A; param.b = b; param.C = C; param.lambda = lambda;
param.Q = Q; param.F = F; param.B = B; param.y = y; param.tf = tf;

%% The intitial condition
x0 = [d1;0;d2;0]; param.x0 = x0;

%% Set the dynamics
% The augmented dynamics
Rtilde = B*inv(R)*B'; param.Rtilde = Rtilde;
Z = [A,-Rtilde;-Q,-A'];
% The controls and running cost
u = @(z) -inv(R)*B'*z(5:8);
dJ = @(z) 1/2*dot(z(1:4),Q*z(1:4)) + 1/2*dot(u(z),R*u(z));
% The dynamics (with running cost)
f = @(t,z) [Z*z(1:8) + [b;zeros(4,1)];dJ(z)];

%% Determine optimal trajectories
find_p = @(p) shoot_forward(x0, f, param, p, 0);
optoptions = optimset('TolFun',1e-6,'TolX',1e-6,'MaxFunEvals',2000);
P = fminsearch(find_p, [-1;4.75;0;0], optoptions); % <- Random intial guess

% Extract out the solutions
u0 = -inv(R)*B'*P;
J = shoot_forward(x0, f, param, P, 1);

%% Determine the landscape
p1 = linspace(-1.5,2.5,200); p2 = linspace(3,7,100);
[P1,P2] = meshgrid(p1,p2); C = P1; Num = C;
for i = 1:length(p1)
    for j = 1:length(p2)
        [C(j,i),Num(j,i)] = shoot_forward(x0, f, param, [p1(i);p2(j);0;0], 0);
        disp([i,j]);
    end
end

figure; h = gca;
surf(P1,P2,C,'EdgeColor','none');
set(h,'zscale','log');
c = colorbar; set(gca,'ColorScale','log');
xlabel('$p_{x_1}(0)$','Interpreter','Latex','FontSize',14);
ylabel('$p_{v_1}(0)$','Interpreter','Latex','FontSize',14);
c.Label.String = 'Cost';
c.Label.FontSize = 14;
view([0,90]);

figure;
surf(P1,P2,Num,'EdgeColor','none'); c = colorbar;
c.Label.String = '# of Jumps';
c.Label.FontSize = 14;
view([0,90]);

%% Given an initial co-state, determine the resulting trajectory
function [Cost,N] = shoot_forward(IC, odefun, param, p0, toggle)
    % Declare the event function
    eventfcn = @(t,z) EventsFcn(t,z,param);
    options = odeset('Events',eventfcn);

    % Aux function
    beta = @(z) dot(param.lambda, param.A*z(1:4)) + dot(param.lambda,param.b);
    H = @(z) 1/2*dot(z(1:4),param.Q*z(1:4)) + dot(z(5:8),param.A*z(1:4)) - ...
        1/2*dot(z(5:8),param.Rtilde*z(5:8)) + dot(z(5:8),param.b);
    C = param.C;
    Cp = [eye(4),zeros(4);zeros(4),C']; Cx = [C,zeros(4);zeros(4),eye(4)];
    gamma = @(z) H(Cp*z) - H(Cx*z);
    epsilon = @(z) -gamma(z)/beta(z);

    % Flow forwards
    [T,Z,te] = ode45(odefun, [0, param.tf], [IC;p0;0], options);
    Cost = Z(end,9); zf = Z(end,1:8);

    % Number of impacts
    N = 0;

    % Does an impact occur?
    while te < param.tf
        p0 = inv(C)'*zf(5:8)' + epsilon(zf')*inv(C)'*param.lambda;
        ic = param.C*zf(1:4)';
        [T1,Z1,te] = ode45(odefun, [te, param.tf], [ic; p0; Cost], options);
        Cost = Z1(end,9); zf = Z1(end,1:8);
        Z = [Z;nan(1,9);Z1]; T = [T;nan;T1];
        N = N + 1;
    end

    % Terminal Cost
    Cost = Cost + 1/2*dot((Z(end,1:4)'-param.y),param.F*(Z(end,1:4)'-param.y));
    
    if toggle == 1
        figure; hold on; grid;
        plot(T,-Z(:,1),'b','LineWidth',2);
        plot(T,Z(:,3),'r','LineWidth',2);
        xlabel('$t$','Interpreter','Latex','FontSize',14);
        ylabel('Displacement','Interpreter','Latex','FontSize',14);
        legend({'$-x_1(t)$','$x_2(t)$'},'Interpreter','Latex','FontSize',14, 'Location','northwest');
    end

end

%% Event detection
function [position, isterminal, direction] = EventsFcn(t,z,param)
    L = param.lambda;
    position = dot(L,z(1:4));
    isterminal = 1;
    direction = -1; % Impact only occurs if the springs are moving together
end