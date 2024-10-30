%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% find_trajectories.m
% By William Clark, 23 October 2024
%
% This function takes in an initial state condition and determines the
% optimal trajectory under the prescribed linear hybrid dynamics.
% The optimal trajectories are found by determining the initial co-states
% that minimize the resulting cost. The minimization algorithm utilized is
% 'fminsearch.'
%
% Inputs: z0 = [x0;y0] <- The initial conditions as a column vector
%         toggle = 0   <- Suppresses the plot of the trajectory (default)
%                = 1   <- Plots the trajectory
%
% Outputs: J  <- The cost
%          u0 <- The initial control value
%          data <- a struct with fields:
%             first: a 1 by 4 vector of either terminal conditions 
%                       (no reset) or of the reset location
%           p_minus: the co-states immedeately before a reset
%                H1: The Hamiltonian immedeately before a reset
%            second: the terminal conditions (if reset)
%            p_plus: The co-states immedeately after a reset
%                H2: The Hamiltonian immedeately before a reset
%                p0: The initial co-states
%                H0: The initial Hamiltonian
%              jump: The difference, p^- - C'*p^+
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [J, u0, data] = find_trajectories(z0, toggle)

if nargin < 2
    toggle = 0;
end

% Set the system parameters
param.A = [0,1;-1,0];
param.B = [1;0];
param.C = [0,0;2,0];
param.L = [0;1];
param.Q = eye(2);
param.F = eye(2);
param.R = 1;
param.tf = 1;

% Make the augmented co-state dynamics
Rtilde = param.B*inv(param.R)*param.B'; param.Rtilde = Rtilde;
Z = [param.A,-Rtilde;-param.Q,-param.A'];

% The controls and running cost
u = @(z) -inv(param.R)*param.B'*z(3:4);
dJ = @(z) 1/2*dot(z(1:2),param.Q*z(1:2)) + 1/2*dot(u(z),param.R*u(z));

% The dynamics (with running cost)
f = @(t,z) [Z*z(1:4);dJ(z)];

% Determine the optimal initial co-states via fminsearch
find_p = @(p) shoot_forward(z0,f,param,p,0);
P = fminsearch(find_p, rand(2,1)-1); % <- Random initial guess

% Extract out the solutions for cost and control
u0 = -inv(param.R)*param.B'*P;
[J,D] = shoot_forward(z0,f,param,P,toggle);

% collecting the interesting data
data = D;
data.p0 = P';
data.H0 = 1/2*dot(z0,param.Q*z0) + dot(P,param.A*z0) - 1/2*dot(P,Rtilde*P);

end

%% Given a random initial co-state, we determine the resulting trajectory
function [Cost,D] = shoot_forward(IC,odefun,param,p0,toggle)
    
    % Declare the event function
    eventfcn = @(t,z) EventsFcn(t,z,param);
    options = odeset('Events',eventfcn);

    % Run the trajectory forward, halting if a reset occurs
    [~,Z,te] = ode45(odefun, [0,param.tf], [IC;p0;0], options);
    % Terminal cost and states
    Cost = Z(end,5); zf = Z(end,1:4); hit = zf(1:2);

    D.first = zf; D.p_minus = zf(3:4);
    D.H1 = 1/2*dot(hit,param.Q*hit') + dot(zf(3:4),param.A*hit') -...
        1/2*dot(zf(3:4),param.Rtilde*zf(3:4)');

    % Does an impact occur?
    if te < param.tf
        param.tf = param.tf - te; % <- Time to go
        ic = param.C*zf(1:2)'; % Post reset state
        % Call the function again and minimize over new co-states
        find_p = @(p) shoot_forward(ic,odefun,param,p,0);
        P = fminsearch(find_p,rand(2,1));
        % Integrate the trajectory
        [~,Z1] = ode45(odefun, [0,param.tf], [ic;P;Cost]);
        % Terminal cost and states
        Cost = Z1(end,5); zf = Z1(end,1:4); Z = [Z;nan(1,5);Z1];
        D.second = zf;
        D.p_plus = P';
        D.H2 = 1/2*dot(zf(1:2),param.Q*zf(1:2)') +...
            dot(zf(3:4),param.A*zf(1:2)') - ...
            1/2*dot(zf(3:4),param.Rtilde*zf(3:4)');
        D.jump = D.p_minus - D.p_plus*param.C;
    end
    
    % Terminal cost
    Cost = Cost + 1/2*dot(Z(end,1:2),param.F*Z(end,1:2)');

    if toggle == 1
        figure; hold on; grid;
        plot([0,1.5],[0,0],'k','LineWidth',2);
        plot([0,0],[0,1.5],'k','LineWidth',2);
        plot(Z(:,1),Z(:,2),'b','LineWidth',2);
        plot(hit(1),hit(2),'k.','MarkerSize',15);
        plot(IC(1),IC(2),'g.','MarkerSize',15);
        plot(Z(end,1),Z(end,2),'r.','MarkerSize',15);
        xlabel('$x$','Interpreter','Latex','FontSize',14);
        ylabel('$y$','Interpreter','Latex','FontSize',14);
        axis([-0.1,1.25,-0.1,1.25]);
    end
end


%% Event detection
function [position,isterminal,direction] = EventsFcn(t,z,param)
  L = param.L;
  position = dot(L,z(1:2)); % The value that we want to be zero
  isterminal = 1;  % Halt integration 
  direction = 0;   % The zero can be approached from either direction
end