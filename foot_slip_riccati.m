%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% foot_slip_riccati.m
% By William Clark, 25 October 2024
%
% This script solves the steady-state periodic temproally triggered Riccati
%   equation. 
% The (stable) periodic orbit is found by propagating an initial guess
%   forward until convergence is achieved
% Once the periodic orbit is found, the state-transition matrix is
%   approximated via the Peano-Baker series
% Finally, the resulting Riccati is solved via "fsolve"
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Script parameters
% Model parameters
alpha = pi/16; delta = pi/8; epsilon = 2; mu = 2;
% Iteration tolerence
tol = 1e-6;
% Number of steps to approximate the state-transition matrix
N_STM = 1000;
% Weight matrices for HLQR
B = [0;1;0]; Q = eye(3); R = 1;

%% Finding the unforces periodic orbit
% The state is z = [theta;xi;eta]
% Initial guess for the initial conditions for the periodic orbit
z0 = [delta; -2; 2]; z_old = [0;0;0];

% Define the full, nonlinear dynamics
C = @(z) mu / (1+mu*sin(z(1))^2);
f = @(t,z) [epsilon*z(2);...
            C(z)*(sin(z(1))*(cos(alpha)+epsilon*z(2)^2*cos(z(1))) - cos(z(1))*(z(3)+z(2)*cos(z(1))));...
            sin(alpha)-z(3)-z(2)*cos(z(1))];
Delta = @(z) [z(1)+2*delta; cos(2*delta)*z(2); z(3) + cos(delta)*(1-cos(2*delta))*z(2)];
% Defining the event function
footEvt = @(t,z) footEvent(t,z,delta);
options = odeset('Events',footEvt);

% Attempt to find a stable periodic orbit
% Iterate until tolerence is achieved
for i = 1:50
    [~,Z,te] = ode45(f, [0,10], z0, options);
    zf = Z(end,:)';
    z_old = z0; z0 = Delta(zf);
    if norm(z_old - z0) < tol
        break;
    end
end
if i == 50
    warning('Failed to converge to a periodic orbit');
end

% With the initial condition found, let us extract the full periodic orbit
sol = ode45(f, [0,10], z0, options);
% The reset time
kappa = sol.xe;

%% Approximating the State Transition Matrix
% The linearized system (about the trajectory)
row1 = [0,epsilon,0];
row3 = @(z) [z(2)*sin(z(1)), -cos(z(1)), -1];
a21 = @(z) (mu*cos(z(1))^2*(- epsilon*cos(alpha)*z(2)^2 + cos(z(1))*z(2) + z(3)))/(mu*(cos(z(1))^2 - 1) - 1) - (mu*(cos(z(1))^2 - 1)*(- epsilon*cos(alpha)*z(2)^2 + 2*cos(z(1))*z(2) + z(3)))/(- mu*cos(z(1))^2 + mu + 1) - (2*mu^2*cos(z(1))^2*(cos(z(1))^2 - 1)*(- epsilon*cos(alpha)*z(2)^2 + cos(z(1))*z(2) + z(3)))/(mu*(cos(z(1))^2 - 1) - 1)^2;
a22 = @(z) -(mu*cos(z(1))*sin(z(1))*(cos(z(1)) - 2*epsilon*z(2)*cos(alpha)))/(mu*sin(z(1))^2 + 1);
a23 = @(z) -(mu*sin(2*z(1)))/(2*(mu*sin(z(1))^2 + 1));
row2 = @(z) [a21(z), a22(z), a23(z)];
A = @(t) [row1; row2(deval(sol,t)); row3(deval(sol,t))];

% Construct the overall matrices
% C is the reset matrix and Z is the (state,co-state) matrix
C = [1, 0, 0; 0, cos(2*delta), 0; 0, cos(delta)*(1-cos(2*delta)), 1];
Z = @(t) [A(t), -B*inv(R)*B'; -Q, -A(t)'];

% Integrate to get the transition matrix via Peano-Baker
times = linspace(0,kappa,N_STM); dt = times(2)-times(1);
Phi = eye(6);
for i = 2:length(times)
    Phi = Phi * expm(dt*Z(times(i)));
end

% Extract out the information
P11 = Phi(1:3,1:3); P12 = Phi(1:3,4:6);
P21 = Phi(4:6,1:3); P22 = Phi(4:6,4:6);

% Attempt to solve the Riccati equation
EQN_matrix = @(S) C'*S*C*(P11+P12*S) - P21 + P22*S;
EQN_vector = @(z) reshape(EQN_matrix(reshape(z,[3,3])), [], 1);

% Initial guess - the identity
s0 = reshape(eye(3), [], 1);

% Solve
s_star = fsolve(EQN_vector, s0);
S_star = reshape(s_star, [3,3]);

disp(S_star);

%% The event function
function [position, isterminal, direction] = footEvent(t,z,delta)
    position = z(1) + delta;
    isterminal = 1;
    direction = 0;
end