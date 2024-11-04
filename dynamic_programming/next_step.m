%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% next_step.m
% By William Clark, 04 November 2024
%
% This function takes in the current state, time step, and control value
% and determines the output state (assuming constant control)
%  
% As the Zeno states/times for the bouncing ball are known analytically,
% the case of Zeno is explicitly handled.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function z_out = next_step(z, dt, u)

% Parameters
g = 1; c = 0.7; x0 = z(1); v0 = z(2);

% There are two main case: u<g and u>g
if u > g % <- In this case, there can be at most a single impact
    vnew = z(2) + (-g+u)*dt;
    xnew = z(1) + z(2)*dt + 1/2*(u-g)*dt^2;
    if xnew < 0 % <- So the single impact should occur
        t_hit = (-v0 - sqrt(v0^2-2*(u-g)*x0))/(u-g);
        if t_hit > dt
            error("We hit, but we don't?");
        end
        t_rem = dt - t_hit;
        vnew = -c^2*(v0 + (-g+u)*t_hit) + (-g+u)*t_rem;
        xnew = -c^2*(v0 + (-g+u)*t_hit)*t_rem + 1/2*(u-g)*t_rem^2;
    end
elseif u == g % <- The middle case as everything is linear
    vnew = z(2);
    xnew = z(1) + z(2)*dt;
    if xnew < 0 % <- So a single impact occurs
        t_hit = -x0/v0;
        vnew = -c^2*z(2);
        xnew = vnew*(dt-t_hit);
    end
else % u<g
    zeno_time = (v0 - (v0^2 + 2*(g-u)*x0)^(1/2))/(g-u) - (2*(v0^2 + 2*(g-u)*x0)^(1/2))/((g-u)*(c^2 - 1));
    if zeno_time <= dt % We go down the drain
        vnew = 0; xnew = 0;
    else
        % this is the garbage part
        % Time until first impact
        t_first = (-z(2)-sqrt(z(2)^2+2*z(1)*(g-u)))/(u-g);
        if dt < t_first
            vnew = z(2)+(u-g)*dt;
            xnew = z(1) + z(2)*dt + 1/2*dt^2*(u-g);
        else
            t_remaining = dt - t_first;
            vnow = -c^2*(z(2) + (u-g)*t_first);
            t_next = 2*vnow/(g-u); % Time until next impact
            while t_remaining > t_next
                vnow = c^2*vnow; t_remaining = t_remaining - t_next;
                t_next = 2*vnow/(g-u);
            end
            vnew = vnow + (u-g)*t_remaining;
            xnew = vnow*t_remaining + 1/2*t_remaining^2*(u-g);
        end
    end
end

z_out = [xnew;vnew];
end