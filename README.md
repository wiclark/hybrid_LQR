# hybrid_LQR
Code for running both the temporally and spatially triggered linear quadratic regulator.

<ol>
  <li>The file "foot_slip_riccati.m" numerically solves the steady-state periodic Riccati equation for the example of legged locomotion with a sliding foot.</li>
  <li>The file "find_trajectories.m" finds optimal arcs to a two-dimensional spatially triggered linear hybrid system.</li>
  <li>The file "plot_cost.m" runs "finds_trajectories.m" over numerous initial conditions to compute the value function.</li>
  <li>The file "double_spring.m" finds arcs to the double spring system.</li>
</ol> 
