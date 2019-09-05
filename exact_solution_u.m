function [u_exact] = exact_solution_u(r,epsilon_v, ri,E)
u_exact = ((ri^3)* epsilon_v )/ (3 * (r^2));
%stress_rr = (- 2*E*epsilon_v*(ri^3))/ (3*(1+neu)*(r^3));
end