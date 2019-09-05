function [stress_rr_exact] = exact_solution_stress(r,epsilon_v, ri, neu, E)
%u_exact = ((ri^3)* epsilon_v )/ (3 * (r^2));
stress_rr_exact = (- 2*E*epsilon_v*(ri^3))/ (3*(1+neu)*(r^3));
end