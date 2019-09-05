function [Kt_e, Fint_e,epsilon_p_return,stress] = elementrout(u_e,element_r,E,neu,sigma_y,lambda,meu,epsilon_p,e)
%%%Z = 0
%%%N1 = 1/2*(1-Z)
%%%N2 = 1/2*(1+Z)

le = element_r(2) - element_r(1);
%disp(le);
%%% epsilon = B*u
%%% u = N*ui
%%% J = le/2
B = [-1/le, 1/le ; 1/(element_r(1) + element_r(2)) , 1/(element_r(1) + element_r(2)) ; 1/(element_r(1) + element_r(2)), 1/(element_r(1) + element_r(2))];
%disp(B);
epsilon_e = B*u_e;
%disp(epsilon_e);
[C_matrix, stress, epsilon_p_return ] = materialrout(lambda,meu, epsilon_e,epsilon_p,sigma_y,e);
%disp(C_matrix);
Kt_e = 2 * transpose(B)* C_matrix * B * (le/2) * (((element_r(1)+element_r(2))/2)^2);
%disp(Kt_e);
Fint_e = 2 * transpose(B) * stress * (le/2) * (((element_r(1)+element_r(2))/2)^2);
%disp(Fint_e);
%disp(Kt_e);
%disp(Fint_e);
end