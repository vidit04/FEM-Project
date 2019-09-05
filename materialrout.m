function [C_matrix, stress, epsilon_p_return ] = materialrout(lambda,meu, epsilon_e,epsilon_p,sigma_y,e)

stress_trail = 2* meu*(epsilon_e-epsilon_p)+ lambda*sum(epsilon_e-epsilon_p) ;
stress_hyd = (1/3)*sum(stress_trail);
stress_dev = stress_trail - stress_hyd;
stress_vm = ((3/2)* transpose(stress_dev) * stress_dev)^(1/2);

if (stress_vm -  sigma_y) < 0  %%% Elastic case
    
    c11 = lambda + 2*meu;
    c12 = lambda;
    c13 = lambda;
    c21 = lambda;
    c22 = lambda + 2*meu;
    c23 = lambda;
    c31 = lambda;
    c32 = lambda;
    c33 = lambda + 2*meu;
    C_matrix = [c11,c12,c13; c21,c22,c23; c31,c32,c33];
    stress = stress_trail;
    epsilon_p_return = epsilon_p;
end

if (stress_vm - sigma_y) >= 0  %%%% Plastic case
    %disp(e*10);
    delta_l = sigma_y/stress_vm;
    C_matrix = ones(3)*(lambda+(2*meu)/3)+delta_l*((((2/3)*meu)*[2,-1,-1;-1,2,-1;-1,-1,2])-(1.5*(2*meu)/(stress_vm^2))*stress_dev*transpose(stress_dev));
    epsilon_p_return = epsilon_p + (1-delta_l)* stress_dev/(2*meu);  %%%% Update plastic strain
    stress = stress_trail - (1-delta_l)* stress_dev;  %%% Update stress
end

end