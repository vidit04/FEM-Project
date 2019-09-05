%%%%%% FEM Assignment Main File  %%%%%%%%%%% Vidit Gupta- 61442

%%%%%%% Please choose
%%%%%%% Question number value
%%%%%%% and Number of element
clear all;
clc;
question =1; %%%% Enter the Choice
             %%%% Question 0 Elastic Covergence for different graph
             %%%% Question 1 Elastic Covergence for single grapgh 
             %%%% Question 2 Plastic Convergence
             %%%% Question 3 Final Load Distribution
             %%%% Question 4 Non Linear Plot  for Radial Stress
nelem = 10; %%%%%% Value can be changed for question 0,1,2
tau =0.01; %%%%%% Value to be inputed for Question 2 
          %%%%%%% (nelem = 10 , 20 , 30 tau = 0.1 , 0.01 , 0.001)
          %%%%%%%
          %%%%%%% 

%%%%%%%%%%%%%%%%  END of INPUT   %%%%%%%%%%%%%%%%%
if question ==0
    
    Linear_Solution_1(nelem,tau);
end
if question ==1

    [x,y ,u_array,x1,stress_rr,stress_array] = Linear_Solution(nelem,tau);
    [x_2,y_2 ,u_array_2,x1_2,stress_rr_2,stress_array_2] = Linear_Solution(nelem-2,tau);
    [x_3,y_3 ,u_array_3,x1_3,stress_rr_3,stress_array_3] = Linear_Solution(nelem-4,tau);
    
    f1= figure;
    plot(x,y,x_2,y_2,x_3,y_3,x,u_array,'b--o');
    ylabel('u(r)');
    xlabel('r in [\mum] (nodes)');
    legend('u FEM elem = 10','u FEM elem = 8','u FEM elem = 6','u Exact');
    title('Linear Elastic Convergence')
    %Linear_Solution(nelem-2);
    %Linear_Solution(nelem-4);
    f2 = figure;
    plot(x1,-stress_rr,x1_2,-stress_rr_2,x1_3,-stress_rr_3,x1,-stress_array,'b--o');
    ylabel('Stress rr');
    xlabel('r in [\mum] (element Gauss Point)');
    legend('Stress rr FEM elem = 10','Stress rr FEM elem = 8','Stress rr FEM elem = 6','Stress rr Exact');
    title('Linear Elastic Convergence')
end

if question ==2
    Plastic_solution(nelem,tau);
end

if question ==3
    final_load_distribution(nelem);
end

if question==4
    stress_history_rr();
end