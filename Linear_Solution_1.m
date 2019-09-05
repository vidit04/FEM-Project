function [x,y ,u_array,x1,stress_rr,stress_array] = Linear_Solution_1(nelem,tau)
%% Generate list of position of nodes according to a geometric series
%    for assignement in "Nonlinear Finite Element Methods" 
%    in summer term 2019
%    lecturer in charge: Dr. Geralf Hütter
%
%Input Parameters
ro=40; %outer radius
ri=10; %radius of inclusion
%nelem=10; %number of elements
meshrefinementfactor=5; %ratio of element sizes at outer and inner radius
%
%ratio between element sizes of subsequent elements for a geometric series
q=meshrefinementfactor^(1./(nelem-1));
%size of first element
lelem=(ro-ri)*(1-q)/(1-meshrefinementfactor*q);
rnode=ri;
rnodes=[ri];
%loop over all elements
for i=1:nelem
        rnode=rnode+lelem;
        rnodes=[rnodes;rnode];
        lelem=lelem*q;
end
%visualize location of nodes
plot(rnodes,zeros(1,nelem+1),'x')
xlabel('r [\mum]')


E = 70000 ;
neu = 0.25;
sigma_y = 70;
ri = 10;
ro = 40;
epsilon_v = 0.01;

lambda = (E*neu)/((1-2*neu)*(1+neu));
%disp(lambda);
meu = E/(2*(1+ neu));
%disp(meu);
element_r = zeros(nelem,2);

for  i = 1:nelem
    element_r(i,1) = rnodes(i);
    element_r(i,2) = rnodes(i+1);
end
%disp(element_r);
delta_u = zeros(nelem,1);
Kt = zeros(nelem+1,nelem+1);
Gk = zeros(nelem+1,1);
epsilon_p = zeros(3,nelem);
stress_ = zeros(3,nelem);
%stress_rr = zeros(1,11);
%tau = 0.01;
u = zeros(nelem+1,1);
%u(1,1) = (1/3)*tau*epsilon_v*ri;
for d = 1:1/tau
    %disp(d*tau*epsilon_v);
    epsilon_init = (((1+neu)*sigma_y)/E);
    %disp(epsilon_init);
    if ((d*tau*epsilon_v) > epsilon_init)
        %disp(d);
        %disp('I am here.')
        break;
    end
    u(1,1) = (1/3)*tau*d*epsilon_v*ri;
    for k = 1:100
        %disp(u)
        %disp(rnodes(1,1))
        Kt = zeros(nelem+1,nelem+1);
        Gk = zeros(nelem+1,1);
        for e = 1:nelem
            Ae  = zeros(2,nelem+1);
            Ae(1,e) = 1;
            Ae(2,e+1) = 1;
            %disp(Ae);
            u_e = Ae*u;
            %disp(u_e);
            %disp(element_r(e,:));
            %disp(epsilon_p(:,e));
            [Kt_e, Fint_e,epsilon_p_return,stress] = elementrout(u_e,element_r(e,:),E,neu,sigma_y,lambda,meu,epsilon_p(:,e),e);
            Kt = Kt + transpose(Ae)*Kt_e*Ae;
            %disp(Kt);
            Gk = Gk + transpose(Ae)*Fint_e;
            %disp(Gk);
            epsilon_p(:,e) = epsilon_p_return;
            stress_(:,e) = stress; 

        end
        stress_rr = stress_(1,:);
        %if e==1
        Kt_red = Kt(2:nelem+1,2:nelem+1);
        %disp(Kt_red);
        Gk_red = Gk(2:nelem+1);
        %disp(Gk_red);
        %%% Kt *delta_u = - G 
        delta_u = inv(Kt_red) *( -Gk_red) ;
        %disp(delta_u);
        u(2:nelem+1) = u(2:nelem+1) + delta_u ;
        %disp(size(Kt_red));
        %disp(size(Gk_red));
        %disp(u);
        %disp(delta_u);
        %disp(epsilon_p)
        %disp(size(u_e))
        if (norm(delta_u,inf)<=0.005* norm(u,inf)) | (norm(Gk_red,inf)<=0.005*norm(Gk,inf))
            X = ['Element = ',num2str(nelem),' Step no ',num2str(d),' at Volumatic strain = ',num2str(d*tau*epsilon_v),' and Uo = ',num2str(u(1,1)),' Number of iteration required to converge is ',num2str(k-1),'.'];
            disp(X);
            break;

        end
    end
end
  
u_array = zeros(nelem+1,1);
stress_array = zeros(nelem,1);
for j = 1:nelem+1
    [u_exact] = exact_solution_u(rnodes(j),((d-1)*tau*epsilon_v), ri, E);
    u_array(j,1) = u_exact;
    %stress_array(j,1) = stress_rr_exact;
end
%disp(u_array);
%disp(stress_array);
x = rnodes;
y = u;
f1= figure;
plot(x,y,x,u_array,'b--o');
ylabel('u(r)');
xlabel('r in [\mum] (nodes)');
legend('u FEM elem = 10','u Exact');
title('Linear Elastic Convergence')

stress1_ = zeros(3,nelem);
for e = 1:nelem
    Ae  = zeros(2,nelem+1);
    Ae(1,e) = 1;
    Ae(2,e+1) = 1;
    u_e = Ae*u;
    [Kt_e_1, Fint_e_1,epsilon_p_return_1,stress_1] = elementrout(u_e,element_r(e,:),E,neu,sigma_y,lambda,meu,epsilon_p(:,e),e);
    stress1_(:,e) = stress_1;
end
stress_rr = stress1_(1,:);
%stress_pp = stress1_(2,:);

r_middle = zeros(nelem,1);
for k=1:nelem
    r_middle(k,1) = (rnodes(k+1)+ rnodes(k))/2;
    %disp(r_middle);
    [stress_rr_exact] = exact_solution_stress(r_middle(k,1),((d-1)*tau*epsilon_v), ri, neu, E);
    stress_array(k,1) = stress_rr_exact;
end
x1 = r_middle;
f2 = figure;
plot(x1,-stress_rr,x1,-stress_array,'b--o');
ylabel('Stress rr');
xlabel('r in [\mum] (element Gauss Point)');
legend('Stress rr FEM elem = 10','Stress rr Exact');
title('Linear Elastic Convergence')

end