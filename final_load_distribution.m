function [] = final_load_distribution(nelem)
%% Generate list of position of nodes according to a geometric series
%    for assignement in "Nonlinear Finite Element Methods" 
%    in summer term 2019
%    lecturer in charge: Dr. Geralf Hütter
%
%Input Parameters
ro=40; %outer radius
ri=10; %radius of inclusion
%nelem= 10; %number of elements
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


E = 70000;
neu = 0.25;
sigma_y = 70;
ri = 10;
ro = 40;
epsilon_v = 0.01;

lambda = E*neu/((1-2*neu)*(1+neu));
meu = E/(2*(1+ neu));
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

tau = 1;
u = zeros(nelem+1,1);
u(1,1) = (1/3)*tau*epsilon_v*ri;

for k = 1:100
    Kt = zeros(nelem+1,nelem+1);
    Gk = zeros(nelem+1,1);
    %disp(u)
    %disp(rnodes(1,1))

    for e = 1:nelem
        Ae  = zeros(2,nelem+1);
        Ae(1,e) = 1;
        Ae(2,e+1) = 1;
        u_e = Ae*u;
        
        [Kt_e, Fint_e,epsilon_p_return,stress] = elementrout(u_e,element_r(e,:),E,neu,sigma_y,lambda,meu,epsilon_p(:,e),e);
        Kt = Kt + transpose(Ae)*Kt_e*Ae;
        Gk = Gk + transpose(Ae)*Fint_e;
        epsilon_p(:,e) = epsilon_p_return;
        stress_(:,e) = stress;
    end
    
    stress_rr = stress_(1,:);
    stress_pp = stress_(2,:);
    Kt_red = Kt(2:nelem+1,2:nelem+1);
    Gk_red = Gk(2:nelem+1);
    %%% Kt *delta_u = - G 
    delta_u = inv(Kt_red) * (- Gk_red) ;
    %disp(delta_u);
    u(2:nelem+1) = u(2:nelem+1) + delta_u ;
    %disp(size(Kt_red));
    %disp(size(Gk_red));
    %disp(u);
    %disp(epsilon_p)
    %disp(Gk)
    if (norm(delta_u,inf)<=0.005*norm(u,inf) | norm(Gk_red,inf)<=0.005*norm(Gk,inf))
        %disp(k);
        X = ['Element = ',num2str(nelem),' Full Load Volumatic strain = ',num2str(tau*epsilon_v),' and Uo = ',num2str(u(1,1)),' Number of iteration required to converge is ',num2str(k-1),'.'];
        disp(X);
        break;
    end
end

%disp(u_array);
%disp(stress_array);
x = rnodes;
y = u;
f1= figure;
plot(x,y,'-o');
ylabel('u(r)');
xlabel('r in [\mum] (nodes)');
legend('u(r) FEM at tau = 1');
title('u(r) at Final Load');
%disp(stress_rr);

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
stress_pp = stress1_(2,:);
r_middle = zeros(nelem,1);

for k=1:nelem
    r_middle(k,1) = (rnodes(k+1)+ rnodes(k))/2;
    %disp(r_middle);
    %[stress_rr_exact] = exact_solution_stress(r_middle(k,1),epsilon_v, ri, neu, E);
    %stress_array(k,1) = stress_rr_exact;
end
x1 = r_middle;
f2 = figure;
plot(x1,-stress_rr,'-o');
ylabel('Stress rr');
xlabel('r in [\mum] (element Gauss Point)');
legend('Stress rr FEM at tau = 1');
title('Stress rr at Final Load');

x1 = r_middle;
f3 = figure;
plot(x1,-stress_pp,'-o');
ylabel('Stress pp');
xlabel('r in [\mum] (element Gauss Point)');
legend('Stress pp FEM at tau = 1');
title('Stress pp at Final Load');
end