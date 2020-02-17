function [t,T,Ts,Vi,Vni] = multistep_TTsViVni(t_min,t_max,N,lambda,d,kappa,k,delta,eta,N_T,c,T0,Ts0,Vi0,Vni0)
% try: multistep_TTsViVni(0,300,1000,100,0.01,0.15,0.00000065,0.39,0.1,3127.5,3,10000,0,0.0000001,0.000000001);

% set step size and initial matrices for values
h = (t_max-t_min)/N;
t = zeros(N+1,1);
T = zeros(N+1,1);
Ts = zeros(N+1,1);
Vi = zeros(N+1,1);
Vni = zeros(N+1,1);

% set initial time and initial conditions
t(1) = t_min;
T(1) = T0;
Ts(1) = Ts0;
Vi(1) = Vi0;
Vni(1) = Vni0;

% use Runge-Kutta method order 4 for initial values for multistep method
for i=1:4

    t(i+1) = t_min+i*h;
    
    K1_T = h*(lambda - d*T(i) - (1-kappa)*k*Vi(i)*T(i));
    K1_Ts = h*((1-kappa)*k*Vi(i)*T(i) - delta*Ts(i));
    K1_Vi = h*((1-eta)*N_T*delta*Ts(i) - c*Vi(i));
    K1_Vni = h*(eta*N_T*delta*Ts(i) - c*Vni(i));
    
    K2_T = h*(lambda - d*(T(i)+(1/2)*K1_T) - (1-kappa)*k*(Vi(i)+(1/2)*K1_Vi)*(T(i)+(1/2)*K1_T));
    K2_Ts = h*((1-kappa)*k*(Vi(i)+(1/2)*K1_Vi)*(T(i)+(1/2)*K1_T) - delta*(Ts(i)+(1/2)*K1_Ts));
    K2_Vi = h*((1-eta)*N_T*delta*(Ts(i)+(1/2)*K1_Ts) - c*(Vi(i)+(1/2)*K1_Vi));
    K2_Vni = h*(eta*N_T*delta*(Ts(i)+(1/2)*K1_Ts) - c*(Vni(i)+(1/2)*K1_Vni));
    
    K3_T = h*(lambda - d*(T(i)+(1/2)*K2_T) - (1-kappa)*k*(Vi(i)+(1/2)*K2_Vi)*(T(i)+(1/2)*K2_T));
    K3_Ts = h*((1-kappa)*k*(Vi(i)+(1/2)*K2_Vi)*(T(i)+(1/2)*K2_T) - delta*(Ts(i)+(1/2)*K2_Ts));
    K3_Vi = h*((1-eta)*N_T*delta*(Ts(i)+(1/2)*K2_Ts) - c*(Vi(i)+(1/2)*K2_Vi));
    K3_Vni = h*(eta*N_T*delta*(Ts(i)+(1/2)*K2_Ts) - c*(Vni(i)+(1/2)*K2_Vni));
    
    K4_T = h*(lambda - d*(T(i)+K3_T) - (1-kappa)*k*(Vi(i)+K3_Vi)*(T(i)+K3_T));
    K4_Ts = h*((1-kappa)*k*(Vi(i)+K3_Vi)*(T(i)+K3_T) - delta*(Ts(i)+K3_Ts));
    K4_Vi = h*((1-eta)*N_T*delta*(Ts(i)+K3_Ts) - c*(Vi(i)+K3_Vi));
    K4_Vni = h*(eta*N_T*delta*(Ts(i)+K3_Ts) - c*(Vni(i)+K3_Vni));
    
    T(i+1) = T(i) + (1/6)*(K1_T + 2*K2_T + 2*K3_T + K4_T);
    Ts(i+1) = Ts(i) + (1/6)*(K1_Ts + 2*K2_Ts + 2*K3_Ts + K4_Ts);
    Vi(i+1) = Vi(i) + (1/6)*(K1_Vi + 2*K2_Vi + 2*K3_Vi + K4_Vi);
    Vni(i+1) = Vni(i) + (1/6)*(K1_Vni + 2*K2_Vni + 2*K3_Vni + K4_Vni);
    
end

% use multistep method to solve system numerically
for i=5:N
    
    t(i+1) = t_min+i*h;
    
    T(i+1) = T(i) + (h/24)*(55*(lambda-d*T(i)-(1-kappa)*k*Vi(i)*T(i)) - 59*(lambda-d*T(i-1)-(1-kappa)*k*Vi(i-1)*T(i-1)) + 37*(lambda-d*T(i-2)-(1-kappa)*k*Vi(i-2)*T(i-2)) - 9*(lambda-d*T(i-3)-(1-kappa)*k*Vi(i-3)*T(i-3)));
    Ts(i+1) = Ts(i) + (h/24)*(55*((1-kappa)*k*Vi(i)*T(i)-delta*Ts(i)) - 59*((1-kappa)*k*Vi(i-1)*T(i-1)-delta*Ts(i-1)) + 37*((1-kappa)*k*Vi(i-2)*T(i-2)-delta*Ts(i-2)) - 9*((1-kappa)*k*Vi(i-3)*T(i-3)-delta*Ts(i-3)));
    Vi(i+1) = Vi(i) + (h/24)*(55*((1-eta)*N_T*delta*Ts(i)-c*Vi(i)) - 59*((1-eta)*N_T*delta*Ts(i-1)-c*Vi(i-1)) + 37*((1-eta)*N_T*delta*Ts(i-2)-c*Vi(i-2)) - 9*((1-eta)*N_T*delta*Ts(i-3)-c*Vi(i-3)));
    Vni(i+1) = Vni(i) + (h/24)*(55*(eta*N_T*delta*Ts(i)-c*Vni(i)) - 59*(eta*N_T*delta*Ts(i-1)-c*Vni(i-1)) + 37*(eta*N_T*delta*Ts(i-2)-c*Vni(i-2)) - 9*(eta*N_T*delta*Ts(i-3)-c*Vni(i-3)));
   
    T(i+1) = T(i) + (h/24)*(9*(lambda-d*T(i+1)-(1-kappa)*k*Vi(i+1)*T(i+1)) + 19*(lambda-d*T(i)-(1-kappa)*k*Vi(i)*T(i)) - 5*(lambda-d*T(i-1)-(1-kappa)*k*Vi(i-1)*T(i-1)) + (lambda-d*T(i-2)-(1-kappa)*k*Vi(i-2)*T(i-2)));
    Ts(i+1) = Ts(i) + (h/24)*(9*((1-kappa)*k*Vi(i+1)*T(i+1)-delta*Ts(i+1)) + 19*((1-kappa)*k*Vi(i)*T(i)-delta*Ts(i)) - 5*((1-kappa)*k*Vi(i-1)*T(i-1)-delta*Ts(i-1)) + ((1-kappa)*k*Vi(i-2)*T(i-2)-delta*Ts(i-2)));
    Vi(i+1) = Vi(i) + (h/24)*(9*((1-eta)*N_T*delta*Ts(i+1)-c*Vi(i+1)) + 19*((1-eta)*N_T*delta*Ts(i)-c*Vi(i)) - 5*((1-eta)*N_T*delta*Ts(i-1)-c*Vi(i-1)) + ((1-eta)*N_T*delta*Ts(i-2)-c*Vi(i-2)));
    Vni(i+1) = Vni(i) + (h/24)*(9*(eta*N_T*delta*Ts(i+1)-c*Vni(i+1)) + 19*(eta*N_T*delta*Ts(i)-c*Vni(i)) - 5*(eta*N_T*delta*Ts(i-1)-c*Vni(i-1)) + (eta*N_T*delta*Ts(i-2)-c*Vni(i-2)));
end

% plot the approximate solution for T and Ts
% in cells per cubic millimeter
subplot(2,1,1)
plot(t,T,'-b');
hold on;
plot(t,Ts,'-g');
legend('T','T*');
title('Numerical Susceptible and Infected','FontSize',12);
xlabel('Time (days)','FontSize',12);
ylabel('Class Size','FontSize',12);

% plot the approximate solution for VI and VNI
% in virus particles per cubic millimeter
subplot(2,1,2);
plot(t,Vi,'-r');
hold on;
plot(t,Vni,'-m');
legend('VI','VNI');
title('Numerical Infectious Virus and Non-Infectious Virus','FontSize',12);
xlabel('Time (days)','FontSize',12);
ylabel('Class Size','FontSize',12);