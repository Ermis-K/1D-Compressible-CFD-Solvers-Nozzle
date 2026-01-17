clear
clc
close all
tic

C=1.1; %stabiltiy
N=81; %number of points
length=3; %length of nozzle
total_steps=5000;

dx=length/(N-1);
x=linspace(0, length, N); %
A = 1 + (2.2 * (x-1.5).^2).*(x <= 1.5) + (0.2223 * (x-1.5).^2).*(x > 1.5);
g=1.4;
step=1;

load handel

%initial conditions
rho=1-0.023*x;
T=1-0.009333*x;
V=0.05+0.11*x;

U1=rho.*A;
U2=rho.*A.*V;
U3=rho.*(T/0.4+g/2.*V.^2).*A;
P = rho.*T;


figure(1)
h4=plot(NaN, NaN, 'b', 'LineWidth', 1.5);
hold on
h5=plot(NaN, NaN, 'r', 'LineWidth', 1.5);
h6=plot(NaN, NaN, 'g', 'LineWidth', 1.5);
legend('\partialU1/\partialt', '\partialU2/\partialt', '\partialU3/\partialt', 'Location', 'best');
title('Residuals')
grid on
xlabel('Time step')



while step<total_steps+1
    step;
    
    
    a=sqrt(T);
    dt=min(C*dx./(a+V));

    U1_old = U1;
    U2_old = U2;
    U3_old = U3;

    %PREDICTOR STEP
    F1=U2;
    F2=((U2.^2)./U1)+((g - 1)/g)*(U3 - (g/2).*((U2.^2)./U1));
    F3=g*(U2.*U3)./U1-g*(g-1)/2*U2.^3./U1.^2;
    %J2=(g-1)/g*(U3-g/2*U2.^2./U1)*
    
    %J2=1/g*rho(1:end-1).*T(1:end-1).*(A(2:end)-A(1:end-1))./dx;



    for i=2:N-1
        J2(i) = (1/g)*rho(i).*T(i).*((A(i+1) - A(i))/dx);

        dU1dt_pre(i)=-(F1(i+1)-F1(i))/dx;
        dU2dt_pre(i)=-(F2(i+1)-F2(i))/dx+J2(i);
        dU3dt_pre(i)=-(F3(i+1)-F3(i))/dx;


        U1(i)=U1(i)+dU1dt_pre(i)*dt;
        U2(i)=U2(i)+dU2dt_pre(i)*dt;
        U3(i)=U3(i)+dU3dt_pre(i)*dt;
    end

    rho=U1./A;
    V=U2./U1;
    T=(g- 1)*(((U3)./(U1)) - (g/2)*((V).^2));


    F1=U2;
    F2=((U2.^2)./U1)+((g - 1)/g)*(U3 - (g/2).*((U2.^2)./U1));
    F3=g*(U2.*U3)./U1-g*(g-1)/2*U2.^3./U1.^2;


    %CORRECTOR STEP
    for i=2:N-1
        J2(i) = (1/g)*rho(i).*T(i).*((A(i) - A(i-1))/dx);

        dU1dt_cor(i)=-(F1(i)-F1(i-1))/dx;
        dU2dt_cor(i)=-(F2(i)-F2(i-1))/dx+J2(i);
        dU3dt_cor(i)=-(F3(i)-F3(i-1))/dx;
    end

    dU1dt_avg=0.5*(dU1dt_pre+dU1dt_cor);
    dU2dt_avg=0.5*(dU2dt_pre+dU2dt_cor);
    dU3dt_avg=0.5*(dU3dt_pre+dU3dt_cor);


    U1(2:N-1)=U1_old(2:N-1)+dU1dt_avg(2:N-1)*dt;
    U2(2:N-1)=U2_old(2:N-1)+dU2dt_avg(2:N-1)*dt;
    U3(2:N-1)=U3_old(2:N-1)+dU3dt_avg(2:N-1)*dt;



    %boundary conditions
    %inflow
    U2(1)=2*U2(2)-U2(3);
    U3(1)=U1(1)*(T(1)/0.4+0.7*V(1)^2);
 
    %outflow
    U1(end)=2*U1(end-1)-U1(end-2);
    U2(end)=2*U2(end-1)-U2(end-2);
    V(end)=U2(end)/U1(end);
    T(end)=0.93/rho(end);
    U3(end)=U1(end)*(T(end)/(g-1)+g/2*V(end)^2);





    rho=U1./A;
    V=U2./U1;
    T=(g-1)*(U3./U1-g/2*V.^2);
    P=rho.*T;


    M = V./(T.^0.5);  % Mach no      
    mass = rho.*V.*A;

    rho_15(step)=rho(15);
    T_15(step)=T(15);
    P_15(step)=P(15);
    M_15(step)=M(15);

    dU1dt_avg_16(step)=dU1dt_avg(16);
    dU2dt_avg_16(step)=dU2dt_avg(16);
    dU3dt_avg_16(step)=dU3dt_avg(16);


    figure(5);
    hold on
    
    if step == 50
        plot(x,mass,'color','g','linewidth',1.5)
    elseif step == 100
        plot(x,mass,'color','m','linewidth',1.5)
    elseif step == 150
        plot(x,mass,'color','b','linewidth',1.5)
    elseif step == 200
        plot(x,mass,'color','r','linewidth',1.5)
    elseif step == 700
        plot(x,mass,'color','c','linewidth',1.5)
        
        ylabel('Mass Flow (Nondimensional)')
        grid on
        legend('50 dt','100 dt','150 dt','200 dt','700 dt')
        title({['Non Conservative'];['Mass flow through nozzle at different dt']})
    end

    hold off

    figure(6);
    hold on
    
    if step == 50
        plot(x,P,'color','g','linewidth',1.5)
    elseif step == 100
        plot(x,P,'color','m','linewidth',1.5)
    elseif step == 150
        plot(x,P,'color','b','linewidth',1.5)
    elseif step == 200
        plot(x,P,'color','r','linewidth',1.5)
    elseif step == 700
        plot(x,P,'color','c','linewidth',1.5)
        
        ylabel('Pressure (Nondimensional)')
        grid on
        legend('50 dt','100 dt','150 dt','200 dt','700 dt')
        title({['Conservative'];['Pressure through nozzle at different dt']})
    end

    hold off


    set(h4, 'YData', dU1dt_avg_16(1,1:step), 'XData', 1:step);
    set(h5, 'YData', dU2dt_avg_16(1,1:step), 'XData', 1:step);
    set(h6, 'YData', dU3dt_avg_16(1,1:step), 'XData', 1:step);

    step=step+1;


    
    drawnow


end


M_analytical = zeros(N, 1); 
Me =sqrt((2 /(g-1)) * ((0.93)^(- (g-1) / g)- 1));
Ae_Astar=(1/ Me) * ((2 / (g+1))*(1+ ((g -1) / 2) *Me^2))^((g+ 1) / (2 * (g- 1)));
Astar=(Ae_Astar)^-1*A(end);
A_Astar=A/Astar;

for i = 1:N
        eqn = @(M) (A_Astar(i)^2) - (1 / M^2) * ((2 / (g + 1)) * (1 + ((g - 1) / 2) * M^2))^((g + 1) / (g - 1));
        M_analytical(i) = fzero(eqn, 0.5); %n
end

M_correct = M_analytical; 



p_analytical = zeros(1, N);
rho_analytical = zeros(1, N);
T_analytical = zeros(1, N);


for i = 1:N
    M_calc = M_correct(i);
    p_analytical(i) = (1 + (g - 1) / 2 * M_calc^2)^(-g / (g - 1)); 
    rho_analytical(i) = (1 + (g - 1) / 2 * M_calc^2)^(-1 / (g - 1)); 
    T_analytical(i) = (1 + (g - 1) / 2 * M_calc^2)^(-1); % 
end


sound(y,Fs)

%hold off
figure(3)
subplot(4,1,1);
plot(rho_15)
title('Non-dimensionalized ρ at point 15');
xlabel('Time step')
grid on;
subplot(4,1,2);
plot(T_15)
title('Non-dimensionalized T at point 15'); 
xlabel('Time step')
grid on;
subplot(4,1,3);
plot(P_15)
title('Non-dimensionalized P at point 15'); 
xlabel('Time step')
grid on;
subplot(4,1,4);
plot(M_15)
title('Non-dimensionalized M at point 15'); 
xlabel('Time step')
grid on;



figure(4)
yyaxis left
plot(x,rho)
hold on
plot(x,rho_analytical,'o')
hold off
title("Steady nondimensional ρ & Ma: " + ...
    "analytical (curves) & numerical results (point).")
xlabel('xlabel')
ylabel('ρ', 'Rotation', 0)
ylim([0 1]); 
yyaxis right
plot(x,M)
hold on
plot(x,M_correct,'x')
hold off
ylabel('Mach  ', 'Rotation', 0)
ylim([0 max(M)]);


results_table = table((1:N)', x', A', rho', V', T', P', M', mass', U1' ,U2', U3', ...
    'VariableNames', {'I', 'x/L', 'A/A*', 'rho/rho0', 'V/a0', 'T/T0', 'P/P0', 'M', 'm_dot', 'U1', 'U2', 'U3'});
disp(results_table);

figure(5)
hold on
plot(x,mass,'color','bla','linewidth',1.5)
legend('50 dt','100 dt','150 dt','200 dt','700 dt',['Final Step: ', num2str(step-1)])
hold off

figure(6)
hold on
plot(x,P,'color','bla','linewidth',1.5)
legend('50 dt','100 dt','150 dt','200 dt','700 dt',['Final Step: ', num2str(step-1)])
hold off

Diff_rho_per=(rho-rho_analytical)*100;
Diff_M_per=(M-M_correct')*100;
results_compare = table((1:N)', x', rho', rho_analytical', Diff_rho_per', M', M_correct, Diff_M_per', ...
    'VariableNames', {'I', 'x/L', 'ρ numerical', 'ρ analytical', 'Diff rho %', 'M numerical', 'M analytical', 'Diff M %'});
disp(results_compare);

fprintf('Final step: %.2f\n', step-1)

toc

rho(A==1)
T(A==1)
P(A==1)
M(A==1)
rho_analytical(A==1)
T_analytical(A==1)
p_analytical(A==1)
M_correct(A==1)