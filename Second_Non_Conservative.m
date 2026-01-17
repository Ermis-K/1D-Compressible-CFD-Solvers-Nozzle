clear
clc
close all

C=1.1; %stabiltiy 
N=81; %number of points
length=3; %length of nozzle
dx=length/(N-1); 
total_steps=7000;

load handel

rho=zeros(1,N); 
V=zeros(1,N);
T=zeros(1,N);

x=linspace(0, length, N); 
A = 1 + (2.2 * (x-1.5).^2).*(x <= 1.5) + (0.2223 * (x-1.5).^2).*(x > 1.5);

%intial conditions t=0
rho=1-0.023*x;
T=1-0.009333*x;
V=0.05+0.11*x;


g=1.4;


rho_15=zeros(1,total_steps);
V_15=zeros(1,total_steps);
T_15=zeros(1,total_steps);
P_15=zeros(1,total_steps);

drdt_avg_16=zeros(1,total_steps);
dVdt_avg_16=zeros(1,total_steps);
dTdt_avg_16=zeros(1,total_steps);


step=1;

figure(55)
h4=plot(NaN, NaN, 'b', 'LineWidth', 1.5);
hold on
h5=plot(NaN, NaN, 'r', 'LineWidth', 1.5);
h6=plot(NaN, NaN, 'g', 'LineWidth', 1.5);
legend('\partial\rho/\partialt', '\partialV/\partialt', '\partialT/\partialt', 'Location', 'best');
title('Residuals')
grid on
xlabel('Time step')
hold off


while step<total_steps+1
    step;
    
    
    a=sqrt(T);
    dt=min(C*dx./(a+V));

    rho_old = rho;
    V_old = V;
    T_old = T;


    %predictor step
    for i=2:N-1
        drdt_pre(i)=-rho(i)*(V(i+1)-V(i))/dx-rho(i)*V(i)*(log(A(i+1))-log(A(i)))/dx...
                -V(i)*(rho(i+1)-rho(i))/dx;

        dVdt_pre(i)=-V(i)*(V(i+1)-V(i))/dx-1/g*((T(i+1)-T(i))/dx+T(i)/rho(i)*(rho(i+1)-rho(i))/dx);

        dTdt_pre(i)=-V(i)*(T(i+1)-T(i))/dx-(g-1)*T(i)*((V(i+1)-V(i))/dx+V(i)*(log(A(i+1))-log(A(i)))/dx);


        rho(i)=rho(i)+drdt_pre(i)*dt;

        V(i)=V(i)+dVdt_pre(i)*dt;

        T(i)=T(i)+dTdt_pre(i)*dt;


    end

%MEXRI EDO EINAI SOSTOS KALA XRISTOUGENNA

    %corrector
    for i=2:N-1
        
        drdt_corr(i)=-rho(i)*(V(i)-V(i-1))/dx-rho(i)*V(i)*(log(A(i))-log(A(i-1)))/dx...
                  -V(i)*(rho(i)-rho(i-1))/dx;

        dVdt_corr(i)=-V(i)*(V(i)-V(i-1))/dx-1/g*((T(i)-T(i-1))/dx...
                  +T(i)/rho(i)*(rho(i)-rho(i-1))/dx);

        dTdt_corr(i)=-V(i)*(T(i)-T(i-1))/dx-(g-1)*T(i)*((V(i)-V(i-1))/dx...
                  +V(i)*(log(A(i))-log(A(i-1)))/dx);

    end

    
    drdt_avg=0.5*(drdt_pre+drdt_corr);
    dVdt_avg=0.5*(dVdt_pre+dVdt_corr);
    dTdt_avg=0.5*(dTdt_pre+dTdt_corr);


    for i=2:N-1

        rho(i)=rho_old(i)+drdt_avg(i)*dt;
        V(i)=V_old(i)+dVdt_avg(i)*dt;
        T(i)=T_old(i)+dTdt_avg(i)*dt;

    end

    %BOUNDARY CONDITIONS
    %EISODOS
    V(1)=2*V(2)-V(3); %floats

    %EKSODOS
    V(end)=2*V(end-1)-V(end-2);         %floats
    T(end)=2*T(end-1)-T(end-2);         %floats
    rho(end)=0.93/T(end);   %floatsfind f


    M = V./(T.^0.5);  % Mach no
    P = rho.*T;       
    mass = rho.*V.*A;


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
        title({['Conservative'];['Mass flow through nozzle at different dt']})
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

    rho_15(step)=rho(15);
    T_15(step)=T(15);
    P_15(step)=P(15);
    M_15(step)=M(15);

    drdt_avg_16(step)=drdt_avg(16);
    dVdt_avg_16(step)=dVdt_avg(16);
    dTdt_avg_16(step)=dTdt_avg(16);

    set(h4, 'YData', drdt_avg_16(1,1:step), 'XData', 1:step);
    set(h5, 'YData', dVdt_avg_16(1,1:step), 'XData', 1:step);
    set(h6, 'YData', dTdt_avg_16(1,1:step), 'XData', 1:step);
    

    drawnow

    step=step+1;
end

M_analytical = zeros(N, 1); 
Me =sqrt((2 /(g-1)) * ((0.93)^(- (g-1) / g)- 1));
Ae_Astar=(1/ Me) * ((2 / (g+1))*(1+ ((g -1) / 2) *Me^2))^((g+ 1) / (2 * (g- 1)));
Astar=(Ae_Astar)^-1*A(end);
A_Astar=A/Astar;



for i = 1:N
        eqn = @(M) (A_Astar(i)^2) - (1 / M^2) * ((2 / (g + 1)) * (1 + ((g - 1) / 2) * M^2))^((g + 1) / (g - 1));
        M_analytical(i) = fzero(eqn, 0.5); % S
end


M_correct = M_analytical; % 



p_analytical = zeros(1, N);
rho_analytical = zeros(1, N);
T_analytical = zeros(1, N);


for i = 1:N
    M_calc = M_correct(i); 
    p_analytical(i) = (1 + (g - 1) / 2 * M_calc^2)^(-g / (g - 1));
    rho_analytical(i) = (1 + (g - 1) / 2 * M_calc^2)^(-1 / (g - 1)); 
    T_analytical(i) = (1 + (g - 1) / 2 * M_calc^2)^(-1);
end


sound(y,Fs)

hold off
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
hold on
plot(x,rho)
plot(x,rho_analytical,'o')
hold off
title("Steady nondimensional ρ & Ma: " + ...
    "analytical (circles) & numerical results (curves).")
xlabel('xlabel')
ylabel('ρ', 'Rotation', 0)
%ylim([0 1]); % Adjust limits for better visualization
yyaxis right
hold on
plot(x,M)
plot(x,M_correct,'x')
hold off
ylabel('Mach  ', 'Rotation', 0)
%ylim([0 max(M)]); % Adjust limits for better visualization


results_table = table((1:N)', x', A', rho', V', T', P', M', mass', ...
    'VariableNames', {'I', 'x/L', 'A/A*', 'rho/rho0', 'V/a0', 'T/T0', 'P/P0', 'M', 'm_dot'});
disp(results_table);


figure(6)
hold on
plot(x,P,'color','bla','linewidth',1.5)
legend('50 dt','100 dt','150 dt','200 dt','700 dt',['Final Step: ', num2str(step-1)])
hold off

figure(5)
hold on
plot(x,mass,'color','bla','linewidth',1.5)
legend('50 dt','100 dt','150 dt','200 dt','700 dt',['Final Step: ', num2str(step-1)])
hold off



Diff_rho_per=(rho-rho_analytical)*100;
Diff_M_per=(M-M_correct')*100;
results_compare = table((1:N)', x', rho', rho_analytical', Diff_rho_per', M', M_correct, Diff_M_per', ...
    'VariableNames', {'I', 'x/L', 'ρ numerical', 'ρ analytical', 'Diff rho %', 'M numerical', 'M analytical', 'Diff M %'});
disp(results_compare);

fprintf('Final step: %.2f\n', step-1)


rho(A==1)
T(A==1)
P(A==1)
M(A==1)
rho_analytical(A==1)
T_analytical(A==1)
p_analytical(A==1)
M_correct(A==1)





















