clear
clc
close all

C=0.2; %stabiltiy
cx=0.05; %viscotiy
N=161; %number of points
length=3; %length of nozzle
dx=length/(N-1); 
total_steps=7000;

load handel


rho=zeros(1,N); 
V=zeros(1,N);
T=zeros(1,N);

x=linspace(0, length, N); 
A=(1+2.2*(x-1.5).^2)/1;

%intial conditions t=0
rho((-1<x)&(x<0.5))=1;
rho((0.4999999999<x)&(x<1.5))=1-0.366*(x((0.4999999999<x)&(x<1.5))-0.5);
rho((1.499999<x)&(x<2.1))=0.634-0.702*(x((1.499999<x)&(x<2.1))-1.5);
rho((2.0999<x)&(x<3.1))=0.5892-0.10228*(x((2.0999<x)&(x<3.1))-2.1);

T((-1<x)&(x<0.5))=1;
T((0.4999999999<x)&(x<1.5))=1-0.167*(x((0.4999999999<x)&(x<1.5))-0.5);
T((1.499999<x)&(x<2.1))=0.833-0.4908*(x((1.499999<x)&(x<2.1))-1.5);
T((2.0999<x)&(x<3.1))=0.93968-0.0622*(x((2.0999<x)&(x<3.1))-2.1);

V=(0.59)./(rho.*A);


p = rho.*T;

g=1.4;


rho_throat=zeros(1,total_steps);
V_15=zeros(1,total_steps);
T_throat=zeros(1,total_steps);
P_throat=zeros(1,total_steps);

drdt_avg_throat=zeros(1,total_steps);
dVdt_avg_throat=zeros(1,total_steps);
dTdt_avg_throat=zeros(1,total_steps);

step=1;


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

        dVdt_pre(i)=-V(i)*(V(i+1)-V(i))/dx-1/g*((T(i+1)-T(i))/dx+T(i)/(rho(i))*(rho(i+1)-rho(i))/dx);

        dTdt_pre(i)=-V(i)*(T(i+1)-T(i))/dx-(g-1)*T(i)*((V(i+1)-V(i))/dx+V(i)*(log(A(i+1))-log(A(i)))/dx);

        S1(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(rho(i+1)-2*rho(i)+rho(i-1));
        S2(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(V(i+1)-2*V(i)+V(i-1));
        S3(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(T(i+1)-2*T(i)+T(i-1));

    end

    for i=2:N-1

        rho(i)=rho(i)+drdt_pre(i)*dt+S1(i);

        V(i)=V(i)+dVdt_pre(i)*dt+S2(i);

        T(i)=T(i)+dTdt_pre(i)*dt+S3(i);

    end

    p = rho.*T;

%MEXRI EDO EINAI SOSTOS KALA XRISTOUGENNA

    %corrector
    for i=2:N-1
        
        drdt_corr(i)=-rho(i)*(V(i)-V(i-1))/dx-rho(i)*V(i)*(log(A(i))-log(A(i-1)))/dx...
                  -V(i)*(rho(i)-rho(i-1))/dx;

        dVdt_corr(i)=-V(i)*(V(i)-V(i-1))/dx-1/g*((T(i)-T(i-1))/dx...
                  +T(i)/(rho(i))*(rho(i)-rho(i-1))/dx);

        dTdt_corr(i)=-V(i)*(T(i)-T(i-1))/dx-(g-1)*T(i)*((V(i)-V(i-1))/dx...
                  +V(i)*(log(A(i))-log(A(i-1)))/dx);


        S1(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(rho(i+1)-2*rho(i)+rho(i-1));
        S2(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(V(i+1)-2*V(i)+V(i-1));
        S3(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(T(i+1)-2*T(i)+T(i-1));
    end

    
    drdt_avg=0.5*(drdt_pre+drdt_corr);
    dVdt_avg=0.5*(dVdt_pre+dVdt_corr);
    dTdt_avg=0.5*(dTdt_pre+dTdt_corr);




    for i=2:N-1

        rho(i)=rho_old(i)+drdt_avg(i)*dt+S1(i);
        V(i)=V_old(i)+dVdt_avg(i)*dt+S2(i);
        T(i)=T_old(i)+dTdt_avg(i)*dt+S3(i);

    end

    %BOUNDARY CONDITIONS
    %EISODOS
    V(1)=2*V(2)-V(3); %floats

    %EKSODOS
    V(end)=2*V(end-1)-V(end-2);         %floats
    T(end)=2*T(end-1)-T(end-2);         %floats
    rho(end)=0.6784/T(end);   %floatsfind f


    M = V./(T.^0.5);  % Mach no
    p = rho.*T;       
    mass = rho.*V.*A;

    




    hold off

    rho_throat_time(step)=rho(find(A == 1));
    T_throat_time(step)=T(find(A == 1));
    P_throat_time(step)=p(find(A == 1));
    M_throat_time(step)=M(find(A == 1));

    drdt_avg_throat(step)=drdt_avg(find(A == 1));
    dVdt_avg_throat(step)=dVdt_avg(find(A == 1));
    dTdt_avg_throat(step)=dTdt_avg(find(A == 1));

    if step==total_steps
        rho_throat=rho(find(A == 1))
        T_throat=T(find(A == 1))
        P_throat=p(find(A == 1));
        M_throat=M(find(A == 1))
        m_throat=rho(find(A == 1)).*V(find(A == 1)).*A(find(A == 1));
        V_throat=V(find(A == 1));
    end

    if step==total_steps
        rho_exit=rho(end)
        T_exit=T(end)
        P_exit=p(end);
        M_exit=M(end)
        m_exit=rho(end).*V(end).*A(end);
        V_exit=V(end);
    end

    if step==1
        rho_first=rho;
        T_first=T;
        P_first=p;
        M_first=M;
        m_first=rho.*V.*A;
        V_first=V;
    end


    



    step=step+1;
end

sound(y,Fs)

%ANALYTICAL CALCULATIONS
pe_p01=0.6748;Ae_At=5.95;

g=1.4;target_value=pe_p01*Ae_At;

eqn=@(Me)(1/Me)*(2/(g+1))^((g+1)/(2*(g-1)))*(1+((g-1)/2)*Me^2)^(-0.5)-target_value;

Me=fzero(eqn,0.1);

pe_p0e=(1+(g-1)/2*Me^2)^(-3.5);

p02_p01=pe_p01/pe_p0e;

eqn=@(M1)((g+1)*M1^2/((g-1)*M1^2+2))^(g/(g-1))*((g+1)/(2*g*M1^2-(g-1)))^(1/(g-1))-p02_p01;

M1=fzero(eqn,2);

A_Athroat=sqrt((1/M1^2)*((2/(g+1))*(1+(g-1)/2*M1^2))^((g+1)/(g-1)));

x=linspace(0,3,N);A=(1+2.2*(x-1.5).^2)/1;

[~,minIndex]=min(abs(A-A_Athroat));

x_shock=x(minIndex);

p2_p1=1+(2*g/(g+1))*(M1^2-1);

M2=sqrt((1+((g-1)/2)*M1^2)/(g*M1^2-(g-1)/2));

M_analytical=zeros(N,2);

for i=1:minIndex
    if abs(A(i)-1)<1e-6
        M_analytical(i,1)=1;
        M_analytical(i,2)=1;
    else
        eqn=@(M)(A(i)^2)-(1/M^2)*((2/(g+1))*(1+((g-1)/2)*M^2))^((g+1)/(g-1));
        M_analytical(i,1)=fzero(eqn,0.5);
        M_analytical(i,2)=fzero(eqn,1.5);
    end
end

M_correct=zeros(1,N);

for i=1:minIndex
    if i<N/2
        M_correct(i)=M_analytical(i,1);
    elseif i>N/2
        M_correct(i)=M_analytical(i,2);
    else
        M_correct(i)=1;
    end
end

for i=1:minIndex
    M_calc=M_correct(i);
end

A_Astar=A*0.6882;

for i=minIndex+1:N
    eqn=@(M)(A_Astar(i)^2)-(1/M^2)*((2/(g+1))*(1+((g-1)/2)*M^2))^((g+1)/(g-1));
    M_sub(i)=fzero(eqn,0.8);
end

for i=minIndex+1:N
    M_correct(i)=M_sub(i);
end

p_analytical=zeros(1,N);rho_analytical=zeros(1,N);T_analytical=zeros(1,N);

for i=1:N
    M_calc=M_correct(i);
    p_analytical(i)=(1+(g-1)/2*M_calc^2)^(-g/(g-1));
    rho_analytical(i)=(1+(g-1)/2*M_calc^2)^(-1/(g-1));
    T_analytical(i)=(1+(g-1)/2*M_calc^2)^(-1);
end

%END OF ANALYTICAL CALCULATIOSN

hold off
figure(3)
subplot(4,1,1);
plot(rho_throat_time)
title('Non-dimensionalized ρ at point 15');
xlabel('Time step')
grid on;
subplot(4,1,2);
plot(T_throat_time)
title('Non-dimensionalized T at point 15'); 
xlabel('Time step')
grid on;
subplot(4,1,3);
plot(P_throat_time)
title('Non-dimensionalized P at point 15'); 
xlabel('Time step')
grid on;
subplot(4,1,4);
plot(M_throat_time)
title('Non-dimensionalized M at point 15'); 
xlabel('Time step')
grid on;



figure(4)
yyaxis left
plot(x,rho)
hold on
plot(x,rho_analytical)
hold off
title("Steady nondimensional ρ & Ma: " + ...
    "analytical (circles) & numerical results (curves).")
xlabel('xlabel')
ylabel('ρ', 'Rotation', 0)
%ylim([0 1]); % Adjust limits for better visualization
yyaxis right
plot(x,M)
hold on
plot(x,M_correct)
hold off
ylabel('Mach  ', 'Rotation', 0)
%ylim([0 max(M)]); % Adjust limits for better visualization


results_table = table((1:N)', x', A', rho', V', T', p', M', mass', ...
    'VariableNames', {'I', 'x/L', 'A/A*', 'rho/rho0', 'V/a0', 'T/T0', 'P/P0', 'M', 'm_dot'});


first_results_table = table((1:N)', x', A', rho_first', V_first', T_first', P_first', M_first', m_first',  ...
    'VariableNames', {'I', 'x/L', 'A/A*', 'rho/rho0', 'V/a0', 'T/T0', 'P/P0', 'M', 'm_dot'});



figure(55)
plot(x,mass)
title('Mass flow')
grid on

figure(6)
plot(x,p)
title('Pressure')
grid on


figure(1)
plot(drdt_avg_throat ,'b')
hold on
plot(dVdt_avg_throat,'r')
plot(dTdt_avg_throat,'g')
hold off
legend('\partialρ/\partialt', '\partialV/\partialt', '\partialT/\partialt', 'Location', 'best');
title('Residuals')
grid on
xlabel('Time step')












