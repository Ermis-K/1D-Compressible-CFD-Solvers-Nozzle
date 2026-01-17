clear
clc
close all
tic

C=0.5; %stabiltiy
total_steps=10000;
cx=0.0; %0.01-0.3
N=161; %number of points

length=3; %length of nozzle
dx=length/(N-1);
x=linspace(0, length, N); 
A=(1+2.2*(x-1.5).^2)/1;
g=1.4;
step=1;


load handel

%initial conditions
rho((-1<x)&(x<0.5))=1;
rho((0.4999999999<x)&(x<1.5))=1-0.366*(x((0.4999999999<x)&(x<1.5))-0.5);
rho((1.499999<x)&(x<2.1))=0.634-0.702*(x((1.499999<x)&(x<2.1))-1.5);
rho((2.0999<x)&(x<3.1))=0.5892-0.10228*(x((2.0999<x)&(x<3.1))-2.1);

T((-1<x)&(x<0.5))=1;
T((0.4999999999<x)&(x<1.5))=1-0.167*(x((0.4999999999<x)&(x<1.5))-0.5);
T((1.499999<x)&(x<2.1))=0.833-0.4908*(x((1.499999<x)&(x<2.1))-1.5);
T((2.0999<x)&(x<3.1))=0.93968-0.0622*(x((2.0999<x)&(x<3.1))-2.1);

V=(0.59)./(rho.*A);

U1=rho.*A;
U2=rho.*A.*V;
U3=rho.*(T/0.4+g/2.*V.^2).*A;
p = rho.*T;

a=sqrt(T);
dt=min(C*dx./(a+V));

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


        S1(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(U1(i+1)-2*U1(i)+U1(i-1));
        S2(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(U2(i+1)-2*U2(i)+U2(i-1));
        S3(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(U3(i+1)-2*U3(i)+U3(i-1));



    end

    for i=2:N-1
        U1(i)=U1(i)+dU1dt_pre(i)*dt+S1(i);
        U2(i)=U2(i)+dU2dt_pre(i)*dt+S2(i);
        U3(i)=U3(i)+dU3dt_pre(i)*dt+S3(i);
    end

    rho=U1./A;
    V=U2./U1;
    T=(g- 1)*(((U3)./(U1)) - (g/2)*((V).^2));
    p=rho.*T;
    

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



    for i=2:N-1
        S1(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(U1(i+1)-2*U1(i)+U1(i-1));
        S2(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(U2(i+1)-2*U2(i)+U2(i-1));
        S3(i)=cx*abs(p(i+1)-2*p(i)+p(i-1))/(p(i+1)+2*p(i)+p(i-1))*(U3(i+1)-2*U3(i)+U3(i-1));
    end


    U1(2:N-1)=U1_old(2:N-1)+dU1dt_avg(2:N-1)*dt+S1(2:N-1);
    U2(2:N-1)=U2_old(2:N-1)+dU2dt_avg(2:N-1)*dt+S2(2:N-1);
    U3(2:N-1)=U3_old(2:N-1)+dU3dt_avg(2:N-1)*dt+S3(2:N-1);

    %boundary conditions
    %U(1)=rho(1)*a(1);
    U2(1)=2*U2(2)-U2(3);
    U3(1)=U1(1)*(T(1)/0.4+0.7*V(1)^2);

    U1(end)=2*U1(end-1)-U1(end-2);
    U2(end)=2*U2(end-1)-U2(end-2);
    V=U2./U1;
    U3(end)=0.6784*A(end)/(g-1)+g/2*U2(end)*V(end);



    rho=U1./A;
    T=(g-1)*(U3./U1-g/2*V.^2);
    p=rho.*T;

    M = V./(T.^0.5);  % Mach no      
    mass = rho.*V.*A;

    if step==total_steps
        rho_throat=rho(find(A == 1))
        T_throat=T(find(A == 1))
        P_throat=p(find(A == 1))
        M_throat=M(find(A == 1))
        m_throat=rho(find(A == 1)).*V(find(A == 1)).*A(find(A == 1))
        V_throat=V(find(A == 1))
        U1_throat=U1(find(A == 1))
        U2_throat=U2(find(A == 1))
        U3_throat=U3(find(A == 1))
    end

    if step==total_steps
        rho_exit=rho(end)
        T_exit=T(end)
        P_exit=p(end)
        M_exit=M(end)
        m_exit=rho(end).*V(end).*A(end)
        V_exit=V(end)
        U1_exit=U1(end)
        U2_exit=U2(end)
        U3_exit=U3(end)
    end

    rho_throat_time(step)=rho(find(A == 1));
    T_throat_time(step)=T(find(A == 1));
    P_throat_time(step)=p(find(A == 1));
    M_throat_time(step)=M(find(A == 1));

    dU1dt_avg_throat(step)=dU1dt_avg(find(A == 1));
    dU2dt_avg_throat(step)=dU2dt_avg(find(A == 1));
    dU3dt_avg_throat(step)=dU3dt_avg(find(A == 1));


    if step==1
        U1_first=U1;
        U2_first=U2;
        U3_first=U3;
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

A_Astar=A*p02_p01;

for i=minIndex+1:N
    eqn=@(M)(A_Astar(i)^2)-(1/M^2)*((2/(g+1))*(1+((g-1)/2)*M^2))^((g+1)/(g-1));
    M_sub(i)=fzero(eqn,0.566);
end

for i=minIndex+1:N
    M_correct(i)=M_sub(i);
end

p_analytical=zeros(1,N);rho_analytical=zeros(1,N);T_analytical=zeros(1,N);

for i=1:N
    p_analytical(i)=(1+(g-1)/2*M_correct(i)^2)^(-g/(g-1));
    rho_analytical(i)=(1+((g-1)/2)*M_correct(i)^2)^(-1/(g-1));
    T_analytical(i)=(1+(g-1)/2*M_correct(i)^2)^(-1);
end

for i=minIndex+1:N
    p_analytical(i)=p_analytical(i)*p02_p01;
    rho_analytical(i)=rho_analytical(i)*p02_p01;
end

%END OF ANALYTICAL CALCULATIOSN





figure(3)
subplot(4,1,1);
plot(rho_throat_time)
title('Non-dimensionalized ρ at throat');
xlabel('Time step')
grid on;
subplot(4,1,2);
plot(T_throat_time)
title('Non-dimensionalized T at throat'); 
xlabel('Time step')
grid on;
subplot(4,1,3);
plot(P_throat_time)
title('Non-dimensionalized P at throat'); 
xlabel('Time step')
grid on;
subplot(4,1,4);
plot(M_throat_time)
title('Non-dimensionalized M at throat'); 
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
%ylim([0 1]); 
yyaxis right
plot(x,M)
hold on
plot(x,M_correct)
hold off
ylabel('Mach  ', 'Rotation', 0)
%ylim([0 max(M)]);

results_table = table((1:N)', x', A', rho', V', T', p', M', mass', U1' ,U2', U3', ...
    'VariableNames', {'I', 'x/L', 'A/A*', 'rho/rho0', 'V/a0', 'T/T0', 'P/P0', 'M', 'm_dot', 'U1', 'U2', 'U3'});

first_results_table = table((1:N)', x', A', rho_first', V_first', T_first', P_first', M_first', m_first', U1_first' ,U2_first', U3_first', ...
    'VariableNames', {'I', 'x/L', 'A/A*', 'rho/rho0', 'V/a0', 'T/T0', 'P/P0', 'M', 'm_dot', 'U1', 'U2', 'U3'});

figure(5)
plot(x,mass)
title('Mass flow')
grid on

figure(6)
plot(x,p)
hold on
plot(x,(p_analytical),'o')
hold off
title('Pressure, numerical: line, analytical: circles')
grid on

rho_analytical(find(A == 1))
T_analytical(find(A == 1))
M_correct(find(A == 1))

rho_analytical(end)
T_analytical(end)
M_correct(end)


figure(1)
plot(dU1dt_avg_throat ,'b')
hold on
plot(dU2dt_avg_throat,'r')
plot(dU3dt_avg_throat,'g')
hold off
legend('\partialU1/\partialt', '\partialU2/\partialt', '\partialU3/\partialt', 'Location', 'best');
title('Residuals')
grid on
xlabel('Time step')














toc