clear variables 

%initialization
first=load('lab4_order1_5.mat');
second=load('lab4_order2_5.mat');

%allocations for tuning on other data
ID_init=31;
StepsNr=100;
SignalNr=3;
ID_fin=ID_init+StepsNr-1;
VAL=ID_init+StepsNr*SignalNr-1;

%1st order allocations
t1=first.t;
FirstImp_id=first.data.InputData(ID_init:ID_fin);
FirstImp_val=first.data.InputData(ID_fin+1:VAL);
FirstOut_id=first.data.OutputData(ID_init:ID_fin);
FirstOut_val=first.data.OutputData(ID_fin+1:VAL);

N1=length(FirstOut_val);

%2nd order allocations
t2=second.t;
SecondImp_id=second.data.InputData(ID_init:ID_fin);
SecondImp_val=second.data.InputData(ID_fin+1:VAL);
SecondOut_id=second.data.OutputData(ID_init:ID_fin);
SecondOut_val=second.data.OutputData(ID_fin+1:VAL);

N2=length(SecondOut_val);

%% First Order

figure('Name','First Order Identfication')
subplot(121)
plot(t1(ID_init:ID_fin),FirstImp_id)
subplot(122)
plot(t1(ID_init:ID_fin),FirstOut_id)

iss=0.5; % read from the graph, uss=u0 in course

yss=sum(FirstOut_id(90:100))/11;  

K=yss/iss;

ymax=0.24; % approx from graph

% Read ymax and read the time constant T at the moment where
% the output decreases to 0.368 of the difference ymax − y0 

y_T=yss+0.368*(ymax-yss);

% => This from graph
t_2=1.8; % where y=y_T
t_1=1.2; % initial t_1
T=t_2-t_1;

% The Model

    % Matrix values 
    A=(-1)*(1/T);
    B=K/T;
    C=1;
    D=0;
    
    % the initial conditions
    x0=yss;

Hss=ss(A,B,C,D);

% lsim(sys,input,time,x0)
y_approx=lsim(Hss,FirstImp_val,t1(ID_fin+1:VAL),x0);

MSE=(1/N1)*sum((y_approx-FirstOut_val).^2);

figure('Name','Validation output and model'),
plot(t1(ID_fin+1:VAL),FirstOut_val,'b')
hold on
plot(t1(ID_fin+1:VAL),y_approx,'r')
title("MSE = "+MSE)

%% Second Order Modelling

figure('Name','Second Order Identfication')
subplot(121)
plot(t2(ID_init:ID_fin),SecondImp_id)
subplot(122)
plot(t2(ID_init:ID_fin),SecondOut_id)

iss_2=0.5; % read from the graph, uss=u0 in course

yss_2=sum(SecondOut_id(90:100))/11;  
K_2=yss_2/iss_2;

% Areas

t00=0.96;
t01=1.69;
t02=2.43;
Ts=t2(2)-t2(1);
r00=round(t00/Ts);
r01=round(t01/Ts);
r02=round(t02/Ts);

Aplus=Ts*sum(SecondOut_id(r00:r01)-yss_2);
Aminus=Ts*sum(yss_2-SecondOut_id(r01:r02));

M=Aminus/Aplus;
zeta=log(1/M)/sqrt(pi^2+log(M)^2);
t3=2.65;
t11=1.21;
T0=t3-t11;
wn=2*pi/(T0*sqrt(1-zeta^2));

%model

A=[0 1;-wn^2 -2*zeta*wn];
B=[0;K_2*wn^2];
C=[1 0];
D=0;
x0_2=[yss_2 0];

Hss_2=ss(A,B,C,D);
y_approx2=lsim(Hss_2,SecondImp_val,t2(ID_fin+1:VAL),x0_2);

MSE=(1/N2)*sum((y_approx2-SecondOut_val).^2);

figure('Name','Validation output and model'),
plot(t2(ID_fin+1:VAL),SecondOut_val,'b')
hold on
plot(t2(ID_fin+1:VAL),y_approx2,'r')
title("MSE = "+MSE)
