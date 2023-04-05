% TRANSIENT ANALYSIS STEP

clear variables 

% data structure has : data(InputData,OutputData) and t

t1=load('lab3_order1_5.mat','t');
data1=load('lab3_order1_5.mat','data');
t2=load('lab3_order2_5.mat','t');
data2=load('lab3_order2_5.mat','data');

%1st order
u1_1=data1.data.InputData(1:100); % identification step
u1_2=data1.data.InputData(201:500); % verification
y1_1=data1.data.OutputData(1:100); % identification out
y1_2=data1.data.OutputData(201:500); % verification out
%2nd order
u2_1=data2.data.InputData(1:100); % identification step
u2_2=data2.data.InputData(201:500); % verification
y2_1=data2.data.OutputData(1:100); % identification out
y2_2=data2.data.OutputData(201:500); % verification out

%first order
t_1=t1.t(1:100);
t_2=t1.t(201:500);
%2nd order
t2_1=t2.t(1:100);
t2_2=t2.t(201:500);

%1st
N=length(y1_2);
%2nd
N2=length(y2_2);

%% First Order
figure
subplot(211)
plot(t_1,u1_1)
subplot(212)
plot(t_1,y1_1)

% yss is an avg of last values to eliminate noise
% Gain,K= (yss-y0)/(uss-u0)

yss=sum(y1_1(90:100))/11; %which is approx 2;
K=yss/4

% Read t0 the start time of the step, t1 the time where the output
% raises 0.632 of the difference. Compute T = t1 − t0 
% Computing T
   
y_T=0.632*yss; % is 1.26 closest is 1.24 which is at t=3.43
T=3.43

% The transfer function
% H = K/(Ts + 1)

H=tf(K,[T 1])

% lsim(sys,input signal,time)
y_approx=lsim(H,u1_2,t_2);

% MSE calculation
MSE=(1/N)*sum((y_approx-y1_2).^2)

figure
plot(t_2,y1_2,'b')
hold on
plot(t_2,y_approx,'r')
hold off
title("MSE = " +MSE)
xlabel('time'),ylabel('output')
legend('output','approx')

%% Second Order
figure
plot(t2_1,y2_1) % uss=0.5

yss2=sum(y2_1(90:100))/11; 

% Gain,K = yss/uss
K2=yss2/0.5

% Overshoot,M = (y(t1)-yss)/yss where y(t1) is peak value
M=(1.5-yss2)/yss2

% zeta = log(1/M) / sqrt( pi^2 + log(M)^2 )
zeta=log(1/M)/sqrt(pi^2+log(M)^2)

% T0=t3-t1 which is time for first 2 peaks
T0=2.36-0.78

% wn= 2pi / (T0 * sqrt( 1 - zeta^2 )
wn=2*pi/(T0*sqrt(1-zeta^2))

% H =( K * wn^2 )/( s^2 + 2*zeta*wn +wn^2 )
H2=tf(K2*(wn^2),[1 2*zeta*wn wn^2])

y2_approx=lsim(H2,u2_2,t2_2);

MSE2=(1/N2)*sum((y2_approx-y2_2).^2)

figure
plot(t2_2,y2_2)
hold on
plot(t2_2,y2_approx)
hold off

