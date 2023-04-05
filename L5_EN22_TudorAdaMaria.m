
clear variables
data5=load('lab5_5.mat');
%%
id=data5.id;
val=data5.val;
tid=data5.tid;
tval=data5.tval;
figure('Name','Identification Data'),plot(id);
figure('Name','Validation Data'),plot(val);

%%
%number of samples
N=length(id.OutputData);

% allocation with detrend as we do NOT have zero-mean values
y=detrend(id.OutputData); %detrend "takes out" the linear trend of the data
u=detrend(id.InputData);

%% Covariance functions
% tau=lags; ryu and ru calculation

ryu=zeros(1,N);
ru=zeros(1,N);

for tau=1:1:N
    for k=1:1:(N-tau)
       ryu(tau)=ryu(tau)+1/N*(y(k+tau-1)*u(k));
       ru(tau)=ru(tau)+1/N*(u(k+tau-1)*u(k));
    end
end

figure('Name','Tau'),plot(1:N,ru)

%% Creating the matrix Ryu, Ru and linear reggresion to find H
T=N;
M=round(T/10);

Ryu=ryu(1:T)';
Ru=zeros(T,M);
for i=1:T
    for j=1:M
        Ru(i,j)=ru(abs(i-j)+1);
    end
end
H=Ru\Ryu;

%% Creating the model
y_hat=conv(u,H);
figure('Name','Test model on identification'),plot(tid,y,tid,y_hat(1:2500))


%% Validation

% allocation with detrend as we do NOT have zero-mean values
y_val=detrend(val.OutputData); %detrend "takes out" the linear trend of the data
u_val=detrend(val.InputData);

%%
N_val=250;
y_hat_val=conv(u_val,H)
MSE=(1/N)*sum((y_hat_val(1:N_val)-y_val).^2);
figure('Name','Test model on validation'),plot(tval,y_val,tval,y_hat_val(1:250))
title("MSE ="+MSE)

