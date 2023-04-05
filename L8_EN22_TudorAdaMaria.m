%% Clear and variable loading
clear variables

load('lab8_5.mat')
figure('Name','Plot ID data'),plot(id)
figure('Name','Plot VAL data'),plot(val)

 N_id=length(id.OutputData);
 y=id.OutputData;
 u=id.InputData;
 Ts=id.Ts;

%% Adjustable variables

 lmax=249; % maximum iteration range in [100, 250]
 alpha=0.01; % stepsize in [0.01,0.5]
 
 Theta=ones(2,lmax); % Theta=[f;b]
 Theta(:,1)=0.491;

 threshold=1e-4; %threshold in [1e-4, 1e-1]

%% Initial parameters for the model

l=1; % initialization of indices
error=zeros(N_id,1); % initialization of error
dError=zeros(2,N_id); % derivative of error of Theta 
cond=1;

%% Cycle as presented in pseudocod
while ((l<lmax) && cond>threshold)
    
    % Recursion formulas
    for k=1:(N_id-1)
        error(k+1)=-Theta(1,l)*error(k)+y(k+1,1)+Theta(1,l)*y(k,1)-Theta(2,l)*u(k,1);
        dError(1,k+1)=-error(k,1)-Theta(1,l)*dError(1,k)+y(k,1);
        dError(2,k+1)=-Theta(1,l)*dError(2,k)-u(k,1);
    end
    % gradient of the objective function
    dV=[0;0];
    for i=1:N_id
        dV=dV+2/N_id*error(i,1)*dError(:,i);
    end

    H=0; % Hessian of the objective function
    TdError=(dError.');
    for i=1:N_id
        H=H+2/N_id*dError(:,i)*TdError(i,:);
    end

    Theta(:,l+1)=Theta(:,l)-alpha*H^(-1)*dV;
    
    cond=norm(Theta(:,l+1)-Theta(:,l));

    l=l+1; % incremented
end

%% Creating the model with idpoly
% idpoly(A,B,C,D,F,0,Ts)
% OE model -> y(t)=B(q)/F(q)*u(t−nk)+e(t)
B=[0 Theta(1,1:l)];
F=[1 Theta(2,1:l)];

model=idpoly(1,B,1,1,F,0,Ts);
figure('Name','Model Comparison')
compare(val,model)

