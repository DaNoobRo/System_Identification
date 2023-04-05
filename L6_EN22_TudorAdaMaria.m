%% Data loading and variable declaration
clear variables

load('lab6_5.mat')
figure('Name','Plot ID data'),plot(id)
figure('Name','Plot VAL data'),plot(val)

% The length of the Data 
N_id=length(id.InputData);
N_val=length(val.InputData);

% The adjustable numbers of a and b parameters 
% for output and input regressors respectevly 
na=8;
nb=8;

% data assignment
y=id.OutputData; 
u=id.InputData;
y_val=val.OutputData;
u_val=val.InputData;

%% Creation of matrix Phi with N rows and na+nb columns
Phi=zeros(N_id,na+nb); 

for i=1:N_id
    % matrix part with output generated regressors
    for j=1:na
        if(i-j)<=0
            Phi(i,j)=0;
        else
            Phi(i,j)=(-1)*y(i-j);
        end
    end
    % matrix part with input generated regressors
    % i-(j-nb) is this an option that works for very different 
    % na and nb ? yes or no ?

    for j=(na+1):na+nb
        if(na+(i-j))<=0
            Phi(i,j)=0;
        else
            Phi(i,j)=u(na+(i-j));
        end
    end
end

Theta=Phi\y;    % Theta-> na+nb rows;

%% Prediction on validation

% PHI FOR VALIDATION DATA
Phi_val=zeros(N_val,na+nb);
% matrix part with output generated regressors
for i=1:N_val
    for j=1:na
        if(i-j)<=0
            Phi_val(i,j)=0;
        else
            Phi_val(i,j)=(-1)*y_val(i-j);
        end
    end
    % matrix part with input generated regressors
    for j=(na+1):na+nb
        if(na+(i-j))<=0
            Phi_val(i,j)=0;
        else
            Phi_val(i,j)=u_val(na+(i-j));
        end
    end
end

y_prediction=Phi_val*Theta; % PREDICTION
MSE=1/N_val*sum((y_prediction-y_val).^2); % MSE

%% PLOT prediction
figure('Name','Prediction on validation'),
plot(y_val),hold on
plot(y_prediction),hold off
title("MSE ="+MSE)

%% Simulation
y_approx=zeros(N_val,1);
holden=zeros(1,na+nb);
for i=2:N_val
    for j=1:na
       if(i-j)>0
            holden(1,j)=(-1)*y_approx(i-j);
       end
    end
    for l=(1+na):(na+nb)
        if(na+(i-j))>0
            holden(1,l)=u_val(na+(i-j));
        end
    end
   y_approx(i)=holden*Theta;
end

MSE_sim=1/N_val*sum((y_approx-y_val).^2); % MSE

%%
figure('Name','Simulation'),
plot(y_val,'b'), hold on
plot(y_approx,'r'),hold off
title("MSE ="+MSE_sim)
