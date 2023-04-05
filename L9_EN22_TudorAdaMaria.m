%% Clear && Load
clear variables 

load('lab9_5.mat');
figure('Name','plotting the given data')
subplot(211),plot(id);
subplot(212),plot(val)

% The length of the Data 
N_id=length(id.InputData);
N_val=length(val.InputData);

% The adjustable numbers of a and b parameters 
% for output and input regressors respectevly 
na=n;
nb=na;
% data assignment
y=id.OutputData; 
u=id.InputData;
y_val=val.OutputData;
u_val=val.InputData;

%% ARX values
% Theta we find on identification data
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
    for j=(na+1):na+nb
        if(na+(i-j))<=0
            Phi(i,j)=0;
        else
            Phi(i,j)=u(na+(i-j));
        end
    end
end
Theta=Phi\y;    % Theta-> na+nb rows;

% Phi on validation
Phi_val=zeros(N_val,na+nb);
for i=1:N_val
    for j=1:na
        if(i-j)<=0
            Phi_val(i,j)=0;
        else
            Phi_val(i,j)=(-1)*y_val(i-j);
        end
    end
    for j=(na+1):na+nb
        if(na+(i-j))<=0
            Phi_val(i,j)=0;
        else
            Phi_val(i,j)=u_val(na+(i-j));
        end
    end
end

%% PLOT prediction
y_prediction=Phi_val*Theta; % PREDICTION
MSE=1/N_val*sum((y_prediction-y_val).^2); % MSE
figure('Name','Prediction on ARX'),
plot(y_val),hold on
plot(y_prediction),hold off
title("MSE ="+MSE)

%% Simulated ARX values

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

%% Ploting simulated ARX
MSE_sim=1/N_val*sum((y_approx-y_val).^2); % MSE

figure('Name','Simulation ARX'),
plot(y_val,'b'), hold on
plot(y_approx,'r'),hold off
title("MSE ="+MSE_sim)

%% IV model

% Z(k)
Z=zeros(N_val,na+nb); 

for i=1:N_val
    % matrix part with output generated regressors
    for j=1:na
        if(i-j)<=0
            Z(i,j)=0;
        else
            Z(i,j)=(-1)*y_approx(i-j);
        end
    end
    for j=(na+1):na+nb
        if(na+(i-j))<=0
            Z(i,j)=0;
        else
            Z(i,j)=u_val(na+(i-j));
        end
    end
end

% PHIMI (na + nb) x (na + nb) matrix & YMI (na + nb) x 1 vector

PHIMI=0;
tPhi_val=(Phi_val.');
for i=1:N_val
    middleman=(Z(i,:).')*(tPhi_val(:,i).');
    PHIMI=PHIMI+1/N_val*middleman;
end
YMI=zeros(na+nb,1);
for i=1:N_val
    YMI=YMI+1/N_val*(Z(i,:).')*(y_val(i,1).');
end

Theta_pebune=PHIMI\YMI;

%% idpoly model
% idpoly(A,B,C,D,F,0,Ts)
A=[1 Theta_pebune(1:na).']; % "leading coeff pf A must be 1"
B=[0 Theta_pebune((na+1):(nb+na)).']; % " put zero from prev lab "
model=idpoly(A,B,[],[],[],0,val.Ts);

figure('Name','Model Comparison')
compare(val,model)



