%% Clear && Load
clear variables 

load('lab10_5.mat');
figure('Name','plotting the given data')
subplot(211),plot(id);
subplot(212),plot(val);

%% Variable declaration
na=3*n; % taken from the data
nb=na;

% initial values
Theta_init=zeros(1,na+nb);
invP_init=100*eye(na+nb);

% The length of the Data 
N_id=length(id.InputData);

tenP=round(N_id*10/100);

%% Function Call : generating Theta

phiRec = phiRecGen(id,na,nb,Theta_init,invP_init);

%% idpoly model 100
% idpoly(A,B,C,D,F,0,Ts)
A=[1 (phiRec(1:na,N_id).')]; % "leading coeff pf A must be 1"
B=[0 (phiRec((na+1):(na+nb),N_id).')]; % " put zero from prev lab "
model=idpoly(A,B,[],[],[],0,val.Ts);

figure('Name','Model Comparison 100% data')
compare(val,model)

%% idpoly model 10
% idpoly(A,B,C,D,F,0,Ts)
A=[1 (phiRec(1:na,tenP).')]; % "leading coeff pf A must be 1"
B=[0 (phiRec((na+1):(na+nb),tenP).')]; % " put zero from prev lab "
model=idpoly(A,B,[],[],[],0,val.Ts);

figure('Name','Model Comparison 10% data')
compare(val,model)

%% Matlab function rarx
mat_matrix=rarx(id,[na nb 1],'ff',1,Theta_init,invP_init);

% idpoly model 100
A=[1 mat_matrix(N_id,1:na)]; 
B=[0 mat_matrix(N_id,(na+1):(na+nb))];
mat_model=idpoly(A,B,[],[],[],0,val.Ts);

figure('Name',' matlab Model Comparison 100% data')
compare(val,mat_model)

% idpoly model 10
A=[1 mat_matrix(tenP,1:na)]; 
B=[0 mat_matrix(tenP,(na+1):(na+nb))];
mat_model=idpoly(A,B,[],[],[],0,val.Ts);

figure('Name',' matlab Model Comparison 10% data')
compare(val,mat_model)

%% Function
function [phiRec] = phiRecGen(id,na,nb,Theta_init,invP_init)

% RETRIEVE u(k),y(k) from identification
u=id.InputData;
y=id.OutputData;
N_id=length(id.InputData);

% DATA INITIALIZATION

% Phi from ARX
Phi=zeros(N_id,na+nb); 
% error vector , in code will work with each element, size doesn't matter
% ;)
error=zeros(N_id,1); 
% will store Theta(k) : (na x nb) column vector 
phiRec=zeros(na+nb,N_id); 
phiRec(:,1)=Theta_init; % initialization
% Inverse of matrix P : (na + nb) x (na x nb) 
invP_prev=invP_init;

% LOOP
for k=2:N_id
    
    % FORM ARX REGRESSOR VECTOR
    % matrix part with output generated regressors
    for j=1:na
        if(k-j)<=0
            Phi(k,j)=0;
        else
            Phi(k,j)=(-1)*y(k-j);
        end
    end
    % matrix part with input generated regressors
    for j=(na+1):na+nb
        if(na+(k-j))<=0
            Phi(k,j)=0;
        else
            Phi(k,j)=u(na+(k-j));
        end
    end

    % Calculating error : scalar = scalar - row(na x nb)*column(na x nb)   
    error(k)=y(k)-Phi(k,:)*phiRec(:,k-1); % aci nu era nevoie sa pastrez eroarea
    % Update Inverse : (na + nb) x (na + nb) matrix
    invP=invP_prev-(invP_prev*(Phi(k,:).')*Phi(k,:)*invP_prev)/...
        (1+Phi(k,:)*invP_prev*(Phi(k,:).'));
    % Compute Weight : (na+nb) column vector
    W=invP*(Phi(k,:).'); 
    % Update Theta
    phiRec(:,k)=phiRec(:,k-1)+W*error(k);
    % Update invP_prev
    invP_prev=invP;

% end k : 2 -> N_id loop
end

% EOFunction
end
