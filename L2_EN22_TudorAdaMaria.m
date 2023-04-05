% loading the data and creating

clear variables

load('lab2_13.mat');

figure('Name','Original Function')
plot(id.X,id.Y,'b')

% with N and N_val how many variables from the datasheet we take
% can be modified
N=41;
N_val=71;

MSE_id=zeros(1,20);
MSE=zeros(1,20);

for k = 1 : 20 % staring loop for testink diff n's

% Linear Model (without error) : y(k) = Phi'(k)*Theta
% n is the degree of the approx polynomial 
% Creating matrix Phi
n=k; 

% y_hat = theta_1 + x*theta_2 + x^2*theta_3 + .....

Phi=zeros(N,n);
for i=1:N
    for j=1:n
        Phi(i,j)=id.X(i)^(j-1);
    end
end


% Calculating Theta

Theta=Phi\id.Y';

% calculating the approximate values and plotting them

y_hat=zeros(1,N);
for i=1:n
    y_hat=y_hat+Theta(i).*id.X.^(i-1);
end

for i=1:N
    MSE_id(k)=MSE_id(k) + (id.Y(i)-y_hat(i))^2;
end

MSE_id(k) = (1/N)*MSE_id(k); 

% Validation Data

y_val=zeros(1,N_val);
for i=1:n
    y_val=y_val+Theta(i).*val.X.^(i-1);
end

% Calculating the mean square error to see accuracy of model

for i=1:N_val
    MSE(k)=MSE(k) + (val.Y(i)-y_val(i))^2;
end

MSE(k) = (1/N_val)*MSE(k); 

end %  END BIG LOOP

% Plot MSE vs n

figure('Name','MSE vs n')
k=1:20;
plot(k,MSE(k))

% By observing MSE we see that the samallest value for the error 
% is at around n = 7 with MSE=40.64; This is the optimal 
% polynomial degree

% Plotting at the optimal polynomial degree

n=7; 
Phi=zeros(N,n);
for i=1:N
    for j=1:n
        Phi(i,j)=id.X(i)^(j-1);
    end
end

Theta=Phi\id.Y';

y_approx=zeros(1,N_val);
for i=1:n
    y_approx=y_approx+Theta(i).*val.X.^(i-1);
end

MSE_ideal = MSE(n); 

% plotting the variables

figure('Name','Validation Data : Comparison')
plot(val.X,val.Y,'b')
hold on
plot(val.X,y_approx,'r') % this adds the the prev. created figure for comp.
hold off
title("MSE="+ MSE_ideal)
xlabel('inputs'),ylabel('outputs')
legend('output','approximation')
