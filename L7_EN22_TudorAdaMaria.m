%% Loading variables and validation data declaration
clear variables
load uval.mat

index=5;
data=system_simulator(index,u);
figure('Name','generated validation data')
subplot(211),plot(u);
subplot(212),plot(data.OutputData)

%% Variables

P=zeros(10,1); % Periods ~ Persistent excitation ?
fit=zeros(10,1); % how well it fits vector, dependent on m
N=300; % sample number
%interval of PRBS
    a=0.5;
    b=1;
    %

%% Generating the identification input PRBS

for m=3:10  % between 3 and 10, this is the number of states : see course
            % for coeff a = 1 depending on m 

P(m)=2^m-1; % PRBS period

uid=PRBS(N,m,a,b);
data_id=system_simulator(index,uid);

figure('Name','generated identification data')
subplot(211),plot(uid);
subplot(212),plot(data_id.OutputData)

%% Curs code for arx model application

% taking different possible values for 
% input and output dependent parameters
na = 1:15;
nb = 1:15;
time_delay = 1:5;
Ns = struc(na, nb, time_delay);
Models = arxstruc(data_id, data, Ns); % creating arx model for each
Nbest = selstruc(Models, 0); % selecting model for the minimum MSE 

% plotting best model
figure('Name',"Comparison with m ="+m)
Best_fit = arx(data_id, Nbest);
compare(Best_fit, data);
%fit(m)=Best_fit.Report.Fit.MSE;
[ymod,fit(m),ic]=compare(Best_fit,data);
end
%% MSE depending on m

figure('Name','MSE depending on m')
plot(3:10,fit(3:10));

%% Function creation, at EOF

function [uid] = PRBS(N,m,a,b)
% Creates a PRBS signal of length N in interval [a,b]
%   N- length of signal; m-number of states; a,b-interval limits

    x=zeros(m,N); % N-> input signal length, m -> number of states 
    x(1,:)=zeros(1,N)+1;
    coeff=zeros(1,10); % line vector of length = 10 
    uid=zeros(N,1);

    % coefficients depending on m
    coeff(1)=int8(ismember(m,[3 4 6 7 8]));
    coeff(2)=int8(ismember(m,[5 8]));
    coeff(3)=int8(ismember(m,[3 10]));
    coeff(4)=int8(ismember(m,[4 9]));
    coeff(5)=int8(ismember(m,5));
    coeff(6)=int8(ismember(m,6));
    coeff(7)=int8(ismember(m,[7 8]));
    coeff(8)=int8(ismember(m,8));
    coeff(9)=int8(ismember(m,9));
    coeff(10)=int8(ismember(m,10));
    
    for j=1:(N-1)
        first=0;
        for i=1:m
            first=first+coeff(i)*x(i,j);
        end
        x(1,j+1)=mod(first,2); % modulo 2 operation
        for i=1:(m-1)
           x(i+1,j+1)=x(i,j);
        end
        uid(j)=x(m,j);
    end
    uid(N)=x(m,N);
    uid=a+(b-a)*uid;
end

