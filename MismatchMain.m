%Script to test the difference between our algorithm with known network
%versus directly plugging in all the known parameters into the definition
%of Hawkes process, under both perfect knowledge of generative model and
%model mistmatch. Replicates Section VII-A of "Tracking Dynamic Point
%Processes on Networks" - IEEE Journal on Information Theory ,
%vol 62, no 7. 2016

clear
clc
close all

saving = false;%Decide to store results
graphing = true; %Generate plots

%% Calculate rates at discrete intervals
no_iters = 4; %Number of independent trials to perform
alpha = exp(-1); %True influence function
T_horizon = 20000; %Length of sensing time
p = 2; %Size of network
mu_bar = ones(p,1)*.005; %Baseline firing rate
W = eye(p)*.75; %True network

delta=.1; %Time discretization parameter
N_delta=ceil(T_horizon/delta); %Number of time bins

eta=10/sqrt(N_delta); %Step size of Algorithm 1
beta=2*alpha; % Incorrect decay rate for model mismatch 1
B=5; %Length of square window for mismatch 2

%Matrices to store instantaneous losses
loss_true_disc=zeros(no_iters,N_delta);
loss_exp=zeros(no_iters,N_delta);
loss_exp_DMD=zeros(no_iters,N_delta);
loss_square=zeros(no_iters,N_delta);
loss_square_DMD=zeros(no_iters,N_delta);

tau_k = cell(p,1);
for iter = 1:no_iters;
    fprintf('Iteration number %2d, ',iter)
    
    %Generate data
    [T, T_index] = HawkesMV(W,mu_bar,alpha,mu_bar,T_horizon);
    N = length(T);
    
    %Vectors to store instantaneous estimates of rates
    lambda_true=zeros(p,N_delta);
    lambda_true(:,1)=mu_bar;
    
    lambda_exp=zeros(p,N_delta);
    lambda_exp(:,1)=mu_bar;
    
    lambda_exp_DMD=zeros(p,N_delta);
    lambda_exp_DMD(:,1)=mu_bar;
    
    lambda_square=zeros(p,N_delta);
    lambda_square(:,1)=mu_bar;
    
    lambda_square_DMD=zeros(p,N_delta);
    lambda_square_DMD(:,1)=mu_bar;
    
    %Convert continuous time data to count data in bins
    y_t=zeros(p,N_delta);
    x_t=zeros(p,N_delta);
    for n=1:N
        idx=ceil(T(n)/delta);
        y_t(T_index(n),idx)=y_t(T_index(n),idx)+alpha^(delta*(idx+1)-T(n));
        x_t(T_index(n),idx)=x_t(T_index(n),idx)+1;
    end
    
    for a = 1:p
        tau_k{a}=T(T_index==a);
    end
    for n=1:N_delta
        pct = N_delta/10;
        if n == N_delta
            fprintf('100%%\n')
        elseif mod(n,pct)<mod(n-1,pct)
            fprintf('%.0f%%, ', n/N_delta*100)
        end
        
        %Calculate Losses
        loss_true_disc(iter,n)   =delta*sum(lambda_true(:,n))        -x_t(:,n)'*log(lambda_true(:,n));
        loss_exp(iter,n)         =delta*sum(lambda_exp(:,n))         -x_t(:,n)'*log(lambda_exp(:,n));
        loss_exp_DMD(iter,n)     =delta*sum(lambda_exp_DMD(:,n))     -x_t(:,n)'*log(lambda_exp_DMD(:,n));
        loss_square(iter,n)      =delta*sum(lambda_square(:,n))      -x_t(:,n)'*log(lambda_square(:,n));
        loss_square_DMD(iter,n)  =delta*sum(lambda_square_DMD(:,n))  -x_t(:,n)'*log(lambda_square_DMD(:,n));
        
        %Update rates
        if n<N_delta;
            lambda_true(:,n+1)  =alpha^delta*lambda_true(:,n)   +W*y_t(:,n)+(1-alpha^delta)*mu_bar;
            lambda_exp(:,n+1)   =beta^delta*lambda_exp(:,n)     +W*y_t(:,n)+(1-beta^delta)*mu_bar;
            lambda_exp_DMD(:,n+1)=(1-eta)*lambda_exp_DMD(:,n)+eta*x_t(:,n);
            lambda_exp_DMD(:,n+1)=beta^delta*lambda_exp_DMD(:,n+1)+W*y_t(:,n)+(1-beta^delta)*mu_bar;
            lambda_square_DMD(:,n+1)=(1-eta)*lambda_square_DMD(:,n)+eta*x_t(:,n);
            
            A_t=zeros(p,1);
            for a=1:p
                t_a=tau_k{a};
                if sum(ceil(t_a/delta)*delta<delta*(n) & (delta*(n)-t_a) < B)==0
                    A_t(a)=.99;
                else
                    numer=sum((ceil(t_a/delta)<(n)) & (delta*(n+1)-t_a<=(B))) ;
                    denom=sum(ceil(t_a/delta)*delta<delta*(n) & delta*(n)-t_a < B);
                    A_t(a)=numer/denom;
                end
                
            end
            
            lambda_square(:,n+1)=A_t.*lambda_square(:,n)+W*y_t(:,n)+(1-A_t).*mu_bar;
            lambda_square_DMD(:,n+1)=A_t.*lambda_square_DMD(:,n+1)+W*y_t(:,n)+(1-A_t).*mu_bar;
            
        end
        
    end
end

if saving;
    save ModelMismatchResults
end

%% Plotting the results
if graphing
    cs_true = cumsum(loss_true_disc,2);
    cs_exp = cumsum(loss_exp,2);
    cs_square = cumsum(loss_square,2);
    cs_exp_DMD = cumsum(loss_exp_DMD,2);
    cs_square_DMD = cumsum(loss_square_DMD,2);
    
    avg_true = cs_true./repmat((1:N_delta)*delta,no_iters,1);
    avg_exp = cs_exp./repmat((1:N_delta)*delta,no_iters,1);
    avg_square = cs_square./repmat((1:N_delta)*delta,no_iters,1);
    avg_exp_DMD = cs_exp_DMD./repmat((1:N_delta)*delta,no_iters,1);
    avg_square_DMD = cs_exp_DMD./repmat((1:N_delta)*delta,no_iters,1);
    
    J=min(floor(N_delta/25),10000);
    
    figure(1)
    plot((1:N_delta)*delta,mean(cs_exp,1),'b-',...
        (1:N_delta)*delta,mean(cs_exp_DMD,1),'r-',...
        (1:N_delta)*delta,mean(cs_true,1),'k-','linewidth',3)
    legend('Direct Calculation with Incorrect Paramters',...
        'Algorithm 1 With Incorrect Paramters',...
        'Direct Calculation with Correct Parameters','Location','Northeast')
    xlabel('Time','fontsize',14)
    ylabel('Loss','fontsize',14)
    title(sprintf('Cumulative Loss, %d Iterations',no_iters),'fontsize',16)
    hold on
    errorbar((J:J:N_delta)*delta, mean(cs_exp(:,J:J:end),1),...
        std(cs_exp(:,J:J:end),[],1),'b.')
    errorbar((J:J:N_delta)*delta, mean(cs_exp_DMD(:,J:J:end),1),...
        std(cs_exp_DMD(:,J:J:end),[],1),'r.')
    errorbar((J:J:N_delta)*delta, mean(cs_true(:,J:J:end),1),...
        std(cs_true(:,J:J:end),[],1),'k.')
    hold off
    
    figure(2)
    plot((1:N_delta)*delta,mean(cs_square,1),'b-',...
        (1:N_delta)*delta,mean(cs_square_DMD,1),'r-',...
        (1:N_delta)*delta,mean(cs_true,1),'k-','linewidth',3)
    legend('Direct Calculation with Incorrect Paramters',...
        'Algorithm 1 With Incorrect Paramters',...
        'Direct Calculation with Correct Parameters','Location','Northeast')
    xlabel('Time','fontsize',14)
    ylabel('Loss','fontsize',14)
    title(sprintf('Cumulative Loss (rect), %d Iterations',no_iters),'fontsize',16)
    hold on
    errorbar((J:J:N_delta)*delta, mean(cs_square(:,J:J:end),1),...
        std(cs_square(:,J:J:end),[],1),'b.')
    errorbar((J:J:N_delta)*delta, mean(cs_square_DMD(:,J:J:end),1),...
        std(cs_square_DMD(:,J:J:end),[],1),'r.')
    errorbar((J:J:N_delta)*delta, mean(cs_true(:,J:J:end),1),...
        std(cs_true(:,J:J:end),[],1),'k.')
    hold off
    
    zeta=5000;
    filt=ones(1,zeta)/zeta;
    filt=repmat(fft(filt,N_delta),no_iters,1);
    ma_exp=ifft(fft(loss_exp,[],2).*filt,[],2);
    ma_exp_DMD=ifft(fft(loss_exp_DMD,[],2).*filt,[],2);
    ma_true=ifft(fft(loss_true_disc,[],2).*filt,[],2);
    ma_square=ifft(fft(loss_square,[],2).*filt,[],2);
    ma_square_DMD=ifft(fft(loss_square_DMD,[],2).*filt,[],2);
    
    figure(3)
    plot(delta*(zeta:N_delta),mean(ma_exp(:,zeta:end),1),'b-',...
        delta*(zeta:N_delta),mean(ma_exp_DMD(:,zeta:end),1),'r-',...
        delta*(zeta:N_delta),mean(ma_true(:,zeta:end),1),'k-','linewidth',3)
    legend('Direct Calculation with Incorrect Paramters',...
        'Algorithm 1 With Incorrect Paramters',...
        'Direct Calculation with Correct Parameters',0)
    xlabel('Time','fontsize',14)
    ylabel('Loss','fontsize',14)
    title(sprintf('Moving average Losses, %d Iterations', no_iters),'fontsize',16)
     
    figure(4)
    plot(delta*(zeta:N_delta),mean(ma_square(:,zeta:end),1),'b-',...
        delta*(zeta:N_delta),mean(ma_square_DMD(:,zeta:end),1),'r-',...
        delta*(zeta:N_delta),mean(ma_true(:,zeta:end),1),'k-','linewidth',3)
    legend('Direct Calculation with Incorrect Paramters',...
        'Algorithm 1 With Incorrect Paramters',...
        'Direct Calculation with Correct Parameters',0)
    xlabel('Time','fontsize',14)
    ylabel('Loss','fontsize',14)
    title(sprintf('Moving average Losses (rect), %d Iterations',no_iters),'fontsize',16)
    
    figure(5)
    excess_exp = ma_exp - ma_exp_DMD;
    plot(delta*(zeta:N_delta),median(excess_exp(:,zeta:N_delta),1),'b-','linewidth',3)
    hold on
    plot(delta*(zeta:N_delta),prctile(excess_exp(:,zeta:N_delta),95,1),'Color',[.6, .6, .6],'linewidth',3)
    plot(delta*(zeta:N_delta),prctile(excess_exp(:,zeta:N_delta),5,1),'k-','linewidth',3)
    xlabel('Time','fontsize',14)
    ylabel('Loss','fontsize',14)
    title(sprintf('Moving Average Loss, OGD - DMD, %d Iterations',no_iters),'fontsize',16)
    legend('Median','95th Percentile','5th Percentile',0)
    
    figure(6)
    excess_square = ma_square - ma_square_DMD;
    plot(delta*(zeta:N_delta),median(excess_square(:,zeta:N_delta),1),'b-','linewidth',3)
    hold on
    plot(delta*(zeta:N_delta),prctile(excess_square(:,zeta:N_delta),95,1),'Color',[.6, .6, .6],'linewidth',3)
    plot(delta*(zeta:N_delta),prctile(excess_square(:,zeta:N_delta),5,1),'k-','linewidth',3)
    xlabel('Time','fontsize',14)
    ylabel('Loss','fontsize',14)
    title(sprintf('Moving Average Loss (rect), OGD - DMD, %d Iterations',no_iters),'fontsize',16)
    legend('Median','95th Percentile','5th Percentile',0)
end