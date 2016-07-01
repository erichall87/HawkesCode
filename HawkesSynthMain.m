% Script to calculate rates and W matrix, which matches the notation
% from the paper.  Will not save estimates of rates or many estimates of W
% for computation efficiency
% Our algorithms are entitled DMD as they are based on our paper "Online
% Convex Optimization in Dynamic Environments," IEEE Journal of Selected
% Topics in Signal Processing - Signal Processing for Big Data,
% vol. 9, no 4. 2015"
% Replications section VII-B of "Tracking Dynamic Point Processes on
% Networks," IEEE Replicates Section VII-A of "Tracking Dynamic Point
%Processes on Networks" - IEEE Journal on Information Theory , 
%vol 62, no 7. 2016.

% Date created: 6/6/2014
% Last modified: 6/30/2016
% Written by Eric C. Hall

clear;clc; close all
graphing = true; %Graph the results
saving = false; %Save the results in a .mat file
savefigs = false; %Save figures as .pdf and .fig

delta=.1; %Time discretization level
T_horizon=1e3; %Length of sensing time
no_iters = 4; %Number of independent trials of experiment

N_delta=ceil(T_horizon/delta);
eta=10/sqrt(N_delta); %Step size associated with rate
rho=.01/sqrt(N_delta); %Step size associated with network
tau=1e-3; %Network regularization parameter

alpha=exp(-1); %Rate of decay of true influence function
beta=.9; %Rate of decay of mismatch influence function
p = 100; %Size of network
mu_bar = rand(p,1)*.009+.001; %Baseline firing rate
mu_0 = mu_bar;

if T_horizon/delta~=N_delta
    warning('Works better when the Time horizon is a multiple of delta')
end

save_iters=min(20000,floor(N_delta/10)); %Save one out of every 'save_iters' rates estimates
no_saves=floor(N_delta/save_iters);

%At various times compute the loss with respect to current estimate
%using full data
batch_loss_DMD=zeros(no_iters,no_saves+1);
batch_loss_OGD=zeros(no_iters,no_saves+1);
batch_loss_DMD_beta=zeros(no_iters,no_saves+1);
batch_loss_OGD_beta=zeros(no_iters,no_saves+1);
loss_cont_true = cell(no_iters,1);
%Instantaneous losses of various estimators
loss_DMD=zeros(no_iters,N_delta);
loss_DMD_W0=zeros(no_iters,N_delta);
loss_DMD_exp=zeros(no_iters,N_delta);
loss_OGD=zeros(no_iters,N_delta);

loss_DMD_beta=zeros(no_iters,N_delta);
loss_DMD_exp_beta=zeros(no_iters,N_delta);
loss_DMD_W0_beta=zeros(no_iters,N_delta);
loss_OGD_beta=zeros(no_iters,N_delta);

loss_true=zeros(no_iters,N_delta);
%Network estimates
W_final_DMD = zeros(p,p,no_iters);
W_final_DMD_beta = zeros(p,p,no_iters);
W_final_OGD = zeros(p,p,no_iters);
W_final_OGD_beta = zeros(p,p,no_iters);
W_true = zeros(p,p,no_iters);

for kk = 1:no_iters
    tic
    fprintf('Iteration number: %d\n',kk)
    %Generate true network
    W = zeros(p);
    for a = 0:p-1
        for b = 0:p-1
            if floor(a/20) == floor(b/20)
                W(a+1,b+1) = rand;
            elseif rand>.8
                W(a+1,b+1) = rand*.3;
            end
        end
    end
    W = W*.8/norm(W,2);
    W_true(:,:,kk) = W;
    %Generate data
    [T, T_index, mu_T] = HawkesMV(W,mu_bar,alpha,mu_0,T_horizon);
    N=length(T);
    
    lambda_DMD=mu_bar; %rate estimate for DMD with known W
    lambda_DMD_W0=mu_bar; %rate estimate for DMD with initial W
    lambda_DMD_exp=mu_bar; %rate estimate for DMD with learned W
    lambda_OGD=mu_bar; %rate estimate for OGD learned W
    
    %Same as above but with incorrect influence function
    lambda_DMD_beta=mu_bar;
    lambda_DMD_exp_beta=mu_bar;
    lambda_DMD_W0_beta=mu_bar;
    lambda_OGD_beta=mu_bar;
    
    lambda_true=mu_bar; %lambda calculated directly from data, known W and alpha
    
    %Network estimates
    W0=zeros(p);
    W_DMD=W0;
    W_DMD_beta=W0;
    W_OGD=W0;
    W_OGD_beta=W0;
    
    %Conversion factors
    K_t=zeros(p,1);
    C_t=zeros(p,1);
    K2_t=zeros(p,1);
    C2_t=zeros(p,1);
    
    %Calculating batch losses
    loss_cont_true{kk}=ContinuousLoss(T,T_index,alpha,mu_bar,W);
    batch_loss_DMD(kk,1)=sum(ContinuousLoss(T,T_index,alpha,mu_bar,W0));
    batch_loss_OGD(kk,1)=batch_loss_DMD(kk,1);
    batch_loss_DMD_beta(kk,1)=batch_loss_DMD(kk,1);
    batch_loss_OGD_beta(kk,1)=batch_loss_DMD(kk,1);
    
    %Our Algorithms
    fprintf('Running algorithms: ')
    J=1; %Counter variable that goes through events to place them into data vectors
    for n=1:N_delta
        pct = N_delta/10;
        if n == N_delta
            fprintf('100%%\n')
        elseif mod(n,pct)<mod(n-1,pct)
            fprintf('%.0f%%, ', n/N_delta*100)
        end
        
        if mod(n,save_iters)==0
            batch_loss_DMD(kk,(n/save_iters)+1)=sum(ContinuousLoss(T,T_index,alpha,mu_bar,W_DMD));
            batch_loss_OGD(kk,(n/save_iters)+1)=sum(ContinuousLoss(T,T_index,alpha,mu_bar,W_OGD));
            batch_loss_DMD_beta(kk,(n/save_iters)+1)=sum(ContinuousLoss(T,T_index,alpha,mu_bar,W_DMD_beta));
            batch_loss_OGD_beta(kk,(n/save_iters)+1)=sum(ContinuousLoss(T,T_index,alpha,mu_bar,W_OGD_beta));
        end
        %Compute data vectors
        x_t=zeros(p,1);
        y_t=zeros(p,1);
        y_t_beta=zeros(p,1);
        if J<=N
            while (J<=N && (ceil(T(J)/delta)<=n))
                x_t(T_index(J))=x_t(T_index(J))+1;
                y_t(T_index(J))=y_t(T_index(J))+alpha^((n+1)*delta-T(J));
                y_t_beta(T_index(J))=y_t_beta(T_index(J))+beta^((n+1)*delta-T(J));
                J=J+1;
            end
        end
        
        %Calculate Losses
        loss_DMD(kk,n)     =delta*sum(lambda_DMD)      -x_t'*log(lambda_DMD);
        loss_DMD_W0(kk,n)  =delta*sum(lambda_DMD_W0)   -x_t'*log(lambda_DMD_W0);
        loss_DMD_exp(kk,n) =delta*sum(lambda_DMD_exp)  -x_t'*log(lambda_DMD_exp);
        loss_OGD(kk,n)     =delta*sum(lambda_OGD)      -x_t'*log(lambda_OGD);
        
        loss_DMD_beta(kk,n)    =delta*sum(lambda_DMD_beta)     -x_t'*log(lambda_DMD_beta);
        loss_DMD_W0_beta(kk,n) =delta*sum(lambda_DMD_W0_beta)  -x_t'*log(lambda_DMD_W0_beta);
        loss_DMD_exp_beta(kk,n)=delta*sum(lambda_DMD_exp_beta) -x_t'*log(lambda_DMD_exp_beta);
        loss_OGD_beta(kk,n)    =delta*sum(lambda_OGD_beta)     -x_t'*log(lambda_OGD_beta);
        
        loss_true(kk,n)     =delta*sum(lambda_true)       -x_t'*log(lambda_true);
        
        %Update rates
        if n<N_delta;
            %Update estimate of W
            W_DMD=W_DMD - rho * (delta * ones(p,1)*K_t' - (x_t./lambda_DMD_exp)*K_t');
            W_DMD=W_DMD - rho*tau;
            W_DMD=W_DMD.*(W_DMD>=0);
            
            W_OGD=W_OGD - rho*(delta * ones(p,1)*C_t' - (x_t./lambda_OGD)*C_t');
            W_OGD=W_OGD - rho*tau;
            W_OGD=W_OGD.*(W_OGD>=0);
            
            W_DMD_beta=W_DMD_beta - rho*(delta*ones(p,1)*K2_t' - (x_t./lambda_DMD_exp_beta)*K2_t');
            W_DMD_beta=W_DMD_beta- rho*tau;
            W_DMD_beta=W_DMD_beta.*(W_DMD_beta>=0);
            
            W_OGD_beta=W_OGD_beta - rho*(delta*ones(p,1)*C2_t'- (x_t./lambda_OGD_beta)*C2_t');
            W_OGD_beta=W_OGD_beta-rho*tau;
            W_OGD_beta=W_OGD_beta.*(W_OGD_beta>=0);
            
            K_t=(1-eta)*alpha^delta*K_t + y_t;
            C_t=alpha^delta*C_t + y_t;
            
            K2_t=(1-eta)*beta^delta*K2_t + y_t_beta;
            C2_t=beta^delta*C2_t + y_t_beta;
            
            %Step in direction of data
            lambda_DMD_W0      =(1-eta)*lambda_DMD_W0     +eta*x_t;
            lambda_DMD_W0_beta =(1-eta)*lambda_DMD_W0_beta+eta*x_t;
            
            %Apply dynamics
            lambda_DMD_W0=(alpha^delta)*lambda_DMD_W0          + W0*y_t        + (1-alpha^delta)*mu_bar;
            lambda_DMD_W0_beta=(beta^delta)*lambda_DMD_W0_beta + W0*y_t_beta   + (1-beta^delta)*mu_bar;
            
            %Convert to different W matrices
            lambda_DMD=lambda_DMD_W0+(W-W0)*K_t;
            lambda_DMD_exp=lambda_DMD_W0+(W_DMD-W0)*K_t;
            lambda_OGD=mu_bar+W_OGD*C_t;
            lambda_true=mu_bar + W*C_t;
            
            lambda_DMD_beta=lambda_DMD_W0_beta+(W-W0)*K2_t;
            lambda_DMD_exp_beta=lambda_DMD_W0_beta+(W_DMD_beta-W0)*K2_t;
            lambda_OGD_beta=mu_bar+W_OGD_beta*C2_t;
            
        end
        W_final_DMD(:,:,kk) = W_DMD;
        W_final_DMD_beta(:,:,kk) = W_DMD_beta;
        W_final_OGD(:,:,kk) = W_OGD;
        W_final_OGD_beta(:,:,kk) = W_OGD_beta;
        
    end
    if saving
        save SynthResults
    end
    toc
end
%% Plotting the results
if graphing
    %Set up filter to do moving average over 'zeta' time points
    %(length zeta*delta time window)
    zeta=5000;
    filt=ones(1,zeta)/zeta;
    filt=repmat(fft(filt,N_delta),no_iters,1);
    
    %Calculate moving average losses
    ma_loss_DMD=ifft(fft(loss_DMD,[],2).*filt,[],2);
    ma_loss_DMD=ma_loss_DMD(:,zeta:end);
    ma_loss_DMD_W0=ifft(fft(loss_DMD_W0,[],2).*filt,[],2);
    ma_loss_DMD_W0=ma_loss_DMD_W0(:,zeta:end);
    ma_loss_DMD_exp=ifft(fft(loss_DMD_exp,[],2).*filt,[],2);
    ma_loss_DMD_exp=ma_loss_DMD_exp(:,zeta:end);
    ma_loss_OGD=ifft(fft(loss_OGD,[],2).*filt,[],2);
    ma_loss_OGD=ma_loss_OGD(:,zeta:end);
    
    ma_loss_true = ifft(fft(loss_true,[],2).*filt,[],2);
    ma_loss_true = ma_loss_true(:,zeta:end);
    
    ma_loss_DMD_beta=ifft(fft(loss_DMD_beta,[],2).*filt,[],2);
    ma_loss_DMD_beta=ma_loss_DMD_beta(:,zeta:end);
    ma_loss_DMD_W0_beta=ifft(fft(loss_DMD_W0_beta,[],2).*filt,[],2);
    ma_loss_DMD_W0_beta=ma_loss_DMD_W0_beta(:,zeta:end);
    ma_loss_DMD_exp_beta=ifft(fft(loss_DMD_exp_beta,[],2).*filt,[],2);
    ma_loss_DMD_exp_beta=ma_loss_DMD_exp_beta(:,zeta:end);
    ma_loss_OGD_beta=ifft(fft(loss_OGD_beta,[],2).*filt,[],2);
    ma_loss_OGD_beta=ma_loss_OGD_beta(:,zeta:end);
    
    J=min(floor(N_delta/100),1000);
    
    figure(1)
    plot(delta*(zeta:J:N_delta),mean(ma_loss_DMD_W0(:,1:J:end),1),...
        delta*(zeta:J:N_delta), mean(ma_loss_DMD(:,1:J:end),1),...
        delta*(zeta:J:N_delta),mean(ma_loss_DMD_exp(:,1:J:end),1),...
        delta*(zeta:J:N_delta),mean(ma_loss_OGD(:,1:J:end),1),'r--',...
        delta*(zeta:J:N_delta),mean(ma_loss_true(:,1:J:end),1),'k--','linewidth',3)
    curr_leg=legend('Algorithm 1:W=0','Alg 1:True W','Algorithm 2','OGD','True Rate','location','east');
    set(curr_leg,'fontsize',18);
    xlabel('Time','fontsize',18)
    ylabel('Loss','fontsize',18)
    title('Moving Average Loss','fontsize',16)
    hold off
    ax = gca;
    FormatFigures(ax, savefigs, 'MALossAlg2')
    
    figure(2)
    plot(delta*(zeta:J:N_delta),median(ma_loss_OGD(:,1:J:end) - ma_loss_DMD_exp(:,1:J:end),1),'b-',...
        delta*(zeta:J:N_delta),prctile(ma_loss_OGD(:,1:J:end) - ma_loss_DMD_exp(:,1:J:end),95,1),'k-',...
        delta*(zeta:J:N_delta),prctile(ma_loss_OGD(:,1:J:end) - ma_loss_DMD_exp(:,1:J:end),5,1),'k-','linewidth',3)
    xlabel('Time','fontsize',18)
    ylabel('Loss','fontsize',18)
    title('Moving Average Difference of OGD and Alg 2','fontsize',16)
    curr_leg = legend('Median','95th Percentile','5th Percentile','location','southeast');
    set(curr_leg,'fontsize',18);
    ax = gca;
    FormatFigures(ax, savefigs, 'MADiffAlg2')
    
    figure(3)
    plot(delta*(zeta:J:N_delta),mean(ma_loss_DMD_W0_beta(:,1:J:end),1),...
        delta*(zeta:J:N_delta),mean(ma_loss_DMD_beta(:,1:J:end),1),...
        delta*(zeta:J:N_delta),mean(ma_loss_DMD_exp_beta(:,1:J:end),1),...
        delta*(zeta:J:N_delta),mean(ma_loss_OGD_beta(:,1:J:end),1),'r--',...
        delta*(zeta:J:N_delta),mean(ma_loss_true(:,1:J:end)),'k--','linewidth',3)
    curr_leg = legend('Algorithm 1:W=0','Algorithm 1:True W','Algorithm 2','OGD','True Rate','location','East');
    set(curr_leg,'fontsize',18)
    xlabel('Time','fontsize',18)
    ylabel('Loss','fontsize',18)
    title('Moving Average Loss, Mismatch','fontsize',16)
    ax = gca;
    FormatFigures(ax, savefigs, 'MALossAlg2Mismatch')
    
    figure(4)
    plot(delta*(zeta:J:N_delta),median(ma_loss_OGD_beta(:,1:J:end) - ma_loss_DMD_exp_beta(:,1:J:end),1),'b-',...
        delta*(zeta:J:N_delta),prctile(ma_loss_OGD_beta(:,1:J:end) - ma_loss_DMD_exp_beta(:,1:J:end),95,1),'k-',...
        delta*(zeta:J:N_delta),prctile(ma_loss_OGD_beta(:,1:J:end) - ma_loss_DMD_exp_beta(:,1:J:end),5,1),'k-','linewidth',3)
    xlabel('Time','fontsize',18)
    ylabel('Loss','fontsize',18)
    title('Moving Average Difference of OGD and Alg 2, Model Mismatch','fontsize',16)
    curr_leg = legend('Median','95th Percentile','5th Percentile','location','southeast');
    set(curr_leg,'fontsize',18);
    ax = gca;
    FormatFigures(ax, savefigs, 'MADiffAlg2Mismatch')
    
    figure(5)
    imagesc(W,[0 1*max(W(:))])
    axis off; axis square;
    title('True W Matrix','fontsize',18)
    colormap hot; colorbar;
    ax = gca;
    FormatFigures(ax, savefigs, 'TrueW')
    
    figure(6)
    imagesc(W_DMD,[0 1*max(W_DMD(:))])
    axis off; axis square;
    title('Alg 2 Estimate of W Matrix','fontsize',16)
    colormap hot; colorbar;
    ax = gca;
    FormatFigures(ax, savefigs, 'Alg2Estimate')
    
    figure(7)
    imagesc(W_OGD,[0 1*max(W_DMD(:))])
    axis off; axis square;
    title('OGD Estimate of W Matrix','fontsize',16)
    colormap hot; colorbar;
    ax = gca;
    FormatFigures(ax, savefigs, 'OGDEstimate')
    
    figure(8)
    imagesc(W_DMD_beta,[0 1*max(W_DMD(:))])
    axis off; axis square;
    title('Alg 2 Estimate of W Matrix, mismatch','fontsize',16)
    colormap hot; colorbar;
    ax = gca;
    FormatFigures(ax, savefigs, 'Alg2EstimateMismatch')

    figure(9)
    imagesc(W_OGD_beta,[0 1*max(W_DMD(:))])
    axis off; axis square;
    title('OGD Estimate of W Matrix, mismatch','fontsize',16)
    colormap hot; colorbar;
    ax = gca;
    FormatFigures(ax, savefigs, 'OGDEstimateMismatch')
    
    figure(10)
    plot(delta*((0:no_saves)*save_iters),mean(batch_loss_DMD,1),'b-',...
        delta*((0:no_saves)*save_iters),mean(batch_loss_OGD,1),'r--',...
        delta*((0:no_saves)*save_iters),mean(batch_loss_DMD_beta,1),'k-',...
        delta*((0:no_saves)*save_iters),mean(batch_loss_OGD_beta,1),'c-',...
        'linewidth',3)
    hold on
    errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_DMD,1),std(batch_loss_DMD,[],1),'b.')
    errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_OGD,1),std(batch_loss_OGD,[],1),'r.')
    errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_DMD_beta,1),std(batch_loss_DMD_beta,[],1),'k.')
    errorbar(delta*((0:no_saves)*save_iters),mean(batch_loss_OGD_beta,1),std(batch_loss_OGD_beta,[],1),'c.')
    hold off
    xlabel('Time','fontsize',18)
    ylabel('Loss','fontsize',18)
    title('Batch Losses vs Time','fontsize',16)
    curr_leg = legend('Alg 2','OGD','Alg 2, Mismatch','OGD, Mismatch',0);
    set(curr_leg,'fontsize',18)
    ax = gca;
    FormatFigures(ax, savefigs, 'BatchLoss')

    no_points = 500;
    thresh = linspace(0,max([W_final_DMD(:);W_final_DMD_beta(:); W_final_OGD(:);W_final_OGD_beta(:)]),no_points);
    
    P_d_DMD=zeros(1,no_points);
    P_f_DMD=zeros(1,no_points);
    P_d_OGD = zeros(1,no_points);
    P_f_OGD = zeros(1,no_points);
    P_d_DMD_beta=zeros(1,no_points);
    P_f_DMD_beta=zeros(1,no_points);
    P_d_OGD_beta=zeros(1,no_points);
    P_f_OGD_beta=zeros(1,no_points);
    
    W_pos = W_true > 0;
    W_sum = sum(W_pos(:));
    W_neg = sum(1 - W_pos(:));
    for t=1:length(thresh);
        W_DMD_temp = W_final_DMD > thresh(t);
        W_DMD_beta_temp = W_final_DMD_beta > thresh(t);
        W_OGD_temp = W_final_OGD > thresh(t);
        W_OGD_beta_temp = W_final_OGD_beta > thresh(t);
        
        P_d_DMD(t)=sum(W_DMD_temp(:).*W_pos(:))/W_sum;
        P_f_DMD(t)=sum(W_DMD_temp(:).*(1-W_pos(:)))/(W_neg);
        
        P_d_DMD_beta(t)=sum(W_DMD_beta_temp(:).*W_pos(:))/W_sum;
        P_f_DMD_beta(t)=sum(W_DMD_beta_temp(:).*(1-W_pos(:)))/(W_neg);
        
        P_d_OGD(t)=sum(W_OGD_temp(:).*W_pos(:))/W_sum;
        P_f_OGD(t)=sum(W_OGD_temp(:).*(1-W_pos(:)))/(W_neg);
        
        P_d_OGD_beta(t)=sum(W_OGD_beta_temp(:).*W_pos(:))/W_sum;
        P_f_OGD_beta(t)=sum(W_OGD_beta_temp(:).*(1-W_pos(:)))/(W_neg);
        
    end
    
    figure(11)
    plot(P_f_DMD,P_d_DMD,'b-',P_f_OGD,P_d_OGD,'r--',...
        P_f_DMD_beta,P_d_DMD_beta,'k-',...
        P_f_OGD_beta,P_d_OGD_beta,'c-','linewidth',3)
    xlabel('P_f','fontsize',18)
    ylabel('P_d','fontsize',18)
    title('ROC for W','fontsize',16)
    curr_leg=legend('Alg 2','OGD','Alg 2, Mismatch','OGD, Mismatch','location','southeast');
    set(curr_leg,'fontsize',18);
    axis square
    ax = gca;
    FormatFigures(ax, savefigs, 'ROC')
    
    W_true_top = zeros(size(W_true));
    for k = 1:no_iters
        if no_iters > 1
            vals = sort(reshape(W_true(:,:,k),[p^2,1]));
            W_true_top(:,:,k) = W_true(:,:,k) > vals(ceil(.9*p^2));
        else
            vals = sort(W(:));
            W_true_top = W_true > vals(ceil(.9*p^2));
        end
    end
    
    no_points = 500;
    thresh = linspace(0,max([W_final_DMD(:);W_final_DMD_beta(:); W_final_OGD(:);W_final_OGD_beta(:)]),no_points);
    
    P_d_DMD=zeros(1,no_points);
    P_f_DMD=zeros(1,no_points);
    P_d_OGD = zeros(1,no_points);
    P_f_OGD = zeros(1,no_points);
    P_d_DMD_beta=zeros(1,no_points);
    P_f_DMD_beta=zeros(1,no_points);
    P_d_OGD_beta=zeros(1,no_points);
    P_f_OGD_beta=zeros(1,no_points);
    
    
    W_sum = sum(W_true_top(:));
    W_neg = sum(1 - W_true_top(:));
    for t=1:length(thresh);
        W_DMD_temp = W_final_DMD > thresh(t);
        W_DMD_beta_temp = W_final_DMD_beta > thresh(t);
        W_OGD_temp = W_final_OGD > thresh(t);
        W_OGD_beta_temp = W_final_OGD_beta > thresh(t);
        
        P_d_DMD(t)=sum(W_DMD_temp(:).*W_true_top(:))/W_sum;
        P_f_DMD(t)=sum(W_DMD_temp(:).*(1-W_true_top(:)))/(W_neg);
        
        P_d_DMD_beta(t)=sum(W_DMD_beta_temp(:).*W_true_top(:))/W_sum;
        P_f_DMD_beta(t)=sum(W_DMD_beta_temp(:).*(1-W_true_top(:)))/(W_neg);
        
        P_d_OGD(t)=sum(W_OGD_temp(:).*W_true_top(:))/W_sum;
        P_f_OGD(t)=sum(W_OGD_temp(:).*(1-W_true_top(:)))/(W_neg);
        
        P_d_OGD_beta(t)=sum(W_OGD_beta_temp(:).*W_true_top(:))/W_sum;
        P_f_OGD_beta(t)=sum(W_OGD_beta_temp(:).*(1-W_true_top(:)))/(W_neg);
        
    end
    
    figure(12)
    plot(P_f_DMD,P_d_DMD,'b-',P_f_OGD,P_d_OGD,'r--',...
        P_f_DMD_beta,P_d_DMD_beta,'k-',...
        P_f_OGD_beta,P_d_OGD_beta,'c-','linewidth',3)
    xlabel('P_f','fontsize',18)
    ylabel('P_d','fontsize',18)
    title('ROC for Largest Elements of W','fontsize',16)
    curr_leg = legend('Alg 2','OGD','Alg 2, Mismatch','OGD, Mismatch','location','southeast');
    set(curr_leg,'fontsize',18)
    axis square
    ax = gca;
    FormatFigures(ax, savefigs, 'ROCTop')
end

