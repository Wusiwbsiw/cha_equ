close all
clear all
clc


%16QAM
M = 16;
x = randi([0 M-1],40000,1);
tx = qammod(x,M);
%Channel
h = [0.05,0.1+0.15*sqrt(-1),0.9-0.15*sqrt(-1),0.15,-0.1+sqrt(-1)*0.2,0.1+sqrt(-1)*0.2,-0.01+0.15*sqrt(-1)];
x1=filter(h,1,tx); 
dB = 25;
n=randn(length(x1),1)+sqrt(-1)*randn(length(x1),1);   % additive white gaussian noise (complex)
n=n/norm(n)*10^(-dB/20)*norm(x);  % scale noise power
rx = x1+n;
%%                                  CMA 

mse_history_CMA_ = zeros(1,1);
ISI_history_CMA_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
 end
MSE = 1;
symbol_r = [];
mse_history_CMA= [];
ISI_history_CMA = [];
 
stability = 0;
phase_offset = 17;
phase_offset_rad = deg2rad(phase_offset);


mu = 3e-5;
for i = 1:39500%(length(x)-N+1)
    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_CMA = [mse_history_CMA,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    
    
	w =w - flip(mu*grad);
    

    
    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_CMA = [ISI_history_CMA ,ISI];  
    
    if i >11
        if(abs(ISI_history_CMA(i-9:i-1)-ISI_history_CMA(i-10:i-2))*ones(9,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end

SER = 1-account_0/length(judge_0);

mse_history_CMA_ = mse_history_CMA_+mse_history_CMA;
ISI_history_CMA_ = ISI_history_CMA_+ISI_history_CMA;
SER_ = SER_+SER;


end

mse_history_CMA_averaged = mse_history_CMA_/experiment_time;
ISI_history_CMA_averaged = ISI_history_CMA_/experiment_time;
SER_CMA_averaged = SER_/experiment_time;


disp('CMA done!')
%%                              CMA with momentum

mse_history_momentum_ = zeros(1,1);
ISI_history_momentum_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
  end
MSE = 1;
symbol_r = [];
mse_history_momentum= [];
ISI_history_momentum = [];

stability = 0;
phase_offset = 15;
phase_offset_rad = deg2rad(phase_offset);


mu = 2.5e-5;
mt = 0.3;
v = 0.2;
for i = 1:39500%(length(x)-N+1)
    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_momentum = [mse_history_momentum,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    
    v = mt*v + mu*grad;
    
	w =w - flip(v);
    

    
    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_momentum = [ISI_history_momentum ,ISI];  
    
    if i >11
        if(abs(ISI_history_momentum(i-9:i)-ISI_history_momentum(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_momentum_ = mse_history_momentum_+mse_history_momentum;
ISI_history_momentum_ = ISI_history_momentum_+ISI_history_momentum;
SER_ = SER_+SER;


end

mse_history_momentum_averaged = mse_history_momentum_/experiment_time;
ISI_history_momentum_averaged = ISI_history_momentum_/experiment_time;
SER_momentum_averaged = SER_/experiment_time;


disp('momentum done!')
%%                              CMA with NAG

mse_history_NAG_ = zeros(1,1);
ISI_history_NAG_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
end
MSE = 1;
symbol_r = [];
mse_history_NAG= [];
ISI_history_NAG = [];

stability = 0;
phase_offset = 20;
phase_offset_rad = deg2rad(phase_offset);


mu = 2.5e-5;
mt = 0.2;
v = 0.8;
for i = 1:39500%(length(x)-N+1)
    y = (flip(w-mt*v)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_NAG = [mse_history_NAG,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    
    v = mt*v + mu*grad;
    
	w =w - flip(v);
    

    
    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_NAG = [ISI_history_NAG ,ISI];  
    
    if i >11
        if(abs(ISI_history_NAG(i-9:i)-ISI_history_NAG(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_NAG_ = mse_history_NAG_+mse_history_NAG;
ISI_history_NAG_ = ISI_history_NAG_+ISI_history_NAG;
SER_ = SER_+SER;


end

mse_history_NAG_averaged = mse_history_NAG_/experiment_time;
ISI_history_NAG_averaged = ISI_history_NAG_/experiment_time;
SER_NAG_averaged = SER_/experiment_time;


disp('NAG done!')
%%                              CMA with Adagrad


mse_history_Adagrad_ = zeros(1,1);
ISI_history_Adagrad_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time

%initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
end
MSE = 1;
symbol_r = [];
mse_history_Adagrad= [];
ISI_history_Adagrad = [];

stability = 0;
phase_offset = 16;
phase_offset_rad = deg2rad(phase_offset);


G = zeros(11,1);
mu = 0.0136;
os=1e-8;
for i = 1:39500%(length(x)-N+1)
    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_Adagrad = [mse_history_Adagrad,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    G = G +grad.*grad;
    
	w =w - flip(  mu*grad./sqrt(G+os)  );

    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_Adagrad = [ISI_history_Adagrad ,ISI];  
    
    if i >11
        if(abs(ISI_history_Adagrad(i-9:i)-ISI_history_Adagrad(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_Adagrad_ = mse_history_Adagrad_+mse_history_Adagrad;
ISI_history_Adagrad_ = ISI_history_Adagrad_+ISI_history_Adagrad;
SER_ = SER_+SER;


end

mse_history_Adagrad_averaged = mse_history_Adagrad_/experiment_time;
ISI_history_Adagrad_averaged = ISI_history_Adagrad_/experiment_time;
SER_Adagrad_averaged = SER_/experiment_time;


disp('Adagrad done!')

%%                              CMA with Adadelta


mse_history_Adadelta_ = zeros(1,1);
ISI_history_Adadelta_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
end
MSE = 1;
symbol_r = [];
mse_history_Adadelta= [];
ISI_history_Adadelta = [];

stability = 0;
phase_offset = 16;
phase_offset_rad = deg2rad(phase_offset);


E_grad2 = zeros(11,1);
E_delta2 = zeros(11,1);

mt = 0.1;
mu = 1e-3;
ws = 30;
grad_history = zeros(11,ws);
os=1e-8;
for i = 1:39500%(length(x)-N+1)

    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_Adadelta = [mse_history_Adadelta,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    
    
	if  i <= ws                      
        grad_history(:,i) = grad;
        for k=1:i
            E_grad2 = mt*E_grad2+(1-mt)*grad_history(:,k).*grad_history(:,k); 
        end
	end
	if  i > ws                      
        grad_history = grad_history(:,2:end);
        grad_history(:,ws) = grad;
        E_grad2 = zeros(11,1);
        for k=1:ws
            E_grad2 = mt*E_grad2+(1-mt)*grad_history(:,k).*grad_history(:,k);
        end
	end

	w =w - flip( sqrt(E_delta2+os).*grad./sqrt(E_grad2)  );
    
	delta = -mu*grad./sqrt(E_grad2+os);
    E_delta2 = mt*E_delta2 + (1-mt)*delta.*delta;

    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_Adadelta = [ISI_history_Adadelta ,ISI];  
    
    if i >11
        if(abs(ISI_history_Adadelta(i-9:i)-ISI_history_Adadelta(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
	if i ==39500&&stability ==1
        disp('Nice')
	end
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_Adadelta_ = mse_history_Adadelta_+mse_history_Adadelta;
ISI_history_Adadelta_ = ISI_history_Adadelta_+ISI_history_Adadelta;
SER_ = SER_+SER;


end

mse_history_Adadelta_averaged = mse_history_Adadelta_/experiment_time;
ISI_history_Adadelta_averaged = ISI_history_Adadelta_/experiment_time;
SER_Adadelta_averaged = SER_/experiment_time;

disp('Adadelta done!')

%%                              CMA with RMSprop


mse_history_RMSprop_ = zeros(1,1);
ISI_history_RMSprop_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
end
MSE = 1;
symbol_r = [];
mse_history_RMSprop= [];
ISI_history_RMSprop = [];

stability = 0;
phase_offset = 14;
phase_offset_rad = deg2rad(phase_offset);


E_grad2 = zeros(11,1);

mu = 1e-3;
mt = 0.2 ;
os=1e-8;
for i = 1:39500%(length(x)-N+1)
    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_RMSprop = [mse_history_RMSprop,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    E_grad2 = mt*E_grad2+(1-mt)*grad.*grad;
	w =w - flip( mu*grad./sqrt(E_grad2+os));
    


    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_RMSprop = [ISI_history_RMSprop ,ISI];  
    
    if i >11
        if(abs(ISI_history_RMSprop(i-9:i)-ISI_history_RMSprop(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
	if i ==39500&&stability ==1
        disp('Nice')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_RMSprop_ = mse_history_RMSprop_+mse_history_RMSprop;
ISI_history_RMSprop_ = ISI_history_RMSprop_+ISI_history_RMSprop;
SER_ = SER_+SER;


end

mse_history_RMSprop_averaged = mse_history_RMSprop_/experiment_time;
ISI_history_RMSprop_averaged = ISI_history_RMSprop_/experiment_time;
SER_RMSprop_averaged = SER_/experiment_time;

disp('RMSprop done!')

%%                              CMA with Adam


mse_history_Adam_ = zeros(1,1);
ISI_history_Adam_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
end
MSE = 1;
symbol_r = [];
mse_history_Adam= [];
ISI_history_Adam = [];

stability = 0;
phase_offset = 14;
phase_offset_rad = deg2rad(phase_offset);

m = zeros(11,1);
v = zeros(11,1);
beta1 = 0.99;
beta2 = 0.86;
mu = 1e-3;
os=1e-8;

for i = 1:39500%(length(x)-N+1)
    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_Adam = [mse_history_Adam,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    
    
    m = beta1*m + (1-beta1)*grad; 
    v = beta2*v + (1-beta2)*grad.*grad;
    
    m_ = m/(1-beta1^i);
    v_ = v/(1-beta2^i);
    
    
	w =w - flip( mu*m_./(sqrt(v_)+os)  );
    

    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_Adam = [ISI_history_Adam ,ISI];  
    
    if i >11
        if(abs(ISI_history_Adam(i-9:i)-ISI_history_Adam(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_Adam_ = mse_history_Adam_+mse_history_Adam;
ISI_history_Adam_ = ISI_history_Adam_+ISI_history_Adam;
SER_ = SER_+SER;


end

mse_history_Adam_averaged = mse_history_Adam_/experiment_time;
ISI_history_Adam_averaged = ISI_history_Adam_/experiment_time;
SER_Adam_averaged = SER_/experiment_time;

disp('Adam done!')

%%                              CMA with Adamax

mse_history_Adamax_ = zeros(1,1);
ISI_history_Adamax_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
end
MSE = 1;
symbol_r = [];
mse_history_Adamax= [];
ISI_history_Adamax = [];

stability = 0;
phase_offset = 12;
phase_offset_rad = deg2rad(phase_offset);

m = zeros(11,1);
v = zeros(11,1);
beta1 = 0.9511;
beta2 = 0.0082;
mu = 0.001;
os=1e-8;

for i = 1:39500%(length(x)-N+1)
    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_Adamax = [mse_history_Adamax,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    
    
    m = beta1*m + (1-beta1)*grad; 
    m_ = m/(1-beta1^i);
    u = max(beta2*v,abs(grad));
	w =w - flip( mu*m_./u  );
    
	v = beta2*v + (1-beta2)*abs(grad).*abs(grad);
    

    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_Adamax = [ISI_history_Adamax ,ISI];  
    
    if i >11
        if(abs(ISI_history_Adamax(i-9:i)-ISI_history_Adamax(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_Adamax_ = mse_history_Adamax_+mse_history_Adamax;
ISI_history_Adamax_ = ISI_history_Adamax_+ISI_history_Adamax;
SER_ = SER_+SER;


end

mse_history_Adamax_averaged = mse_history_Adamax_/experiment_time;
ISI_history_Adamax_averaged = ISI_history_Adamax_/experiment_time;
SER_Adamax_averaged = SER_/experiment_time;


disp('Adamax done!')


%%                              CMA with Nadam

mse_history_Nadam_ = zeros(1,1);
ISI_history_Nadam_ = zeros(1,1);
SER_ = 0;
experiment_time = 100;

for p = 1:experiment_time
    %initial position
w = ones(11,1)*0;
w(5) = 1;   %w权重矩阵的初始延时为
w_init = w;
N=11;
Rx = [];
err = [];
for i = 1:39500%(length(x)-N+1)
    Rx = [Rx,rx(i:(i+N-1))];%考虑延时后的信号分为组
end
MSE = 1;
symbol_r = [];
mse_history_Nadam= [];
ISI_history_Nadam = [];

stability = 0;
phase_offset = 12;
phase_offset_rad = deg2rad(phase_offset);

m = zeros(11,1);
v = zeros(11,1);
beta1 = 0.00986;
beta2 = 0.0082;
mu = 0.00097;
os=1e-8;

for i = 1:39500%(length(x)-N+1)
    y = (flip(w)).'*Rx(:,i);%考虑均衡器每个抽头都有信号
    e_CMA = abs(y)^2-13.2;
    MSE = 0.99*MSE + (1-0.99)*(abs(e_CMA));
    mse_history_Nadam = [mse_history_Nadam,MSE];
    grad = e_CMA*y*conj(Rx(:,i)); 
    
    
    m = beta1*m + (1-beta1)*grad; 
    v = beta2*v + (1-beta2)*grad.*grad;
    
    m_ = m/(1-beta1^i);
    v_ = v/(1-beta2^i);
    
    
	w =w - flip( mu*(beta1*m_+(1-beta1)*grad/(1-beta1^i))./(os+sqrt(v_))   );
    

    abs_h2 = abs(conv(h,w)).*abs(conv(h,w));
    [a,index] = max(abs_h2);%信号冲激响应和均衡器系数卷积后最大值的下标为index，得延时为index-1
    ISI = (ones(1,length(abs_h2))*abs_h2-max(abs_h2))/max(abs_h2);
    ISI_history_Nadam = [ISI_history_Nadam ,ISI];  
    
    if i >11
        if(abs(ISI_history_Nadam(i-9:i)-ISI_history_Nadam(i-10:i-1))*ones(10,1))<1e-3&&stability == 0
            stability = 1;
            stability_number = i;
        end
    end
    
    if i ==39500&&stability ==0
        disp('DAMN IT')
    end
    
    if stability==1
        symbol_r = [symbol_r, qamdemod(y*exp(phase_offset_rad*1i),16)];
    end
    
end

symbol_t = x(stability_number+4:stability_number+3+length(symbol_r)).';
judge_0 = symbol_r-symbol_t;
account_0 = 0;

for i=1:length(judge_0)
    if judge_0(i)==0
        account_0 = account_0+1;
    end
end
SER = 1-account_0/length(judge_0);
mse_history_Nadam_ = mse_history_Nadam_+mse_history_Nadam;
ISI_history_Nadam_ = ISI_history_Nadam_+ISI_history_Nadam;
SER_ = SER_+SER;


end

mse_history_Nadam_averaged = mse_history_Nadam_/experiment_time;
ISI_history_Nadam_averaged = ISI_history_Nadam_/experiment_time;
SER_Nadam_averaged = SER_/experiment_time;

  disp('Nadam done!')

%%                              plot

leng = 39500;
it = 0:min(leng)-1;
figure(1);
%plot(it,10*log10(ISI_history_CMA),'o',it,10*log10(ISI_history_momentum),'+',it,10*log10(ISI_history_NAG),'*',it,10*log10(ISI_history_Adagrad),'.',it,10*log10(ISI_history_Adadelta),'square',it,10*log10(ISI_history_RMSprop),'diamond',it,10*log10(ISI_history_Adam),'^',it,10*log10(ISI_history_Adamax),'pentagram',it,10*log10(ISI_history_Nadam));
plot(it,10*log10(ISI_history_CMA_averaged),it,10*log10(ISI_history_momentum_averaged),it,10*log10(ISI_history_NAG_averaged),it,10*log10(ISI_history_Adagrad_averaged),it,10*log10(ISI_history_Adadelta_averaged))%,;

xlabel('iteration time');
ylabel('ISI history(dB)');
legend('CMA','momentum','NAG','Adagrad','Adadelta');
figure(2);
plot(it,10*log10(ISI_history_RMSprop_averaged),it,10*log10(ISI_history_Adam_averaged),it,10*log10(ISI_history_Adamax_averaged),it,10*log10(ISI_history_Nadam_averaged));
xlabel('iteration time');
ylabel('ISI history(dB)');
legend('RMSprop','Adam','Adamax','Nadam');

figure(3);
plot(it,10*log10(mse_history_CMA_averaged),it,10*log10(mse_history_momentum_averaged),it,10*log10(mse_history_NAG_averaged),it,10*log10(mse_history_Adagrad_averaged),it,10*log10(mse_history_Adadelta_averaged))%,;
xlabel('iteration time');
ylabel('mse history(dB)');
legend('CMA','momentum','NAG','Adagrad','Adadelta');

      
figure(4);
plot(it,10*log10(mse_history_RMSprop_averaged),it,10*log10(mse_history_Adam_averaged),it,10*log10(mse_history_Adamax_averaged),it,10*log10(mse_history_Nadam_averaged));
xlabel('iteration time');
ylabel('mse history(dB)');
legend('RMSprop','Adam','Adamax','Nadam');

