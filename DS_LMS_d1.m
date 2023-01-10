function [e,w_hat,update_ratio,detection,false_alarm] = DS_LMS_d1(x,d,w_hat,P_up,var_noise,imp)

% Calculate threshold
v = 5;
mu = 1/(v*(length(w_hat)+1)*var(x));
alpha = (P_up/v)/(1-P_up/v);
tau = (1+alpha)*(qfuncinv(P_up/2))^2;

% Patameter
N = length(d);
L = length(w_hat);
w_hat_all = zeros(L,N);
y = zeros(1,N);
e = zeros(1,N);
x_vec = zeros(length(w_hat),1);

time = 0;
time_d = 0;
time_f = 0;

% Average Power Estimation Parameter
Nw = 11;
lambda = 0.9;
C1 = 1.483*(1+5/(Nw-1));
power_av = zeros(N,0);
f = zeros(Nw,1);

% DS_LMS_d
for i = 1:N
    x_vec = [x(i); x_vec(1:end-1)];
    y(i) = x_vec.'*w_hat;
    e(i) = d(i) - y(i);
    W = (x_vec.'*x_vec)^(-1);
    
    if i==1
        power_av(i) = (d(i)^2);
    else
        f = [(d(i)^2); f(1:end-1)];
        power_av(i) = lambda*power_av(i-1) + C1*(1-lambda)*median(f);
    end
    
    % Update Equation
    if abs(e(i))^2 < tau*var_noise
        w_hat = w_hat;
    else
        if (d(i)^2) >= (L+1)*power_av(i) %判斷有imp
            w_hat = w_hat;
            if imp(i) ~= 0 %實際有:detection
                time_d = time_d + 1;
            else %實際無:false alarm
                time_f = time_f + 1;
            end
        else %判斷沒有
            w_hat = w_hat + mu*x_vec*e(i)';
            time = time+1;
        end
    end
    w_hat_all(:,i) = w_hat;
end

w_hat = w_hat_all;
update_ratio = time/N;
detection = (time_d)/(sum(imp~=0));
false_alarm = (time_f)/(sum(imp==0));