              function    [e,w_hat,update_ratio] = DS_LMS_Diniz(x,d,w_hat,P_up,var_noise)

% Calculate threshold
v = 5;
mu = 1/(v*(length(w_hat)+1)*var(x));
alpha = (P_up/v)/(1-P_up/v);
tau = (1+alpha)*(qfuncinv(P_up/2))^2;

% Patameterie
N = length(d);
L = length(w_hat);
w_hat_all = zeros(L,N);
y = zeros(1,N);
e = zeros(1,N);
x_vec = zeros(length(w_hat),1);
time = 0;

% DS_LMS
for i=1:N
    x_vec = [x(i);x_vec(1:end-1)];
    y(i) = x_vec.'*w_hat;
    e(i) = d(i) - y(i);
    
    % Update Equation
    if abs(e(i))^2 < tau*var_noise
        w_hat = w_hat;
    else
        w_hat = w_hat + mu*x_vec*e(i)';
        time = time + 1;
    end
    w_hat_all(:,i) = w_hat;
end

w_hat = w_hat_all;
update_ratio = time/N;