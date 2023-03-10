function NMSD = Normalized_Mean_Squaue_Deviation2(w,w_hat)

N = size(w_hat,2);
L = length(w);

for i = 1:N
    if i<N/2
        NMSD(i) = norm(w(:,1)-w_hat(:,i))^2/norm(w(:,1))^2;
    else
        NMSD(i) = norm(w(:,2)-w_hat(:,i))^2/norm(w(:,2))^2;
    end
end