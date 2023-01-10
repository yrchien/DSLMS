% This code requires Statistics and Machine Learning Toolbox.
% You can change the setting of INPUT and IN_MODE to see the resulting
% performance (NMSD and averaged updating ratio)
% The detailed information can be found in the following paper
% The codes were implemented by Chih-Hsiang YU
% https://www.jstage.jst.go.jp/article/transfun/E105.A/2/E105.A_2021EAL2046/_article
% Ying-Ren CHIEN, Chih-Hsiang YU, Impulse-Noise-Tolerant Data-Selective LMS Algorithm, 
% IEICE Transactions on Fundamentals of Electronics, Communications and Computer Sciences, 
% 2022, vol. E105.A, no. 2 , pp. 114-117, 2022/02/01
% https://doi.org/10.1587/transfun.2021EAL2046, 
% https://www.jstage.jst.go.jp/article/transfun/E105.A/2/E105.A_2021EAL2046/_article/-char/ja
% Any comments are welcomed.

clear all;
close all;
clc;
ensamble = 100;

INPUT=1; %1 white Gaussian 2==>AR1 3==>AR2 4==>AR4
IN_MODE=1; %1 weak 2 mild 3 strong

for k = 1:ensamble
    % White Noise
    N = 1e4;
    n = 1:N;
    x = randn(1,N);
    
    switch(INPUT)
        case 1
            disp(['WHITE input:',num2str(k)]);
        case 2
            disp(['AR 1 input:',num2str(k)]);
            for q=2:N
                x(q)=x(q)+0.88*x(q-1);
            end
        case 3
            disp(['AR 2input:',num2str(k)]);
            x(2)=-0.55*x(1)+x(2);          
            for q=3:N
                x(q)=-0.55*x(q-1)-0.221*x(q-2)+x(q);       
            end
        case 4
            disp(['AR 4 input:',num2str(k)]);
            x(2)=-0.55*x(1)+x(2);
            x(3)=-0.55*x(2)-0.221*x(1)+x(3);
            x(4)=-0.55*x(3)-0.221*x(2)-0.49955*x(1)+x(4);
            for q=5:N
                x(q)=-0.55*x(q-1)-0.221*x(q-2)-0.49955*x(q-3)-0.4536*x(q-4)+x(q);
            end
        otherwise
            disp('ERROR INPUT');
            exit();
    end
    
    w=[0.1010 0.3030 0 -0.2020 -0.4040 -0.7071 -0.4040 -0.2020].';
    %   w = [0.3, 0.7, -0.5, -0.09, 0.19, 0.19, -0.28, -0.32, 0.51, 0.71, -0.77, -0.165, 0.86, 0.39, -0.32, -0.46]';
    %    w = randn(64,1);
    w = [w,-w];
    %w = [w,randn(length(w),1)];
    d = zeros(1,N);
    x_vec = zeros(length(w),1);
    for i = 1:N
        x_vec = [x(i); x_vec(1:end-1)];
        if i < N/2
            d(i) = x_vec.'*w(:,1);
        else
            d(i) = x_vec.'*w(:,2);
        end
    end
    
    % Background Noise
    var_d = var(d);
    SNRdB = 30;
    var_v = var_d./(10^(SNRdB/10));
    v = sqrt(var_v).*randn(1,N);
    dn1 = d + v;
    
  
    % (2)Impulse noise
    switch(IN_MODE)
        case 1 %weak
            disp('Weakly IN');
            GINR = 0.1;
            pb = 0.01;
        case 2 %mild
            disp('Mild IN');
            GINR = 0.01;
            pb = 0.05;
        case 3 %strong
            disp('Strongly IN');
            GINR = 0.005;
            pb = 0.08;
        otherwise
            disp('ERROR IN MODE');
            exit();
    end
    sigma = sqrt(var_v);
    imp = BG_Noise(pb,sigma,GINR,length(d));
    dn2 = dn1 + imp;
    
    % Adaptive filter
    P_up = 0.3;
    w1_hat = zeros(length(w),1);
    w2_hat = zeros(length(w),1);
    w_d_hat = zeros(length(w),1);
    w_e_hat = zeros(length(w),1);
    w_m_hat = zeros(length(w),1);
   [e1,w1_hat,update_ratio1(k,1)] = DS_LMS_Diniz(x,dn1,w1_hat,P_up,var_v);
    [e2,w2_hat,update_ratio2(k,1)] = DS_LMS_Diniz(x,dn2,w2_hat,P_up,var_v);
    [e_d,w_d_hat,update_ratio_d(k,1),detection_d(k,1),false_alarm_d(k,1)] = DS_LMS_d1(x,dn2,w_d_hat,P_up,var_v,imp);
    [e_e,w_e_hat,update_ratio_e(k,1),detection_e(k,1),false_alarm_e(k,1)] = DS_LMS_e1(x,dn2,w_e_hat,P_up,var_v,imp);
     
    NMSD1(k,:) = Normalized_Mean_Squaue_Deviation2(w,w1_hat);
    NMSD2(k,:) = Normalized_Mean_Squaue_Deviation2(w,w2_hat);
    NMSD_d(k,:) = Normalized_Mean_Squaue_Deviation2(w,w_d_hat);
    NMSD_e(k,:) = Normalized_Mean_Squaue_Deviation2(w,w_e_hat);
    % NMSD_m(k,:) = Normalized_Mean_Squaue_Deviation2(w,w_m_hat);
end

NMSD1 = sum(NMSD1,1)/ensamble;
NMSD2 = sum(NMSD2,1)/ensamble;
NMSD_d = sum(NMSD_d,1)/ensamble;
NMSD_e = sum(NMSD_e,1)/ensamble;

NMSD1 = 10*log10(NMSD1);
NMSD2 = 10*log10(NMSD2);
NMSD_d = 10*log10(NMSD_d);
NMSD_e = 10*log10(NMSD_e);


% plot
% 
% figure
% subplot('position',[0.1 0.1 0.2 0.7]);
% boxplot([update_ratio1, update_ratio2, update_ratio_d, update_ratio_e],'Labels',{'Diniz''s DS-LMS(no IN)','Diniz''s DS-LMS','Jeong''s IN Detection','Proposed Method'})
% grid on;set(gca,'FontSize',14);
% title('Update ratio')
% 
% subplot('position',[0.4 0.1 0.2 0.7]);
% boxplot([detection_d, detection_e],'Labels',{'Jeong''s IN Detection','Proposed Method'})
% set(gca,'FontSize',14);grid on;
% title('Detection rate')
% 
% subplot('position',[0.65 0.1 0.3 0.7]);
% boxplot([false_alarm_d, false_alarm_e],'Labels',{'Jeong''s IN Detection','Proposed Method'});grid on;
% set(gca,'FontSize',14);
% title('False alarm rate')
disp('***************************************');
disp(['Diniz DS_LMS(no IN) Update Ratio:   ',num2str(mean(update_ratio1))]);
disp(['Diniz DS_LMS        Update Ratio:   ',num2str(mean(update_ratio2))]);
disp(['Jeong               Update Ratio:   ',num2str(mean(update_ratio_d))]);
disp(['Proposed            Update Ratio:   ',num2str(mean(update_ratio_e))]);
disp('***************************************');
disp(['Jeong             Detection Rate:   ',num2str(mean(detection_d))]);
disp(['Proposed          Detection Rate:   ',num2str(mean(detection_e))]);
disp('***************************************');
disp(['Jeong             False Alarm   :   ',num2str(mean(false_alarm_d))]);
disp(['Proposed          False Alarm   :   ',num2str(mean(false_alarm_e))]);
disp('***************************************');

figure
plot(n,NMSD1, '-.k', 'linewidth', 2);
hold on;
plot(n,NMSD2, 'm', 'linewidth', 2);
plot(n,NMSD_d, ':r', 'linewidth', 2);
plot(n,NMSD_e, '--b', 'linewidth', 2);
set(gca,'FontSize',18);

title(['P_{up}=', num2str(P_up),'; SNR=',num2str(SNRdB),'; \Gamma=',num2str(GINR),'; P_b=',num2str(pb)]);
xlabel('Number of iterations (n)');
ylabel('NMSD (dB)');
grid;
legend('Diniz DS-LMS(no IN)','Diniz DS-LMS','Jeong IN Detection','Proposed Method');
