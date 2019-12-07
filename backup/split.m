clc;
clear all;
alpha=3.7;
c=1;
pd=1;                      %发送功率
d0=300;
shaVar=0;
SNR=50;
bandwidth=10^7;
q_u_dB=23;
q_u=power(10,q_u_dB/10)/1000;
q_d_dB=30;
q_d=power(10,q_d_dB/10)/1000;

% sigma=q_d./power(10,SNR/10);
sigma=10^-10;
sigma_si_dB=[-100:10:-50];
sigma_si=power(10,sigma_si_dB/10);
d=50;                    %小区半径
M=8;%中心基站天线数
U=4;% 上行用户数
D=6;%下行用户数
user_number=U+D;   %总用户数量
loop=600;            %仿真循环次数

h_u=zeros(M,U);
h_d=zeros(M,D);
kk=10^(-120/10);
loop=600;
looo=1;
for m=1:loop
    for k=1:D
        userloc_d(k,:)=Ud6(d);    %产生随机用户位置
        dis_d(k)=norm(userloc_d(k,:));
        lamda(k)=c*dis_d(k)^(-alpha);
        %lamda(k)=1;                 %%计算大尺度衰落
        h_d(:,k) = lamda(k)*sqrt(1/2)*(randn(1,M) +1i*randn(1,M));      %生成信道向量
    end
    for k=1:U
        userloc_u(k,:)=Ud6(d);    %产生随机用户位置
        dis_u(k)=norm(userloc_u(k,:));
        lamda(k)=c*dis_u(k)^(-alpha);                %%计算大尺度衰落
        %lamda(k)=1;
        h_u(:,k) = lamda(k)*sqrt(1/2)*(randn(1,M) +1i*randn(1,M));      %生成信道向量
    end
    
    for k=1:D
        for s=1:U
            dis_d_u(k,s)=norm(userloc_d(k,:)-userloc_u(s,:));
            lamda(k)=c*dis_d_u(k,s)^(-alpha);                %%计算大尺度衰落
            %lamda(k)=1;
            h_u_d(k,s) = lamda(k)*sqrt(1/2)*(randn(1,1) +1i*randn(1,1));      %生成信道向量
        end
    end
    
    
    var_c=kk*q_u;
    symbol_d=randn(1,D)';
    H_si=sqrt(0.5*sigma_si(1))+diag(sqrt(sigma_si(1)/2)*(randn(1,M) +1i*randn(1,M)));
    
    
    %天线配置
    %split
    x_d=zeros(1,M);
    x_d(1:M/2)=ones(1,M/2);
    x_u=ones(1,M)-x_d;
    
    
    %大天线配置
    X_u=diag(x_u);
    X_d=diag(x_d);
    
    h_u_real=X_u*h_u;
    h_d_real=X_d*h_d;
    H_si_real=X_u*H_si*X_d;
    
    %噪声
    % c_d_power=diag(w_d*w_d');
    % c_u_power=kk*q_u;
    eta_u_power=sigma;
    % c_u=sqrt(c_u_power).*randn(U,1);
    % c_d=sqrt(c_d_power).*randn(M,1);
    eta_u=sqrt(eta_u_power).*randn(M,1);
    eta_d=sqrt(eta_u_power).*randn(D,1);
    
    
    %波束赋形
    
    
    h_inv=pinv(h_d_real');
    for i=1:D
        w_d(:,i)=h_inv(:,i)/norm(h_inv(:,i));
    end
    %     H=h_d_real';
    %     w_d=H'*inv(H*H');
    %     w_d=4/trace(w_d*w_d')*w_d;
    %下行用户频谱效率
%     for i=1:D
%         ss=w_d*w_d'-w_d(:,i)*w_d(:,i)';
%         SINR_D(i)=(q_d*(h_d_real(:,i)'*w_d(:,i))*(w_d(:,i)'*h_d_real(:,i)))/(sigma+q_d*real(h_d_real(:,i)'*ss*h_d_real(:,i))+q_u*norm(h_u_d(i,:))^2);
%         SE_D(i)=log2(1+SINR_D(i));
%     end
    
    SE_D_Fi(m)= Cal_SE_D(w_d,q_d,q_u,h_d_real,sigma,h_u_d);
    %上行用户频谱效率
    %         for i=1:U
    %             [row,clo,mid]=Cblank(q_u*(h_u_real*h_u_real')+H_si_real*(w_d*w_d')*H_si_real'+sigma*X_u);
    %             mid=pinv(mid);
    %             Z=zeros(M,M);
    %             Z(row,clo)=mid;
    %             r_u(:,i)=sqrt(q_u)*Z*h_u_real(:,i);
    %             inner=0;
    %
    %             for j=1:length(u_list)
    %                 inner=inner+q_d*H_si_real(u_list(j),:)*(w_d*w_d')*H_si_real(u_list(j),:)';
    %             end
    %             SINR_U(i)=q_u*norm(r_u(:,i)'*h_u_real(:,i))^2/(q_u*r_u(:,i)'*(h_u_real*h_u_real'-h_u_real(:,i)*h_u_real(:,i)')*r_u(:,i)+inner+sigma);
    %             SE_U(i)=log2(1+SINR_U(i));
    %         end
    
    [SE_U_Fi(m),r_u]=Cal_SE_U(U,x_u,M,q_u,q_d,h_u_real,H_si_real,sigma,w_d);
    
    Gamma=q_u*h_u*h_u';
    Lamda_u=zeros(M,M);
    Lamda_d=zeros(M,M);
    a_u=zeros(M,1);
    a_d=zeros(M,1);
    for lam=1:U
        Lamda_u=Lamda_u+diag(r_u(:,lam)')*Gamma*diag(r_u(:,lam));
        a_u=a_u+sqrt(q_u)*real(diag(r_u(:,lam)*h_u(:,lam)'))-sigma/2*diag(r_u(:,lam)*r_u(:,lam)');
    end
    
%     for lan=1:D
%         Lamda_d
%     end
    
    
    
    
    
    
    
end
cdfplot(SE_D_Fi+SE_U_Fi);
hold on;