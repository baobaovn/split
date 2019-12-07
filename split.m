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
sigma=10e-11;
sigma_si_dB=[-100:10:-50];
sigma_si=power(10,sigma_si_dB/10);
d=50;                    %小区半径
M=8;%中心基站天线数
U=4;% 上行用户数
D=6;%下行用户数
user_number=U+D;   %总用户数量
loop=100;            %仿真循环次数

h_u=zeros(M,U);
h_d=zeros(M,D);
kk=10^(-120/10);
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
            h_u_d(k,s) = lamda(k)*sqrt(1/2)*(randn(1,1) +1i*randn(1,1));     
            %生成信道向量
        end
    end
    
    
    var_c=kk*q_u;
    symbol_d=randn(1,D)';
    H_si=sqrt(0.5*sigma_si(1))+diag(sqrt(sigma_si(1)/2)*(randn(1,M)+...
        1i*randn(1,M)));
    
    %预编码
    %zf
%     h_inv=pinv(h_d');
%     for i=1:D
%         w_d(:,i)=h_inv(:,i)/norm(h_inv(:,i));
%     end
    %mrc
    for i=1:D
       w_d(:,i)=h_d(:,i)/norm(h_d(:,i));
    end
    
    %天线配置
    %split
    x_d=zeros(1,M);
    x_d(1:M/2)=ones(1,M/2);
    x_u=ones(1,M)-x_d;
    
    X_u=diag(x_u);
    X_d=diag(x_d);
    
    h_u_real=X_u*h_u;
    h_d_real=X_d*h_d;
    H_si_real=X_u*H_si*X_d;
    
    %波束赋形
    
    
    %split下行用户频谱效率
%     h_inv=pinv(h_d_real');
%     for i=1:D
%         w_d_spl(:,i)=h_inv(:,i)/norm(h_inv(:,i));
%     end
    SE_D_Fi_spl(m)= Cal_SE_D(w_d,q_d,q_u,h_d_real,sigma,h_u_d);
    %split上行用户频谱效率
    r_u = Cal_ru(q_u,h_u,H_si,w_d,sigma,X_u,X_d,U);
    SE_U_Fi_spl(m)=Cal_SE_U(U,x_u,M,q_u,q_d,h_u_real,...
        H_si_real,sigma,w_d,r_u);
    
    
    
    %文中算法
    [SE_U_Fi_pro(m),SE_D_Fi_pro(m),x_u_prop] = Rlx(h_d,h_u,w_d,h_u_d,...
        H_si,M,D,U,sigma,q_u,q_d);
    
    SE_Fi_spl(m)=SE_D_Fi_spl(m)+SE_U_Fi_spl(m);
    SE_Fi_pro(m)=SE_U_Fi_pro(m)+SE_D_Fi_pro(m);
    %穷举法
    [SE_Fi_exh(m),x_so] = exh(h_d,h_u,w_d,h_u_d,H_si,M,D,U,sigma,q_u,q_d);
   
    
    
    m
end
cdfplot(SE_Fi_spl);
hold on;
cdfplot(SE_Fi_pro);
hold on;
cdfplot(SE_Fi_exh);
legend('split','prop','exh');