function [SE_U_pro,SE_D_pro,x_u_prop] = Rlx(h_d,h_u,w_d,h_u_d,H_si,M,D,U,sigma,q_u,q_d)
alpha=1;
for dif=1:3
    x_u_old=zeros(M,1);
    x_u_new=rand(M,1);
    loop=1;
    rou=0.5;
    while norm(x_u_old-x_u_new)>10e-3
        x_u_old=x_u_new;
        Xp_u=diag(x_u_old);
        Xp_d=eye(M)-Xp_u;
        r_d=Cal_rd(h_d,Xp_d,w_d,h_u_d,sigma,D);
        r_u = Cal_ru(q_u,h_u,H_si,w_d,sigma,Xp_u,Xp_d,U);
        
        Gamma=q_u*h_u*h_u';
        Lamda_u=zeros(M,M);
        Lamda_d=zeros(M,M);
        a_u=zeros(M,1);
        a_d=zeros(M,1);
        for lam=1:U
            Lamda_u=Lamda_u+diag(r_u(:,lam)')*Gamma*diag(r_u(:,lam));
            a_u=a_u+sqrt(q_u)*real(diag(r_u(:,lam)*h_u(:,lam)'))-sigma/2*diag(r_u(:,lam)*r_u(:,lam)');
        end
        theta=cell(1,D);
        r_dd=r_d*r_d';
        for lan=1:D
            theta{lan}=r_dd*h_d(:,lan)*h_d(:,lan)';
            Lamda_d=Lamda_d+diag(w_d(:,lan))*theta{lan}*diag(w_d(:,lan))';
            a_d=a_d+real(r_d(lan)'*diag(h_d(:,lan)*w_d(:,lan)'));
        end
        Lamda=real(Lamda_u+Lamda_d);
        
        big_R=r_u*r_u';
        big_W=w_d*w_d';
        Grad=(big_R.'*Xp_u*conj(H_si)*(eye(M)-Xp_u)*(big_W.')*(eye(M)-Xp_u)*(H_si.')+conj(H_si)*(eye(M)-Xp_u)*(big_W.')*(eye(M)-Xp_u)*(H_si.')*Xp_u*(big_R.'))-(H_si.'*Xp_u*big_R.'*Xp_u*conj(H_si)*(eye(M)-Xp_u)*big_W.'+big_W.'*(eye(M)-Xp_u)*H_si.'*Xp_u*big_R.'*Xp_u*conj(H_si));
        bb=a_u+real(Lamda*ones(M,1))-a_d-0.5*diag(Grad);
        b=real(bb');
        
        cvx_begin quiet
        variables x(M,1);
        minimize (x.'*Lamda*x-2*b*x+0.5*alpha*norm(x-x_u_old));
        subject to
        x>=0;
        x<=1;
        cvx_end
        x_u_new=x_u_old+rou*(x-x_u_old);
        tongji(loop,:)=x_u_new;
        loop
        loop=loop+1;
        if loop>=10
            break;
        end
    end
    x_u_prop=round(x_u_new);
    X_u=diag(x_u_prop);
    X_d=eye(M)-diag(x_u_prop);
    h_u_real=X_u*h_u;
    h_d_real=X_d*h_d;
    H_si_real=X_u*H_si*X_d;
    
    %ÖØÐÂ¼ÆËãÔ¤±àÂë
    % if rank(X_d)~=0
    %     h_inv=pinv(h_d_real');
    %     for i=1:D
    %         w_d(:,i)=h_inv(:,i)/norm(h_inv(:,i));
    %     end
    %     r_u = Cal_ru(q_u,h_u,H_si,w_d,sigma,X_u,X_d,U);
    %     SE_D_Fi_pro= Cal_SE_D(w_d,q_d,q_u,h_d_real,sigma,h_u_d);
    % else
    %     SE_D_Fi_pro=0;
    % end
    
    
    SE_D_Fi_pro(dif)=Cal_SE_D(w_d,q_d,q_u,h_d_real,sigma,h_u_d);
    SE_U_Fi_pro(dif)=Cal_SE_U(U,x_u_prop,M,q_u,q_d,h_u_real,H_si_real,sigma,w_d,r_u);
    SE_Fi_pro(dif)=SE_D_Fi_pro(dif)+SE_U_Fi_pro(dif);
end
[value,index]=max(SE_Fi_pro);
SE_U_pro=SE_U_Fi_pro(index);
SE_D_pro=SE_D_Fi_pro(index);
end

