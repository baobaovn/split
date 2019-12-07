function [SE_Fi,x_so] = exh(h_d,h_u,w_d,h_u_d,H_si,M,D,U,sigma,q_u,q_d)
pii=zeros(1,M);
SE_Fi=0;
solo=zeros(1,M);
for i=0:2^M-1
    x=i;
    for j=1:M
        pii(j)=fix(x/2^(M-j));
        x=mod(x,2^(M-j));
    end
    x_d=pii;
    x_u=ones(1,M)-x_d;
    X_u=diag(x_u);
    X_d=diag(x_d);
    h_u_real=X_u*h_u;
    h_d_real=X_d*h_d;
    H_si_real=X_u*H_si*X_d;
 

%ÖØÐÂ¼ÆËãÔ¤±àÂë
%     if rank(X_d)~=0
%         h_inv=pinv(h_d_real');
%         for j=1:D
%             w_d(:,j)=h_inv(:,j)/norm(h_inv(:,j));
%         end
%     end
    
    
    
    
    
    
    SE_D=Cal_SE_D(w_d,q_d,q_u,h_d_real,sigma,h_u_d);
    r_u = Cal_ru(q_u,h_u,H_si,w_d,sigma,X_u,X_d,U);
    SE_U=Cal_SE_U(U,x_u,M,q_u,q_d,h_u_real,H_si_real,sigma,w_d,r_u);
    if SE_D+SE_U>=SE_Fi
        SE_Fi=SE_D+SE_U;
        x_so=x_u;
    end
    solo=[solo;pii];
end
end

