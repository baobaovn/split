function [r_u] = Cal_ru(q_u,h_u,H_si,w_d,sigma,X_u,X_d,U_num)
h_u_real=X_u*h_u;
H_si_real=X_u*H_si*X_d;
for i=1:U_num
    Z=zeroinv(q_u*(h_u_real*h_u_real')+H_si_real*(w_d*w_d')*H_si_real'+sigma*X_u);
    r_u(:,i)=sqrt(q_u)*Z*h_u_real(:,i);
end
end

