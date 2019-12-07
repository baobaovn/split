function [r_d] = Cal_rd(h_d,X_d,w_d,h_u_d,sigma,D_num)
h_d_real=X_d*h_d;
for i=1:D_num
    r_d(i)=h_d_real(:,i)'*w_d(:,i)*zeroinv(h_d_real(:,i)'*(w_d*w_d')*h_d_real(:,i)+h_u_d(i,:)*h_u_d(i,:)'+sigma);
end
end

