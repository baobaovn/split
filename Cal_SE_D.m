function SE_D_Fi = Cal_SE_D(w_d,q_d,q_u,h_d_real,sigma,h_u_d)
D=size(h_d_real,2);
for i=1:D
    ss=w_d*w_d'-w_d(:,i)*w_d(:,i)';
    SINR_D(i)=(q_d*(h_d_real(:,i)'*w_d(:,i))*(w_d(:,i)'*h_d_real(:,i)))/(sigma+q_d*real(h_d_real(:,i)'*ss*h_d_real(:,i))+q_u*norm(h_u_d(i,:))^2);
    SE_D(i)=log2(1+SINR_D(i));
end

SE_D_Fi=sum(SE_D);
end

