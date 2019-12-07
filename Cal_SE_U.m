function [SE_U_Fi] = Cal_SE_U(U_num,x_u,M,q_u,q_d,h_u_real,H_si_real,sigma,w_d,r_u)
%�����������û���������������ţ��������������й��ʣ����й��ʣ������ŵ���Ч���󣬵�Ч�����ž����������������ߵ�Ч����Ԥ�������,���ջ�
u_list=find(x_u~=0);
X_u=diag(x_u);

for i=1:U_num
%     Z=zeroinv(q_u*(h_u_real*h_u_real')+H_si_real*(w_d*w_d')*H_si_real'+sigma*X_u);
% %     [row,clo,mid]=Cblank(q_u*(h_u_real*h_u_real')+H_si_real*(w_d*w_d')*H_si_real'+sigma*X_u);
% %     mid=pinv(mid);
% %     Z=zeros(M,M);
% %     Z(row,clo)=mid;
%     r_u(:,i)=sqrt(q_u)*Z*h_u_real(:,i);
    inner=0;
    
    for j=1:length(u_list)
        inner=inner+q_d*H_si_real(u_list(j),:)*(w_d*w_d')*H_si_real(u_list(j),:)';
    end
    SINR_U(i)=q_u*norm(r_u(:,i)'*h_u_real(:,i))^2/(q_u*r_u(:,i)'*(h_u_real*h_u_real'-h_u_real(:,i)*h_u_real(:,i)')*r_u(:,i)+inner+sigma);
    SE_U(i)=log2(1+SINR_U(i));
end
SE_U_Fi=sum(SE_U);
end

