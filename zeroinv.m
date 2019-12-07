function [Z] = zeroinv(X)
[row,clo,mid]=Cblank(X);
M=size(X,1);
mid=pinv(mid);
Z=zeros(M,M);
Z(clo,row)=mid;
end

