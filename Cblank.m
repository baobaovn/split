function [hang,lie,X] = Cblank(X)
hang=all(X==0,2);
lie=all(X==0,1);
X(hang,:)=[];
X(:,lie)=[];
hang=logical(1-hang);
lie=logical(1-lie);
end

