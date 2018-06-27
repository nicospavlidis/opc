function [v,fval,idx] = lda2m(X)
%Interface to ldakmeans that allows it to be used as generic projection pursuit function to perform divisive clustering
%[V,FVAL,IDX] = LDA2M(X)
[idx,v,fval,C,iter] = ldakmeans(X,2);
end
