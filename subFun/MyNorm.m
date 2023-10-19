function  [X, normal] =MyNorm(x)
%%
x = double(x);

n=size(x,1);

normal.xm=mean(x);

x=x-repmat(normal.xm,n,1);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
normal.xscale=sqrt(sum(sum(x.^2,2))/n);

X=x/normal.xscale;


%%
% D = double(x);
% 
% n=size(D,1);
% 
% normal.xm=mean(D);
% 
% D=D-repmat(normal.xm,n,1);
% 
% 
% normal.xscale=sqrt(sum(D.^2)/n);
% 
% 
% X=D./repmat(normal.xscale,n,1) ;
