function [f,g]=nnloss(w,par)

% problem dimensions and parameters
d=par.dimensions;
td=par.traindata;
cl=par.classdata;

% construct weight matrices
nd=length(d)-1;
M=cell(1,nd);
b=0;
for k=1:nd
    a=b+1;
    b=a+d(k)*d(k+1)-1;
    M{k}=reshape(w(a:b),d(k+1),d(k));
end

% forward computation
L=cell(1,nd);
t=M{1}*td;
L{1}=1./(1+exp(-t));
for k=2:nd
    t=M{k}*L{k-1};
    L{k}=1./(1+exp(-t));
end

% classify result
if par.classify
    f=L{nd};
    return
end

% loss function
f=0.5*sum((L{end}-cl).^2);

% now the gradient by backpropagation
if nargout>1
    g=[];
    t=L{nd}-cl;
    for k=nd:-1:1
        h=t.*(L{k}.*(1-L{k}));
        t=M{k}'*h;
        if k>1
            G=h*L{k-1}';
        else
            G=h*td';
        end
        g=[G(:) ; g ];
    end
end

return