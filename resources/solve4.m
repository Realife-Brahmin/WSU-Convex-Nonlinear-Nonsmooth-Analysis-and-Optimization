
% liver data analysis

% load data
A=textread('liver.csv','%s','delimiter',',');          %#ok
A=reshape(A,11,[]);
[r,c]=size(A);

% ceate numerical data arrays
%   data is 10 by ~580 maximum normalized data values
%   class is 1 by ~580 0-1 classifier
data=zeros(r-1,c);
class=zeros(1,c);
idx=[];
for k=1:r
    if k==2
        for j=1:c
            if ~isempty(A{2,j})
                data(k,j)=(length(A{2,j})-4)/2;
            else
                idx=[idx j];                           %#ok
            end
        end
    elseif k==11
        for j=1:c
            try
                class(j)=str2num(A{11,j})-1;           %#ok
            catch 
                idx=[idx j];                           %#ok
            end
        end
    else
        for j=1:c
            try
                data(k,j)=str2num(A{k,j});             %#ok
            catch
                idx=[idx j];                           %#ok
            end
        end
    end
end

% remove instances with missing data
idx=unique(idx);
data(:,idx)=[];
class(idx)=[];

% normalize data to [0,1]
MinValues=min(data,[],2);
MaxValues=max(data,[],2);
data=(data-MinValues)./(MaxValues-MinValues);
[r,c]=size(data);

% subsample features
features=(1:10)';
data=data(features,:);
[r,c]=size(data);

% construct training data 
m=100;
cdx=find(class==1);
ddx=find(class==0);
par.traindata=[data(:,cdx(1:m)) data(:,ddx(1:m))];
sh=1/3;
par.classdata=[class(cdx(1:m)) class(ddx(1:m))]*(1-2*sh)+sh;

% set up NN parameters and initial weights
par.dimensions=[10 10 10 1];
q=length(par.dimensions);
matrixsizes=par.dimensions(1:q-1).*par.dimensions(2:q);
numweights=sum(matrixsizes);
sz=1/sqrt(par.dimensions(1));
w=2*sz*rand(numweights,1)-sz;

% call the optimization
pr.objective=@nnloss;
pr.par=par;
pr.x0=w;
pr.method='BFGS';
pr.linesearch='StrongWolfe';
pr.maxiter=100000;
pr.dftol=1E-9;
pr.ngtol=1E-9;
pr.dxtol=1E-9;
pr.maxcond=1000;
pr.progress=1000;

% This is the call to the optimizer
pr.par.classify=false;
out=optimize(pr);

% This is a call to classify the training data
pr.par.classify=true;
res=pr.objective(out.x(:,end),pr.par);

% Show the confusion matrix for the training data
NN=sum(res<0.5&par.classdata<0.5);
NP=sum(res<0.5&par.classdata>=0.5);
PN=sum(res>=0.5&par.classdata<0.5);
PP=sum(res>=0.5&par.classdata>=0.5);
ConfusionMatrix=[NN PN ; NP PP]                %#ok

% Classify the test data and show the confusion matrix
pr.par.traindata=[data(:,cdx(m+1:end)) data(:,ddx(m+1:end))];
pr.par.classdata=[class(cdx(m+1:end)) class(ddx(m+1:end))]*1/3+1/3;
test=pr.objective(out.x(:,end),pr.par);
NN=sum(test<0.5&pr.par.classdata<0.5);
NP=sum(test<0.5&pr.par.classdata>=0.5);
PN=sum(test>=0.5&pr.par.classdata<0.5);
PP=sum(test>=0.5&pr.par.classdata>=0.5);
ConfusionMatrix=[NN PN ; NP PP]               %#ok


