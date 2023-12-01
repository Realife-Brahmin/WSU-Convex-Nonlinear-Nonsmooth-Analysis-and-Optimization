
% solve the minimum time path problem
% use sinusoidal parameterization to order n
n=4;

% read in the velocity data array defined on 
% [0,1]x[0,1] and set the path end points
pathpar=[];
pathpar.v=readmatrix('SpeedData.csv');
[my,mx]=size(pathpar.v);
pathpar.A=[.05 .05];
pathpar.B=[.95 .95];

% add obstruction
%[yo,xo]=ndgrid(1:my,1:mx);
%yc=my/2+55;xc=mx/2+45;
%pathpar.v((yo-yc).^2+(xo-xc).^2<1600)=0.01;

% set optimization parameters
pr.objective=@pathtime;
pr.par=pathpar;
pr.x0=0.1*randn(2*n,1);  
pr.method='BFGS';
pr.linesearch='StrongWolfe';
pr.dftol=1E-8;
pr.ngtol=1E-8;
pr.dxtol=1E-8;
pr.c1=0.001;
pr.c2=0.9;
pr.m=5;
pr.maxiter=999;
pr.progress=10;

% call the optimization routine
out=optimize(pr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the optimal path overlayed on the velocity map
% first, set up graphics parameters
FontSize=24;        % size of figure fonts
LineWidth=3;        % width of path line
PointSize=100;      % path end point area
LineColor=[1 1 1];  % path color
ColorMap=jet(512);  % velocity color map
NumPathPoints=1000; % number of points (s) defining the path
FigureScale=1.4;    % figure size scale on screen

% start the figure drawing
figure('position',FigureScale*[200 200 950 800]);
% draw the velocity map
imagesc(pathpar.v)
colormap(ColorMap)
% set axes parameters
set(gca,'xtick',[],'ytick',[],'box','on')
hc=colorbar('fontsize',FontSize);
TL=get(hc,'ticklabels');
for k=1:length(TL)
    TL{k}=num2str(str2double(TL{k}),'%.1f');
end
set(hc,'ticklabels',TL);
% and add the starting and ending positions to the plot
hold on
scatter([pathpar.A(1) pathpar.B(1)]*mx,[pathpar.A(2) pathpar.B(2)]*my,...
    PointSize,LineColor,'filled')
% compute the path and add to the plot (this must match
% the computation of the objective function)
s=linspace(0,1,NumPathPoints)';
xx=(1-s)*pathpar.A(1)+s*pathpar.B(1);
yy=(1-s)*pathpar.A(2)+s*pathpar.B(2);
for k=1:n
    S=sin(k*pi*s);
    xx=xx+out.x(k,end)*S;
    yy=yy+out.x(k+n,end)*S;
end
xxs=1+xx*(mx-1);
yys=1+yy*(my-1);
plot(xxs,yys,'color',LineColor,'linewidth',LineWidth)

%saveas(gcf,'minpathex.png')

% start the figure drawing
figure('position',FigureScale*[200 200 950 800]);
% draw the velocity map
imagesc(pathpar.v)
colormap(ColorMap)
% set axes parameters
set(gca,'xtick',[],'ytick',[],'box','on')
hc=colorbar('fontsize',FontSize);
TL=get(hc,'ticklabels');
for k=1:length(TL)
    TL{k}=num2str(str2double(TL{k}),'%.1f');
end
set(hc,'ticklabels',TL);
% and add the starting and ending positions to the plot
hold on
scatter([pathpar.A(1) pathpar.B(1)]*mx,[pathpar.A(2) pathpar.B(2)]*my,...
    PointSize,LineColor,'filled')