

pr.objective=@rosenbrock;
pr.par=10;
pr.x0=[0;0;0;0;0;0;0;0;0];
pr.x0=zeros(18,1);

pr.method='CG';
pr.linesearch='Armijo';
pr.maxiter=1000;
pr.progress=50;

out=optimize(pr);

ShowResults(out);