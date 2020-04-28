disp(' '); disp('clear all, close all');
clear all, close all
disp(' ');

disp('meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];');
meanfunc = {@meanSum, {@meanLinear, @meanConst}}; hyp.mean = [0.5; 1];
disp('covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);');
covfunc = {@covMaterniso, 3}; ell = 1/4; sf = 1; hyp.cov = log([ell; sf]);
disp('likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);');
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);
disp(' ');

disp('n = 20;');
n = 20;
disp('x = gpml_randn(0.3, n, 1);');
x = gpml_randn(0.3, n, 1);
disp('K = feval(covfunc{:}, hyp.cov, x);');
K = feval(covfunc{:}, hyp.cov, x);
disp('mu = feval(meanfunc{:}, hyp.mean, x);');
mu = feval(meanfunc{:}, hyp.mean, x);
disp('y = chol(K)''*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);');
y = chol(K)'*gpml_randn(0.15, n, 1) + mu + exp(hyp.lik)*gpml_randn(0.2, n, 1);
 
figure(1)
set(gca, 'FontSize', 24)
disp(' '); disp('plot(x, y, ''+'')');
plot(x, y, '+', 'MarkerSize', 12)
axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')
%print -depsc f0.eps
disp(' '); disp('Hit any key to continue...'); pause

disp(' ');
disp('nlml = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y)');
nlml = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y)
disp(' ')

disp('z = linspace(-1.9, 1.9, 101)'';');
z = linspace(-1.9, 1.9, 101)';
disp('[m s2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, z);');
[m s2 mm ss2] = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, y, z);

figure(2)
set(gca, 'FontSize', 24);
disp(' ');
disp('f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];'); 
f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)]; 
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8);');
fill([z; flipdim(z,1)], f, [7 7 7]/8);

disp('hold on; plot(z, m); plot(x, y, ''+'')')
hold on; plot(z, m, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')
%print -depsc f1.eps
disp(' '); disp('Hit any key to continue...'); pause

figure(3)
set(gca, 'FontSize', 24);
disp(' ');
disp('f = [m+2*sqrt(s2); flipdim(m-2*sqrt(s2),1)];'); 
f = [mm+2*sqrt(ss2); flipdim(mm-2*sqrt(ss2),1)]; 
disp('fill([z; flipdim(z,1)], f, [7 7 7]/8);');
fill([z; flipdim(z,1)], f, [7 7 7]/8);

disp('hold on; plot(z, m); plot(x, y, ''+'')')
hold on; plot(z, mm, 'LineWidth', 2); plot(x, y, '+', 'MarkerSize', 12)
axis([-1.9 1.9 -0.9 3.9])
grid on
xlabel('input, x')
ylabel('output, y')
%print -depsc f1.eps
disp(' '); disp('Hit any key to continue...'); pause