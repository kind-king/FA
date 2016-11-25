function p=ftaylor(f,a,n)
%TAYLORPLOT: MATLAB function M-file that takes as input
%a function in inline form, a center point a, a left endpoint,
%a right endpoint, and an order, and plots
%the Taylor polynomial along with the function at the given order.
syms x
p = vectorize(taylor(f(x),n+1,a));
x=linspace(left,right,100);
f=f(x);
p=eval(p);
figure('Name','Taylor approximations','NumberTitle','off');
plot(x,f,x,p,'r')
end