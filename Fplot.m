function Fplot(x,f1,f2,f3,discrete,Number,H)
%slCharacterEncoding('windows-1251');
%syms x
%x=linspace(p1,p2,100);
%p = vectorize(taylor(@(x)func(x),n+1,a));
%x=linspace(left,right,100);
%f=f(x);
%p=eval(p);

Mt=sum( (func(x)-f1).^2 )./discrete
M0=sum( (func(x)-f2).^2 )./discrete

subplot(3,2,1);
plot(x,func(x),'-b',x,f1,'-r')
title(['Degree ',num2str(Number)],'Color','b')
grid on;

legend('func(x)',...
       'Taylor''s approximation of func(x)',...
       'Location','Best');

subplot(3,2,3);
plot(x, func(x)-f1,'.',x,0*x,'*');  ylabel('mistake');
grid on;

%subplot(3,2,5);
%semilogy(x,abs(func(x)-f1),'.');  ylabel('log mistake');
%title(['Квадратична похибка рівна: ',num2str(Mt)]);
%grid on

%legend('approximation of sin(x)/x up to O(x^6)',...
%       'approximation of sin(x)/x up to O(x^8)',...
%       'approximation of sin(x)/x up to O(x^{10})',...
%       'sin(x)/x','Location','Best')
%title('Taylor Series Expansion')

%title(['Temperature is ',num2str(c),'C'])
%title(['Case number #',int2str(n)],'Color','y')
%title({'First line';'Second line'})

%C0(n+2-r)
%['Degree ',num2str(Number)]

subplot(3,2,2); 
plot(x,func(x),'-b',x,f2,'-r');  
title(['Degree ',num2str(Number)],'Color','b')
ylabel('sin(x)');
grid on;

legend('func(x)',...
       'Functional analysis approximation of func(x)',...
       'Location','Best');
  
subplot(3,2,4);
plot(x, func(x)-f2,'.',x,0*x,'*',x,f3,'r');  ylabel('mistake');
title(['Amplitude ',num2str(H)],'Color','b');
grid on;
  
subplot(3,1,3);
semilogy(x,abs(func(x)-f1),'*b',x,abs(func(x)-f2),'*r');  ylabel('log mistake');
%title({['Середня квадратична похибка Тейлора рівна: ',num2str(Mt)];['Середня квадратична похибка функціонально аналізу рівна: ',num2str(M0)]});
title({['The average square error Taylor''s equal: ',num2str(Mt)];['The average square error FA''s equal: ',num2str(M0)]});
 
%legend('похибка Тейлора',...
%       'похибка функціонально аналізу',...
%       'Location','Best');
legend('Taylor''s error',...
       'FA''s error',...
       'Location','southwest');

grid on
end