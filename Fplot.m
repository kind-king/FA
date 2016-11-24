function Fplot(x,f1,f2,C)
%syms x
%x=linspace(p1,p2,100);
%p = vectorize(taylor(@(x)func(x),n+1,a));
%x=linspace(left,right,100);
%f=f(x);
%p=eval(p);

Mt=sum( (func(x)-f1).^2 )
M0=sum( (func(x)-f2).^2 )

subplot(3,2,1);
plot(x,func(x),'-',x,f1,'r')
title({'First line';'Second line'})
grid on;

legend('sin(x)',...
       'approximation of sin(x)',...
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


subplot(3,2,2); 
plot(x,f2,'-',x, func(x),'r');  
title({'First line';'Second line'})
%ylabel('sin(x) "-"              func(x) "-"');
grid on;
  
subplot(3,2,4);
plot(x, func(x)-f2,'.',x,0*x,'*');  ylabel('mistake');
title(['Amplitude ',num2str(C)],'Color','b');
grid on;
  
subplot(3,1,3);
semilogy(x,abs(func(x)-f1),'.',x,abs(func(x)-f2),'.');  ylabel('log mistake');
title({['Квадратична похибка Тейлора рівна: ',num2str(Mt)];['Квадратична похибка функціонально аналізу рівна: ',num2str(M0)]});
legend('похибка Тейлора',...
       'похибка функціонально аналізу','Best');
grid on
end