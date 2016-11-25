clear;      % звільнити пам'ять MATLAB від усіх змінних
close all;  % закрити всі графічні вінка
%slCharacterEncoding('ibm-5347_P100-1998'); 

% налаштування параметрів графічних вікон для виведення графіків
  scrsz = get(0,'ScreenSize');                                              % отримати розміри екрану
  set(0,'defaultFigurePosition',[scrsz(3)*0/4+9 48 scrsz(3)*4/4-18 scrsz(4)-124])  % отримати розміри екрану
  set(0,'defaultTextFontName','MS Sans Serif')                              % встановити фонт для відображення кирилиці
%currentCharacterEncoding = slCharacterEncoding();
  format short g 
  
  Nf = 15;                       % порядок Функціонального аналізу
  nf = 0;                       % нижній порядок Функціонального аналізу
  n=Nf-nf;
  % константи для базису
  %       w=1:n+1;
  %d=10.^(n+4-w); M0=2.682e-020
  % [1e+013       1e+012       1e+011       1e+010       1e+009
  % 1e+008       1e+007       1e+006       1e+005        1e+004
  % 1e+003]; M0=2.682e-020
         w(1:n+1)=2;
         %u=[0.400000000000000;1.09544511501033;1.41421356237311;1.67332005306808;1.89736659609680;2.09761769609492;2.28035085438274;2.44948974968094;2.60767756751292;2.75707700133031;2.90183169193382;]
         % коефіцієнти нормованого базису
         d= 1000+0+10.^((3) + w);
  %R=[1e+000       1.e+000       1e+000       1e+000       1e+000       1e+000       1e+000       1e+000       1e+000        1.0e+000         1.e+000];
  %H(1:n+1)=1;
  H(1:n+1)=d; %d                % довільні коефіцієнти (нормуючі) для ортогонального базису
  %q=1:n;
  %H(n+1-q)=d(n+1-q);
  t=0.3;                        % затримка для паузи
  discrete=100;                 % кількість дискретів на графіку
  tol=1.0e-3;                   % похибка інтегралів
  
  Nt = 3;                       % порядок Тейлора
  nt = 0;                       % нижній порядок Тейлора
  %Nch = 6;                     % порядок Чебушева
  %nch = 0;                     % нижній порядок Чебушева
  p1=0;                         % нижні межі промужку апроксимації
  p2=(4/4)*pi*1;                % верхні межі промужку апроксимації
  p=(p2-p1)/2+p1;               % середина проміжку апроксимації
  
  %global kf, jf
  %kf=nf:Nf;                     
  %jf=nf:Nf;
  %function U=Basis2(x)
  %      global kf, jf
  %      U=basis(x,kf).*basis(x,jf);
  %end
  %function U=BasisFun(x)
  %      global kf
  %      U=func(x).*basis(x,kf);
  %end
  % отримаємо коефіцієнти полінома ортогонального базису
  
  %for k=nf:Nf
  %    for j=nf:Nf
          %Y(k+1,j+1)=k+j % відлагодження програми
          % інтеграл базисних функцій
  %        Y(Nf+1-k,Nf+1-j)=quad(@(x)Basis2(x,k,j),p1,p2);    
          %Y(k-nf+1,j-nf+1)
          %E(k+1)=k % відлагодження програми
          % інтеграл функції на базисну функцію
  %        E(Nf+1-k)=quad(@(x)BasisFun(x,k),p1,p2);
          %E(k-nf+1)
  %    end
  %end
  
  P=FullBasis(nf,Nf,p1,p2,H);
  
  'wait... ... ...'
  for k=0:n
      %for j=0:n
      %    'wait...'
      %    Y(n+1-k,n+1-j)=quad(@(x)Basis2Poly(x,P(k+1,1:end),P(j+1,1:end)),p1,p2);
      %    'wait... ...'
      %    Y1(n+1-k,n+1-j)=quad(@(x)Basis2(x,k+nf,j+nf),p1,p2);  
      %end
      'wait... ... ...'
      Y0(n+2-k)=ConvInt(P(k+1,1:end),p1,p2); %=1
      %conv(P(k+1,1:end),P(k+1,1:end));       %P(k+1,1:end)
      'wait... ... ... ...'
      E(n+2-k)=quad(@(x)BasisPolyFunc(x,P(k+1,1:end)),p1,p2,tol);
      %syms x ;
      %E1(n+2-k)=quadgk(@(x)BasisPolyFunc(x,P(k+1,1:end)),p1,p2,tol);
      E2(n+2-k)=quadl(@(x)BasisPolyFunc(x,P(k+1,1:end)),p1,p2,tol);
      E3(n+2-k)=quadv(@(x)BasisPolyFunc(x,P(k+1,1:end)),p1,p2,tol);
      %'wait... ... ... ... ...'
      %E1(n+1-k)=quad(@(x)BasisFun(x+nf,k+nf),p1,p2);
  end
  
%  Det=det(Y) % перевірка чи матриця добре обумовлена
%  C1=E/Y; % поліноміальні коефіцієнти
  C0=E./Y0; % похибка матриці відсутня через метод обрахування без застосування матриці
  C0(1)=0;
  %C0=E; % формула для нормованого базису
  %Mis0=0
%  options.TRANSA=true;
%  C=linsolve(Y,E',options); % поліноміальні коефіцієнти
%  Mis=(C'*Y-E)'  % перевірка похибки C (не ідеальна матриця)
%  Mis1=(C1*Y-E)'  % перевірка похибки C1 (не ідеальна матриця)
  %v=length(C) %Nf-nf+1
  %w = conv(C,v)
  %C2=E1/Y1; % поліноміальні коефіцієнти
%  C3=linsolve(Y1,E1',options); % поліноміальні коефіцієнти
  %C1=polyint(C)
  
  %for
      %w
  %end
  %C=C1; % швижка заміна змінної % результат індентичний
  
  %syms x
  %tf=taylor(@(x)func(x),Nt+1,x,p);
  %taylor('sin(x)',Nt,x,0)
  %pretty(tf)
  
  %taylorplot(f,a,left,right,n)
  
  %f=inline('sin(x)')
  %taylorplot(@(x)func(x),p,p1,p2,Nt); ylabel('sin(x) "-"              Tfunc(x) "-"');
  %grid on;
  
  %pause()
  %close
  %taylorplot(@(x)func(x),Nt,pi,1)
  
  %x = p1:.001:p2;
  %syms x
  %f=@(x)func(x);
  %Taylor=ftaylor(@(x)func(x),p,Nt); % формула Тейлора
  %tp = vectorize(taylor(f(x),n+1,a)); 
  
  x=linspace(p1,p2,discrete);
  
  %f=@(x)func(x);
  %tp=eval(Taylor);
  P(n+2,1:end)=0;
  %  y=0;
  y0=0;
%  y1=0;
  figure('Name','Порівняльна характеристика апроксимації функції за формулами функціонального аналізу','NumberTitle','off')

  for r=1:n+1
      Taylor=ftaylor(@(x)func(x),p,nt+r-1); % формула Тейлора
      tp=eval(Taylor);
%      y=C(n+2-r)*basisPoly(x,P(r,1:end))+y;
      
      y0=C0(n+3-r)*basisPoly(x,P(r,1:end))+y0;
%      y1=C3(n+2-r)*basis(x,r+nf-1)+y1;
      %y1=C2(n+2-r)*basis(x,r+nf-1)+y1;
%      plot(x,C(n+2-r)*basisPoly(x,P(r,1:end)));  ylabel('Basis1');
      %hold on;
%      grid on;
%      pause()
   %   plot(x,C0(n+2-r)*basisPoly(x,P(r,1:end)));  ylabel('Basis2');
      %hold on;
   %   grid on;
   %   pause(t)
   
   %v=C0(n+2-r)*basisPoly(x,P(r,1:end));
   
    %title(['Temperature is ',num2str(c),'C'])
    %title(['Case number #',int2str(n)],'Color','y')
    %title({'First line';'Second line'})
    %figure(r);
    %subplot(3,2,1);
    %plot(x,func(x),'-',x,tp1,'r')
    %grid on;
    
    Fplot(x,tp,y0,C0(n+2-r)*basisPoly(x,P(r+1,1:end)),discrete,r-1,H(n+2-r));
    
    %figure('Name','Simulation Plot Window','NumberTitle','off')
 %    subplot(3,1,1); 
 % plot(x,y0,'-',x, sin(x),'-');  ylabel('sin(x) "-"              FAfunc(x) "-"');
 % grid on;
  
 % subplot(3,1,2);
 % plot(x, sin(x)-y0,'.',x,0*x,'*');  ylabel('mistake');
 % grid on;
  
 % subplot(3,1,3);
 % semilogy(x,abs(sin(x)-y0),'.');  ylabel('log mistake');
 % grid on
  %pause()
  %close
  
  
  end
  close
%  C
  C0'
%  C1'
  
%  M=sum( (sin(x)-y).^2 ) 
  %y=polyval(C,x);
  %plot(y);
  
  
  
%  subplot(3,1,1); 
  %plot(x,y,'.');  %ylabel('func(x) .');
  %grid on;
  %hold on;
%  plot(x,y,'.',x, sin(x),'-');  ylabel('sin(x) "-"              func(x) "."');
%  grid on;
  
%  subplot(3,1,2);
%  plot(x, sin(x)-y,'.',x,0*x,'*');  ylabel('mistake');
%  grid on;
  
%  subplot(3,1,3);
%  semilogy(x,abs(sin(x)-y),'.');  ylabel('log mistake');
%  grid on
%  pause()
%  close
  
   

  
  
  
  
%  for r=1:n+1
%      plot(x,C3(n+2-r)*basis(x,r+nf-1));  ylabel('Basis3');
%      grid on;
%      pause(t)
%  end
%  close
%  M1=sum( (sin(x)-y1).^2 ) 
  
  
%  subplot(3,1,1); 
%  plot(x,y1,'.',x, sin(x),'-');  ylabel('sin(x) "-"              func(x) "."');
%  grid on;
  
%  subplot(3,1,2);
%  plot(x, sin(x)-y1,'.',x,0*x,'*');  ylabel('mistake');
%  grid on;
  
%  subplot(3,1,3);
%  semilogy(x,abs(sin(x)-y1),'.');  ylabel('log mistake');
%  grid on
%  pause()
%  close
  

%legend('approximation of sin(x)/x up to O(x^6)',...
%       'approximation of sin(x)/x up to O(x^8)',...
%       'approximation of sin(x)/x up to O(x^{10})',...
%       'sin(x)/x','Location','Best')
%title('Taylor Series Expansion')


  
  %plot(x,sin(x),'-',x,y,'.');  ylabel('sin(x) "-"     func(x) "."');
  %grid on;
  %pause()
  %close
  %x = logspace(-2,1);
  
  %abs(sin(x)-y);
  
  
  
  %Y(0,1)=1;
  %Y(kf+1,jf+1)=quad(@(x)Basis2(x,kf,jf),p1,p2);% інтеграл базисних функцій
  %E(kf+1)=quad(@(x)BasisFun(x,kf),p1,p2);      % інтеграл функції на базисну функцію
  %U=basis(x,kf);
  
  %function f = fint(x)
  %f=sin(x);
  
  
  %[x,y] = basis(x,y);
  %f(x)=sin(x);
  
  %147 тейлор
  %https://www.youtube.com/watch?v=JvZo5HM2e5g&list=PLnbT_JFJiDUrXD0RKZCYmh
  %yWntVy_PFgi&index=147
  %f=sin(x);
  %tf=taylor(f);
  %pretty(tf)
  
  %taylortool
  
  %syms k
  %s*symsum(f, k , 1 Inf)
  
  %diff(I,x,l);
  %I=int(f,x);
  %pretty(I);
  
  %коефіцієнти полінома
  
  %x1=fzero(func, x0);
  %x2=fzero(func, [3 4]);
  
  %149 символьне рішення
  
  %x3=fminbnd(func, -1, 1);
  %[x4,f,flag]=fminsearch(@func, [-0.5]);
  %fminbnd();
  
  %format long
  %integ=quad(@func, -1,1);
  %fun=inline('x.^3-0.9');
   
  %p=[1 0 1 0 0 1];
  %x5=polyder(p);
  %q=[1 2 3];
  %[n, x6]=polyder(p,q);
  
  %options.TRANSA=true;
  
  %x7=linsolve(a,b,options); %74
  
  %plot(x_inp);
  %pause;
  %close
  
  %figure(5);
  %n=1:50;
  %subplot(2,1,1); plot(t(n), cos_Wc_t(n),'-o');  ylabel('cos(Wc*t)');
  %subplot(2,1,2); plot(t(n), sin_Wc_t(n),'-o');  ylabel('-sin(Wc*t)');
  %pause;
  %close