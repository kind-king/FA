clear;      % �������� ���'��� MATLAB �� ��� ������
close all;  % ������� �� �������� ����

% ������������ ��������� ��������� ���� ��� ��������� �������
  scrsz = get(0,'ScreenSize');                                              % �������� ������ ������
  set(0,'defaultFigurePosition',[scrsz(3)/4 50 scrsz(3)*3/4 scrsz(4)-130])  % �������� ������ ������
  set(0,'defaultTextFontName','MS Sans Serif')                              % ���������� ���� ��� ����������� ��������

  format short g 
  
  Nf = 10;                       % ������� ��������������� ������
  nf = 0;                       % ������ ������� ��������������� ������
  n=Nf-nf;
  % ��������� ��� ������
  %        1:n+1
  H(1:n+1)=100;
  t=0.3;                        % �������� ��� �����
  
  Nt = 6;                       % ������� �������
  %nt = 0;                      % ������ ������� �������
  %Nch = 6;                     % ������� ��������
  %nch = 0;                     % ������ ������� ��������
  p1=0;                         % ����� ��� �������� ������������
  p2=(1/4)*pi;                  % ������ ��� �������� ������������
  p=(p2-p1)/2+p1;               % �������� ������� ������������
  
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
  % �������� ����������� ������� �������������� ������
  
  %for k=nf:Nf
  %    for j=nf:Nf
          %Y(k+1,j+1)=k+j % ������������ ��������
          % �������� �������� �������
  %        Y(Nf+1-k,Nf+1-j)=quad(@(x)Basis2(x,k,j),p1,p2);    
          %Y(k-nf+1,j-nf+1)
          %E(k+1)=k % ������������ ��������
          % �������� ������� �� ������� �������
  %        E(Nf+1-k)=quad(@(x)BasisFun(x,k),p1,p2);
          %E(k-nf+1)
  %    end
  %end
  P=FullBasis(nf,Nf,p1,p2,H);
  'wait... ... ...'
  for k=0:n
      for j=0:n
          'wait...'
          Y(n+1-k,n+1-j)=quad(@(x)Basis2Poly(x,P(k+1,1:end),P(j+1,1:end)),p1,p2);
          'wait... ...'
          Y1(n+1-k,n+1-j)=quad(@(x)Basis2(x,k+nf,j+nf),p1,p2);  
      end
      'wait... ... ...'
      Y0(n+1-k)=ConvInt(P(k+1,1:end),p1,p2);
      %conv(P(k+1,1:end),P(k+1,1:end));       %P(k+1,1:end)
      'wait... ... ... ...'
      E(n+1-k)=quad(@(x)BasisPolyFunc(x,P(k+1,1:end)),p1,p2);
      'wait... ... ... ... ...'
      E1(n+1-k)=quad(@(x)BasisFun(x+nf,k+nf),p1,p2);
  end
  
  Det=det(Y) % �������� �� ������� ����� ����������
  C1=E/Y; % ����������� �����������
  C0=E./Y0; % ������� ������� ������� ����� ����� ����������� ��� ������������ �������
  Mis0=0
  options.TRANSA=true;
  C=linsolve(Y,E',options); % ����������� �����������
  Mis=(C'*Y-E)'  % �������� ������� C (�� �������� �������)
  Mis1=(C1*Y-E)'  % �������� ������� C1 (�� �������� �������)
  %v=length(C) %Nf-nf+1
  %w = conv(C,v)
  C2=E/Y; % ����������� �����������
  C3=linsolve(Y1,E1',options); % ����������� �����������
  %C1=polyint(C)
  
  %for
      %w
  %end
  %C=C1; % ������ ����� ����� % ��������� �����������
  
  x = p1:.001:p2;
  y=0;
  y0=0;
  y1=0;
  for r=1:n+1
      y=C(n+2-r)*basisPoly(x,P(r,1:end))+y;
      y0=C0(n+2-r)*basisPoly(x,P(r,1:end))+y0;
      y1=C3(n+2-r)*basis(x,r+nf-1)+y1;
      %y1=C2(n+2-r)*basis(x,r+nf-1)+y1;
      plot(x,C(n+2-r)*basisPoly(x,P(r,1:end)));  ylabel('Basis1');
      %hold on;
      grid on;
      pause()
      plot(x,C0(n+2-r)*basisPoly(x,P(r,1:end)));  ylabel('Basis2');
      %hold on;
      grid on;
      pause()
  end
  close
  C
  C0'
  C1'
  M0=sum( (sin(x)-y0).^2 )
  M=sum( (sin(x)-y).^2 ) 
  %y=polyval(C,x);
  %plot(y);
  
  figure(5);
  
  subplot(3,1,1); 
  %plot(x,y,'.');  %ylabel('func(x) .');
  %grid on;
  %hold on;
  plot(x,y,'.',x, sin(x),'-');  ylabel('sin(x) "-"              func(x) "."');
  grid on;
  
  subplot(3,1,2);
  plot(x, sin(x)-y,'.',x,0*x,'*');  ylabel('mistake');
  grid on;
  
  subplot(3,1,3);
  semilogy(x,abs(sin(x)-y),'.');  ylabel('log mistake');
  grid on
  pause()
  close
  
   
  subplot(3,1,1); 
  plot(x,y0,'.',x, sin(x),'-');  ylabel('sin(x) "-"              func(x) "."');
  grid on;
  
  subplot(3,1,2);
  plot(x, sin(x)-y0,'.',x,0*x,'*');  ylabel('mistake');
  grid on;
  
  subplot(3,1,3);
  semilogy(x,abs(sin(x)-y0),'.');  ylabel('log mistake');
  grid on
  pause()
  close
  
  
  
  for r=1:n+1
      plot(x,C3(n+2-r)*basis(x,r+nf-1));  ylabel('Basis3');
      grid on;
      pause(t)
  end
  close
  M1=sum( (sin(x)-y1).^2 ) 
  
  
  subplot(3,1,1); 
  plot(x,y1,'.',x, sin(x),'-');  ylabel('sin(x) "-"              func(x) "."');
  grid on;
  
  subplot(3,1,2);
  plot(x, sin(x)-y1,'.',x,0*x,'*');  ylabel('mistake');
  grid on;
  
  subplot(3,1,3);
  semilogy(x,abs(sin(x)-y1),'.');  ylabel('log mistake');
  grid on
  pause()
  close
  
  
  %plot(x,sin(x),'-',x,y,'.');  ylabel('sin(x) "-"     func(x) "."');
  %grid on;
  %pause()
  %close
  %x = logspace(-2,1);
  
  %abs(sin(x)-y);
  
  %syms x
  %tf=taylor('func(x)',Nt,x,0)
  %taylor('sin(x)',Nt,x,0)
  %pretty(tf)
  
  %Y(0,1)=1;
  %Y(kf+1,jf+1)=quad(@(x)Basis2(x,kf,jf),p1,p2);% �������� �������� �������
  %E(kf+1)=quad(@(x)BasisFun(x,kf),p1,p2);      % �������� ������� �� ������� �������
  %U=basis(x,kf);
  
  %function f = fint(x)
  %f=sin(x);
  
  
  %[x,y] = basis(x,y);
  %f(x)=sin(x);
  
  %147 ������
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
  
  %����������� �������
  
  %x1=fzero(func, x0);
  %x2=fzero(func, [3 4]);
  
  %149 ��������� ������
  
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