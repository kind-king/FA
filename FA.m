clear;      % звільнити пам'ять MATLAB від усіх змінних
close all;  % закрити всі графічні вінка

%currentCharacterEncoding = slCharacterEncoding();
slCharacterEncoding('windows-1251'); 

% налаштування параметрів графічних вікон для виведення графіків
  scrsz = get(0,'ScreenSize');                                              % отримати розміри екрану
  set(0,'defaultFigurePosition',[scrsz(3)*0/4+9 48 scrsz(3)*4/4-18 scrsz(4)-124])  % отримати розміри екрану
  set(0,'defaultTextFontName','MS Sans Serif')                              % встановити фонт для відображення кирилиці

  format short g 
  
  method = 0;                     % 1 2 3 0
  
  Nf = 10;                      % порядок Функціонального аналізу
  nf = 0;                       % нижній порядок Функціонального аналізу
  n=Nf-nf;
  w(1:n+1)=2;                   % константи для базису
  d= 1+0+10.^((0) + w);         % коефіцієнти нормованого базису
  H(1:n+1)=10000; %d            % довільні коефіцієнти (нормуючі) для ортогонального базису
  %q=1:n;
  %H(n+1-q)=d(n+1-q);
  t=0.3;                        % затримка для паузи
  discrete=100;                 % кількість дискретів на графіку
  tol=1.0e-3;                   % похибка інтегралів
  AbsTol=1e-12;
  RelTol=1e-8;                   
  %'RelTol',1e-8,'AbsTol',1e-12
  discreteInt=10000;            % кільксть інтервалів всього проміжку для інтегрування 
  
  Nt = 3;                       % порядок Тейлора
  nt = 0;                       % нижній порядок Тейлора
  %Nch = 6;                     % порядок Чебушева
  %nch = 0;                     % нижній порядок Чебушева
  p1=-1;                        % нижні межі проміжку апроксимації
  p2=(4/4)*pi*4;                % верхні межі проміжку апроксимації
  p=(p2-p1)/2+p1;               % середина проміжку апроксимації

  % отримаємо коефіцієнти полінома ортогонального базису
  disp('wait')
  P=FullBasis(nf,Nf,p1,p2,H);
  for k=1:n+1
      for j=1:n+1
          F=(@(x)Basis2Poly(x,P(k,1:end),P(j,1:end)));
          switch method
              case(1)
                  disp('wait...')
                  [Y1(k,j),Y_1(k,j)]=quad(F,p1,p2,tol); 
              case(2)
                  disp('wait... ...')
                  [Y2(k,j),Y_2(k,j)]=quadl(F,p1,p2,tol);
              case(3)
                  disp('wait... ... ...')
                  [Y3(k,j),Y_3(k,j)]=quadv(F,p1,p2,tol);
          end
          % перевірка похибки базису
      end
      F=(@(x)BasisPolyFunc(x,P(k,1:end)));
      switch method
          case (1)   
              disp(['...',num2str(n-k+1)]) 
              num2str(n-k+1)
              [E1(1,k),E1(2,k)]=quad(F,p1,p2,tol);
          case (2)
                    disp(['... ...',num2str(n-k+1)]) 
                    num2str(n-k+1)
                    [E2(1,k),E2(2,k)]=quadl(F,p1,p2,tol);
          case (3)
                    disp(['... ... ... ',num2str(n-k+1)]) 
                    num2str(n-k+1)
                    [E3(1,k),E3(2,k)]=quadv(F,p1,p2,tol);
          otherwise
              disp(['wait... ... ... ... ',num2str(n-k+1)]) 
              Y0(k)=ConvInt(P(k,1:end),P(k,1:end),p1,p2); %=1
              disp(['... ... ... ... ',num2str(n-k+1)])
              
              [E0(1,k),E0(2,k)]=quadgk(F,p1,p2,'RelTol',RelTol,'AbsTol',AbsTol,'MaxIntervalCount',discreteInt); 
              %'AbsTol' 'RelTol' 'Waypoints' 'MaxIntervalCount'
      end
  end
  
  options.TRANSA=true;
  switch method
      case (1)
          Det=det(Y1); % перевірка чи матриця добре обумовленна
          C=linsolve(Y1,(E1(1,1:end))',options); % поліноміальні коефіцієнти
          Mis1=(C'*Y1-E1(1,1:end))'  % перевірка похибки C1 (не ідеальна матриця)
      case (2)
          Det=det(Y2); % перевірка чи матриця добре обумовленна
          C=linsolve(Y2,(E2(1,1:end))',options); % поліноміальні коефіцієнти
          Mis2=(C'*Y2-E2(1,1:end))'  % перевірка похибки C1 (не ідеальна матриця)
      case (3)
          Det=det(Y3); % перевірка чи матриця добре обумовленна
          C=linsolve(Y3,(E3(1,1:end))',options); % поліноміальні коефіцієнти
          Mis3=(C'*Y3-E3(1,1:end))'  % перевірка похибки C1 (не ідеальна матриця)
      otherwise
          Det=1;
          C=E0(1,1:end)./Y0; % похибка матриці відсутня через метод обрахування без застосування матриці, але при умові що функції базсу повністю ортогональні
          %C=E(1,1:end); % формула для нормованого ортогонально базису
          %Mis0=0
  end
  if (Det<1e-3)
      disp('Матриця погано обумовленна! Визначник матриці меньше мілі одниці.')
  end
  
  x=linspace(p1,p2,discrete); % дискретизація змінної
  P(n+2,1:end)=0; % пустий поліном для запобігання конфлікту побудови графіка останьої базисної функції
  C(n+2)=0;       % додатковий зайвий коефіцієнт для запобіганю конфлікта побудови графіків
  H(n+2)=0;       % додатковий зайвий коефіцієнт для запобіганю конфлікта побудови графіків
  y=0;
  D(1:Nf+1)=0;
  figure('Name','Порівняльна характеристика апроксимації функції за формулами функціонального аналізу','NumberTitle','off')

  for r=1:n+1
      Taylor=ftaylor(@(x)func(x),p,nt+r-1); % формула Тейлора
      tp=eval(Taylor);
      y=C(r)*basisPoly(x,P(r,1:end))+y;
      Fplot(x,tp,y,C(r+1)*basisPoly(x,P(r+1,1:end)),discrete,r-1,H(r+1));
      pause()
      %close
      D(1:Nf+1)=P(r,1:end).*C(r)+D(1:Nf+1);
  end
  close
  C'
  
  syms x
  U=basisPoly(x,D)
  pretty(U);
