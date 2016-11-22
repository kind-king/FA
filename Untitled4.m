clear;      % звільнити пам'ять MATLAB від усіх змінних
close all;  % закрити всі графічні вінка
y1=1;
y2=2;
p1=0;
p2=1;
options.TRANSA=true;
%U=zeros(1, y2-y1+1);
%P(y1+1,y2+1-y1:y2+1)=1;
%P(y1+1,y2+1-y1)=1;
%P(y1+1:y2+1,y2+1-y1)=1;
%P(1:y2-y1+1,y2+1-y1)=1;
% задаємо довжину поліному
P(1,y2+1)=0;
% ведемо додаткову змінну для зручності оперування
y=y2-y1;
% зазначаємо першу константу від якої можна буде відштовхуватися у рішені 
% системи рівнянь
P(1,y+1)=1;
for n=2:y+1
    %P(n,1:end-n+1)=1;%a
    %pc=conv(P(n,1:end),P(1:n-1,1:end));
    %pi=polyint(pc);
    % обнулення матриці A та вектора B (з запасом)
    for k=1:n-1
        for m=1:n-1
            A(k,m)=0;
        end
        B(k)=0;
    end
    %A=0;
    %B=0;
    % зазначаємо наступні константи...
    P(n,y+1)=1;
    for k=1:n-1
        for v=1:k
            for m=1:n-1
                % поміняємо місцями індекси матрці і таким чином транспортуємо її
                A(m,k)=P(k,y+2-v).*delta(p1,p2,m+v+y1.*2)+A(m,k);
                %A(k,m)=P(n-1,y2+2-y1-v).*delta(p1,p2,m+v)+A(k,m);
                %B(k)=-(P(n-1,y2+2-y1-v).*delta(p1,p2,v)+B(k))./(n-1);
            end
            B(k)=-(P(n,y+1).*P(k,y+2-v).*delta(p1,p2,v+y1.*2))+B(k);
        end
    end
    % Відладка
    A
    B
    %zeros(1, n+2);
    %B=zeros(1, n+2);
    %delta(p1,p2,);
    %B(1:n+2)=-delta(p1,p2,1);
    % запишемо коефіцієнти в поліном
    %P(n,y2-y1-n+2:y2-y1)=linsolve(A,B',options);
    % перевернемо вектор рішення та запишемо його коефіцієнти в поліном
    p(1:n-1)=linsolve(A,B',options)
    for z=0:n-2
        P(n,y-z)=p(z+1)
    end
    %for k=0:n
    %    U(y2+1-n)=P(n,n+1-k)*x^k+U(y2+1-n);
    %end
end

%y=0;
  %for n=1:Nf
  %    y=C(Nf+1-n)*basis(x,n)+y;
  %end
%function U=Basis2(x)
%    U=basis(x,kf).*basis(x,jf);
%    U=quad(@Basis2, p1,p2);
%end
%function U=BasisFun(x)
%global kf
%U=func(x).*basis(x,kf);
%end
%Y(kf,jf)=quad(@Basis2, p1,p2);% інтеграл базисних функцій
%E(kf)=quad(@BasisFun, p1,p2); % інтеграл функції на базисну функцію
%U=basis(x,kf);