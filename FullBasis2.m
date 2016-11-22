clear;      % звільнити пам'ять MATLAB від усіх змінних
close all;  % закрити всі графічні вінка
y1=0;
y2=3;
p1=0;
p2=1;
options.TRANSA=true;
P(1,y2+1)=0;
P(1,y2+1-y1)=1;
for n=2:y2-y1+1
    P(n,y2+1-y1)=1;
    % обнулення матриці A та вектора B (з запасом)
    for k=1:n-1
        for m=1:n-1
            A(k,m)=0;
        end
        B(k)=0;
    end
    for k=1:n-1
        for v=1:n-1
            for m=1:n-1
                A(m,k)=P(k,y2+2-y1-v).*delta(p1,p2,m+v+y1^2)+A(m,k);
            end
            B(k)=-(P(k,y2+2-y1-v).*delta(p1,p2,v=y1^2))+B(k);
        end
    end
    % перевернемо вектор рішення та запишемо його коефіцієнти в поліном
    %P(n,y2-y1-n+2:y2-y1)=linsolve(A,B',options)
    p(1:n-1)=linsolve(A,B',options);
    for z=0:n-2
        P(n,y2-y1-z)=p(z+1)
    end
end
