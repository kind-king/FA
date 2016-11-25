function P=FullBasis(y1,y2,p1,p2,H)
options.TRANSA=true;
%U=zeros(1, y2-y1+1);
%P(y1+1,y2+1-y1:y2+1)=1;
%P(y1+1,y2+1-y1)=1;
%P(y1+1:y2+1,y2+1-y1)=1;
%P(1:y2-y1+1,y2+1-y1)=1;
if (p2==-p1)
n0=1;
% ������ ������ ���� :
% 0 0 0 0 1
% 0 0 1 0 1
% 1 0 1 0 1

P(1,y2+1)=0;
y=y2-y1;
P(1,y+1)=H(1);

for n=2+n0:1+n0:y+1
    % ��������� ������� A �� ������� B 
    clear A B
    A((n-1)/(1+n0),(n-1)/(1+n0))=0;
    B((n-1)/(1+n0))=0;
    % ��������� �������� ���������...
    P(n,y+1)=H(n);
    for k=1:1+n0:n-1
        for v=1:1+0:k
            for m=1+n0:1+n0:n-1
                % �������� ������ ������� ������ � ����� ����� ������������ ��
                A(m/(1+n0),(k+n0)/(1+n0))=P(k,y+2-v).*intPoly(p1,p2,m+v+y1.*2)+A(m/(1+n0),(k+n0)/(1+n0))
                %A(k,m)=P(n-1,y2+2-y1-v).*delta(p1,p2,m+v)+A(k,m);
                %B(k)=-(P(n-1,y2+2-y1-v).*delta(p1,p2,v)+B(k))./(n-1);
            end
            B((k+n0)/(1+n0))=-(P(n,y+1).*P(k,y+2-v).*intPoly(p1,p2,v+y1.*2))+B((k+n0)/(1+n0))
        end
    end
    % ����������� ������ ������ �� �������� ���� ����������� � ������
    p(1+n0:1+n0:n-1)=linsolve(A,B',options);
    %linsolve(A,B',options)
    for z=0+n0:n-2
        P(n,y-z)=p(z+1);
    end
    '������� ���� �� �� ������ � ������� �����';
end

else
% ������ ������� �������
P(1,y2+1)=0;
% ������ ��������� ����� ��� �������� ����������
y=y2-y1;
% ��������� ����� ��������� �� ��� ����� ���� �������������� � ����� 
% ������� ������
P(1,y+1)=H(1)
%u(1)=H(1)*intPoly(p1,p2,1+y1.*2); % ���������� ����������
for n=2:y+1
    %n % ��������� ������������� (������� ����� �� �� ������� ��������)
    %P(n,1:end-n+1)=1;%a
    %pc=conv(P(n,1:end),P(1:n-1,1:end));
    %pi=polyint(pc);
    % ��������� ������� A �� ������� B (� �������)
    clear A B
    A(n-1,n-1)=0;
    B(n-1)=0;
    %for k=1:n-1
    %    for m=1:n-1
    %        A(k,m)=0;
    %    end
    %    B(k)=0;
    %end
    %A=0;
    %B=0;
    % ��������� �������� ���������...
    P(n,y+1)=H(n)
    for k=1:n-1
        for v=1:k
            for m=1:n-1
                % �������� ������ ������� ������ � ����� ����� ������������ ��
                A(m,k)=P(k,y+2-v).*intPoly(p1,p2,m+v+y1.*2)+A(m,k);
                %A(k,m)=P(n-1,y2+2-y1-v).*delta(p1,p2,m+v)+A(k,m);
                %B(k)=-(P(n-1,y2+2-y1-v).*delta(p1,p2,v)+B(k))./(n-1);
            end
            B(k)=-(P(n,y+1).*P(k,y+2-v).*intPoly(p1,p2,v+y1.*2))+B(k);
        end
    end
    
    %for k=0:y
    %Y0(y+1-k)=ConvInt(P(k+1,1:end),p1,p2);
    
    %function T=ConvInt(P,p1,p2)
    %T=0;
    %P=conv(P(1:end),P(1:end));
    %for k=0:length(P)-1
    %    T=P(length(P)-k).*intPoly(p1,p2,k+1)+T;
    %    %U=P(length(P)-n).*x.^n+U;
    %end
    
    % ³������
    %A
    %B
    %zeros(1, n+2);
    %B=zeros(1, n+2);
    %delta(p1,p2,);
    %B(1:n+2)=-delta(p1,p2,1);
    % �������� ����������� � ������
    %P(n,y2-y1-n+2:y2-y1)=linsolve(A,B',options);
    % ����������� ������ ������ �� �������� ���� ����������� � ������
    p(1:n-1)=linsolve(A,B',options);
    %
    %u=1./(sqrt( ConvInt(P(k+1,1:end),p1,p2) ))
    %
    for z=0:n-2
        P(n,y-z)=p(z+1);
    end
    %u(n)=(sqrt( ConvInt(P(n,1:end),p1,p2) )); % ���������� ����������
        %for k=0:n
    %    U(y2+1-n)=P(n,n+1-k)*x^k+U(y2+1-n);
    %end
end

%for n=1:y+1
%    P(n,1:end)=P(n,1:end)./u(n); % ����������
%end

end
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
%Y(kf,jf)=quad(@Basis2, p1,p2);% �������� �������� �������
%E(kf)=quad(@BasisFun, p1,p2); % �������� ������� �� ������� �������
%U=basis(x,kf);