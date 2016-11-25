clear;      % �������� ���'��� MATLAB �� ��� ������
close all;  % ������� �� ������� ����

Nf = 10;                      % ������� ��������������� ������
nf = 0;                       % ����� ������� ��������������� ������
n=Nf-nf;
% ��������� ��� ������
         w=1:n+1;
d=10.^(n+4-w);
H(1:n+1)=1; %d
p1=-1;                        % ���� ��� �������� ������������
p2=1;                         % ����� ��� �������� ������������
p=(p2-p1)/2+p1;               % �������� ������� ������������
  
r=0:4;
P=FullBasis(nf,Nf,p1,p2,H);

%basisPoly(x,P(r,1:end));

syms x y
%x=linspace(-1.5,1.5,100);

subplot(2,1,1);
fplot(chebyshevT(0, x))
axis([-1.5 1.5 -2 2])
grid on

ylabel('T_n(x)')
legend('T_0(x)', 'T_1(x)', 'T_2(x)', 'T_3(x)', 'T_4(x)', 'Location', 'Best')
title('Chebyshev polynomials of the first kind')

subplot(2,1,2);
fplot(basisPoly(x,P(r,1:end)))
axis([-1.5 1.5 -2 2])
grid on

ylabel('F_n(x)')
legend('F_0(x)', 'F_1(x)', 'F_2(x)', 'F_3(x)', 'F_4(x)', 'Location', 'Best')
title('Functional Analiz polynomials')

pause()
close