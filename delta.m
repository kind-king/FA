function d=delta(a,b,n)
% �����������!
%d=( b.^(n)-a.^(n) )./n; ������������� �������, �� ���������� �� ��������...
% � ����������� a=-x �� b=x ��� ������� ������� �������������� ������....
% ����� ����� ���� ������� ���������� ���� ������� ���������� �� ��
% ������� (a-b) �� ��������� ������ ����������� �������� ��� ���� ����
% ��������: a^n-b^n=( a-b )*( a^(n-1)+a^(n-2)*b^1+a^(n-3)*b^(2)+...
% ...+a^2*b^(n-3)+a^1*b^(n-2)+b^(n-1) )
Sum=0;
for k=0:n-1
    Sum=b.^(n-1-k).*a.^k+Sum;
end
d=Sum./n;
end
