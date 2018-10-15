function [x, area] = DevideX(height, h)
%DEVIDEX ���ֿռ�����
% ��ȡ����Ϊÿһ����ϲ���������ͬ��������ϳ����޷���������Ĳ�����
% ��������ĳ��������һ������ϲ���
%
% ����ֵ��
% x λ������
% area ��־�����������Ĳ��ϲ�
%
% Gavin <www.bigbugs.cn>

x = [];
height = reshape(height, 1, length(height));
x(1) = 0;num = 1;deep = 0;
L = cumsum(height);
for i = 1:length(height)
    pos = deep+h(i):h(i):L(i); % ���ո����Ĳ�����һ����ϻ���
    if(isempty(pos)) % ��ʱӦ���ǲ������ڲ��ϳ���
        error('�������ڲ��ϳ���')
    end
    if(pos(end) ~= L(i)) % ���ϳ����޷�������������ϲ����϶��೤��
        pos(end) = L(i);
    end
    x(num+1:num+length(pos)) = pos;
    num = num + length(pos);
    deep = deep + height(i);
end
n = length(x);
area = zeros(1, n);
L = [0 L];
% ��x(i-1)��x(i)Ϊ��area(i)�����
area(1) = 1;
for i = 1:length(height)
    mask = x > L(i) & x <= L(i+1);
    area(mask) = i;
end
end
