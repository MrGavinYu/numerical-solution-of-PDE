function [x, area] = DevideX(height, h)
%DEVIDEX 划分空间网格
% 采取策略为每一层材料步长基本相同，如果材料长度无法整除建议的步长，
% 则将最后多余的长度与最后一个网格合并。
%
% 返回值：
% x 位置坐标
% area 标志着网格所属的材料层
%
% Gavin <www.bigbugs.cn>

x = [];
height = reshape(height, 1, length(height));
x(1) = 0;num = 1;deep = 0;
L = cumsum(height);
for i = 1:length(height)
    pos = deep+h(i):h(i):L(i); % 按照给定的步长对一层材料划分
    if(isempty(pos)) % 此时应该是步长大于材料长度
        error('步长大于材料长度')
    end
    if(pos(end) ~= L(i)) % 材料长度无法整除步长，则合并材料多余长度
        pos(end) = L(i);
    end
    x(num+1:num+length(pos)) = pos;
    num = num + length(pos);
    deep = deep + height(i);
end
n = length(x);
area = zeros(1, n);
L = [0 L];
% 从x(i-1)至x(i)为第area(i)层材料
area(1) = 1;
for i = 1:length(height)
    mask = x > L(i) & x <= L(i+1);
    area(mask) = i;
end
end
