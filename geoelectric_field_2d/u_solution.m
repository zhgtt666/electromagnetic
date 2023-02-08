function [u] = u_solution(x,z,xs,zs,sigma1,I)
% 计算二维均匀半空间中，一条无线长导线在电流I作用下，产生的瞬变电场
% x，z：坐标位置
% xs,zs: 线源位置
% sigma1: 导电率
% u: 返回的电位
% I: 电流强度
% 参考《稳定地电场三维有限差分正演模拟》中3c公式
% 注意：根据当前的计算，发现解析解的值过大，因此，目前在计算边界条件使用零值。
if x==xs && z==zs
  u=0;
else
 u= I/(2*pi*sigma1*sqrt((zs-z)^2+(xs-x)^2));
end
end

