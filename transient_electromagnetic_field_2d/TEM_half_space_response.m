function E=TEM_half_space_response(x,z,t,sigma,I)
% 计算二维均匀半空间中，一条无线长导线在电流I作用下，产生的瞬变电场
% x，z：坐标位置 z>0
% t：时间 t>0
% sigma: 导电率
% I: 电流强度
% E: 返回的电场强度值
% 参考《辛会翠硕士论文》书中2-28公式
if t<=0
    errordlg('t must be positive.');
    E=inf;
    return;
end
mu=4*pi*1e-7;
T=4*t/(mu*sigma);
if z<0
    errordlg('Z must be non-negative.');
    E=inf;
    return;
elseif z==0
    if x==0
        E=I/(pi*sigma*T);
    else
        E=I/(pi*sigma*x*x)*(1-exp(-x^2/T)); % z=0的情况
    end
else 
    %z>0的情况
    a3=37/84;
    a5=1/7;
    a7=13/105;
    b2=31/28;
    b4=43/70;
    b6=17/105;
    b8=26/105;

    RR=x^2+z^2;
    u=x/(sqrt(T));
    F=(u + a3*u^3 + a5*u^5 + a7*u^7)/(1 + b2*u^2 + b4*u^4 + b6*u^6 + b8*u^8);

    E=I/(pi*sigma)*(((z^2-x^2)/RR + 2*z^2/T)*exp(-RR/T)/RR...
        - 2*z*exp(-z^2/T)/(sqrt(pi)*RR)*(1/sqrt(T) - 2*x*F*(1/T + 1/RR)))...
        + I/(pi*sigma)*(x^2 - z^2)/(RR*RR)*erfc(z/sqrt(T));
end