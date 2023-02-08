TEM - dt=1e-7s.gif :  由于初始时间dt=1e-7s过大，出现问题
TEM - dt=1e-8s.gif : 初始时间dt=1e-8s 改小后，消除问题
TEM - Riemann boundary.gif ： dt=1e-7s，但是顶边界用Riemann边界，说明上面TEM - dt=1e-7s.gif 中的问题，是上边界采用从地表往空气中延拓一个网格引起的。
TEM - 第二层未延展.gif ：在利用傅里叶变换延拓时，对其两端进行补零延拓。