function [p,e,t]=preprocess_rec(a,b,epsilon)
ellipse=[3,4,0,sqrt(pi*a*b/4),sqrt(pi*a*b/4),0,0,0,sqrt(pi*a*b/4),sqrt(pi*a*b/4)]';
%rectangular=[3,4,0,epsilon*a,epsilon*a,0,0,0,epsilon*b,epsilon*b]';
%列向量含义：3代表矩形，4代表4个点，0,epsilon*a,epsilon*a,0代表四个点的横坐标，0,0,epsilon*b,epsilon*b代表四个点的纵坐标
rectangular=[3,4,0,5,5,0,0,0,epsilon*a/2,epsilon*a/2]';
g_m=[ellipse,rectangular];
%将几何列向量组合构建几何矩阵
sf='rectangular-ellipse+rectangular&ellipse';
%set formula,此式中定义几何为矩形板减去椭圆夹杂再加上二者的交集
ns = char('ellipse','rectangular');
%ns：Name-space matrix
ns = ns';
g = decsg(g_m,sf,ns);
%Decompose constructive solid 2-D geometry into minimal regions,  g_m- Geometry description matrix，sf - Set formula，ns - Name-space matrix
%g：Decomposed geometry matrix
model = createpde;
geometryFromEdges(model,g);
figure(1)
pdegplot(model,'EdgeLabels','on','FaceLabels','on')
%画几何图形
[p,e,t] = initmesh(g);
%网格剖分函数,p为网格节点，e为网格边缘，t为网格单元
%p的两列分别为节点的x、y坐标
%t的前三列为节点编号，最后一列为1则为方板单元，为2则为夹杂单元
%e的前两行为边缘单元中的边缘上的两个点的编号，第5行判断是什么边缘（共7条边缘），第6行判断为是否是夹杂单元，第7行判断边缘是否为曲线（1为曲线）
figure(2)
pdemesh(p,e,t);
%绘制网格图
 [p1,e1,t1]=refinemesh(g,p,e,t);
 [p2,e2,t2]=refinemesh(g,p1,e1,t1);
 [p,e,t]=refinemesh(g,p2,e2,t2);
 %网格剖分连续细化三次并获得细化后的网格信息p2、e2、t2
 figure(3)
pdemesh(p,e,t);
%绘制细化后的网格图
end