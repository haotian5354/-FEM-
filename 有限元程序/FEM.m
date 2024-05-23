%定义参数
a=1;            %椭圆长轴
E_M=1   ;       %基体模量
E_I=[0.001,1:10];     %夹杂模量
mu_M=0.3;       %基体泊松比
mu_I=0.2;       %夹杂泊松比
sigma_0=1;      %应力边条
b=[1/3,1,3];    %三种椭圆短轴
epsilon=2;     %epsilon为平面矩形宽度和椭圆短半轴的比值
t_0=1;          %平板的厚度
%% 有限元前处理
%[p,e,t]=preprocess(a,b(1),epsilon);  %椭圆夹杂
%网格剖分函数,p为网格节点，e为网格边缘，t为网格单元
%p的两列分别为节点的x、y坐标
%t的前三列为节点编号，最后一列为1则为方板单元，为2则为夹杂单元
[p,e,t]=preprocess_rec(a,b(1),epsilon);  %矩形夹杂
%% 有限元计算
K=zeros(length(p)*2,length(p)*2);
%初始化总体刚度矩阵,二维问题DOF为2*节点数
En=zeros(length(t),1);
En_0=zeros(length(t),1);
%单元应变能矩阵

for i=1:length(t)
    if t(4,i)==1
        mu=mu_M;
        E=E_M;
    elseif t(4,i)==2
        mu=mu_I;
        E=E_I(3);
    end
%根据单元类型确定材料参数，t(4,i)为1即基体，t(4,i)为2即夹杂
    K_i=element_stiffness_matrix(i,p,t,E,mu,t_0);
%生成i单元的单元刚度矩阵
    for j=1:3
        for k=1:3
        K(t(j,i)*2-1:t(j,i)*2,t(k,i)*2-1:t(k,i)*2)=K(t(j,i)*2-1:t(j,i)*2,t(k,i)*2-1:t(k,i)*2)+K_i(j*2-1:j*2,k*2-1:k*2);
        %将单刚按节点顺序组装进去，j表示行，k表示列
        end
    end
%单元刚度阵组装到总体刚度阵
end
%F=equivalent_nodal_load(p,e,sigma_0,t_0);  %椭圆等效结点荷载
F=equivalent_nodal_load_rectangular(p,e,sigma_0,t_0);  %矩形等效节点荷载
for i=1:length(e)
    %if (e(5,i)==5)||(e(5,i)==6)  %椭圆夹杂
        %满足任意一条即为y对称约束边界,u_x=0
    if (e(5,i)==8)||(e(5,i)==7)    %矩形夹杂
        K((2*e(1,i)-1),(2*e(1,i)-1))=1e5*K((2*e(1,i)-1),(2*e(1,i)-1));
        K((2*e(2,i)-1),(2*e(2,i)-1))=1e5*K((2*e(2,i)-1),(2*e(2,i)-1));
        F((2*e(1,i)-1),1)=0;
        F((2*e(2,i)-1),1)=0;
        %乘大数法
    %elseif (e(5,i)==3)||(e(5,i)==4) %椭圆夹杂
    elseif (e(5,i)==5)||(e(5,i)==6)   %矩形夹杂
        %满足任意一条即为x对称约束边界，u_y=0
        K((2*e(1,i)),(2*e(1,i)))=1e5*K((2*e(1,i)),(2*e(1,i)));
        K((2*e(2,i)),(2*e(2,i)))=1e5*K((2*e(2,i)),(2*e(2,i)));
        F(2*e(1,i),1)=0;
        F(2*e(2,i),1)=0;
    end
end
%约束条件的引入，消除奇异性
u=K\F;
%求解位移向量
%% 有限元后处理
u_node=zeros(size(p,2),1);
for i=1:length(u_node)
    u_node(i)=sqrt(u(2*i-1)^2+u(2*i)^2);
end
%节点位移magnitude
fill([p(1,t(1,:));p(1,t(2,:));p(1,t(3,:))],[p(2,t(1,:)); ...
    p(2,t(2,:));p(2,t(3,:))],[u_node(t(1,:))'; ...
    u_node(t(2,:))';u_node(t(3,:))'])
colorbar
colormap(jet)
axis equal

%dE_total_ellipse=strain_energy_ellipse(mu_M,E_M,mu_I,E_I(3),p,t,u,e);  %计算椭圆夹杂引起的应变能改变量
dE_total_rectangular=strain_energy_rectangular(mu_M,E_M,mu_I,E_I(3),p,t,u,e);  %计算矩形夹杂引起的应变能改变量




    



