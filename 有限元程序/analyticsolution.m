%定义参数
a=1;            %椭圆长轴
E_M=1   ;       %基体模量
E_I=10;          %夹杂模量
mu_M=0.3;       %基体泊松比
mu_I=0.2;       %夹杂泊松比
sigma_0=1;      %应力边条
b=[1/3,1,3];    %三种椭圆短轴
epsilon=10;     %epsilon为平面矩形宽度和椭圆短半轴的比值

lambda_I=E_I*mu_I/(1+mu_I)/(1-2*mu_I);
lambda_M=E_M*mu_M/(1+mu_M)/(1-2*mu_M);
%lame常数
G_M=E_M/(2*(1+mu_M));
G_I=E_I/(2*(1+mu_I));
%剪切模量
u_r=sigma_0*a*(1-mu_M)/G_M*(1-2*mu_I)/(1-2*mu_I+G_I/G_M);
delta_U=pi*sigma_0*a*u_r*(1-(1+mu_M)*E_I/(1+mu_I) ...
    *E_M-2*(mu_I-mu_M)*(1+mu_M)*E_I/(1+mu_I)/(1-2*mu_I)/E_M);
