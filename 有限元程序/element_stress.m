function sigma=element_stress(mu_M,E_M,mu_I,E_I,p,t,u)
sigma=zeros(3,size(t,2));
for i=1:size(t,2)
    if t(4,i)==1
        mu=mu_M;
        E=E_M;
    elseif t(4,i)==2
        mu=mu_I;
        E=E_I;
    end
    %根据单元类型确定材料参数，t(4,i)为1即基体，t(4,i)为2即夹杂
    A=1/2*det([1,p(1,t(1,i)),p(2,t(1,i));1,p(1,t(2,i)),p(2,t(2,i));1,p(1,t(3,i)),p(2,t(3,i))]);
    %A为三角形单元的面积
    B=zeros(3,6);
    %应变-位移矩阵初始化
    rotate=[1,2,3,1,2,3];
    %轮换辅助矩阵
    for n=1:3
        B(1,2*n-1)=1/(A*2)*(p(2,t(rotate(n+1),i))-p(2,t(rotate(n+2),i)));
        B(1,2*n)=0;
        B(2,2*n-1)=0;
        B(2,2*n)=1/(A*2)*(-p(1,t(rotate(n+1),i))+p(1,t(rotate(n+2),i)));
        B(3,2*n-1)=B(2,2*n);
        B(3,2*n)=B(1,2*n-1);
    end
    %应变-位移矩阵，已进行数值验证，编写正确
    E_0=E/(1-mu^2);
    mu_0=mu/(1-mu);
    D_0=E_0/(1-mu_0^2);
    D=D_0*[1,mu_0,0;mu_0,1,0;0,0,(1-mu_0)/2];
    S=D*B;
    sigma(:,i)=S*[u(2*t(1,i)-1);u(2*t(1,i));u(2*t(2,i)-1);u(2*t(2,i));u(2*t(3,i)-1);u(2*t(3,i))];
    %计算应力矩阵，从上到下分别是sigma_x,sigma_y,sigma_xy
end
end