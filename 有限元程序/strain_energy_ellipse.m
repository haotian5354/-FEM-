function E_total=strain_energy_ellipse(mu_M,E_M,mu_I,E_I,p,t,u,e)
%p的两列分别为节点的x、y坐标
%t的前三列为节点编号，最后一列为1则为方板单元，为2则为夹杂单元
%e的前两行为边缘单元中的边缘上的两个点的编号，第5行判断是什么边缘（共7条边缘），
% 第6行判断为是否是夹杂单元，第7行判断边缘是否为曲线（1为曲线）
sigma=element_stress(mu_M,E_M,mu_I,E_I,p,t,u);
%求单元应力矩阵
E_total=0;
En=zeros(size(e,2),1);
for i=1:size(e,2)
    if e(7,i)==1
        %所有曲线进行积分
        n=[p(2,e(2,i))-p(2,e(1,i)),p(1,e(1,i))-p(1,e(2,i))]';
        n=n/norm(n);
        %边界处法向向量
        En(i)=1/2*(((1-(1+mu_M)*E_I/(1+mu_I)*E_M)*([1,0;0,1]*n)'-(mu_I-mu_M)*(1+mu_M)*E_I ...
            /(1+mu_I)/(1-2*mu_I)/E_M*2*n'))...
        *1/2*([u(2*e(1,i)-1),u(2*e(1,i))]'+[u(2*e(2,i)-1),u(2*e(2,i))]')*norm([p(1,e(1,i)),p(2,e(1,i))]-[p(1,e(2,i)),p(2,e(2,i))]);
    elseif e(7,i)~=1
        En(i)=0;
    end
    E_total=E_total+En(i);
end
E_total=4*E_total;
end
