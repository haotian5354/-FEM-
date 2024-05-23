%�������
a=1;            %��Բ����
E_M=1   ;       %����ģ��
E_I=[0.001,1:10];     %����ģ��
mu_M=0.3;       %���岴�ɱ�
mu_I=0.2;       %���Ӳ��ɱ�
sigma_0=1;      %Ӧ������
b=[1/3,1,3];    %������Բ����
epsilon=2;     %epsilonΪƽ����ο�Ⱥ���Բ�̰���ı�ֵ
t_0=1;          %ƽ��ĺ��
%% ����Ԫǰ����
%[p,e,t]=preprocess(a,b(1),epsilon);  %��Բ����
%�����ʷֺ���,pΪ����ڵ㣬eΪ�����Ե��tΪ����Ԫ
%p�����зֱ�Ϊ�ڵ��x��y����
%t��ǰ����Ϊ�ڵ��ţ����һ��Ϊ1��Ϊ���嵥Ԫ��Ϊ2��Ϊ���ӵ�Ԫ
[p,e,t]=preprocess_rec(a,b(1),epsilon);  %���μ���
%% ����Ԫ����
K=zeros(length(p)*2,length(p)*2);
%��ʼ������նȾ���,��ά����DOFΪ2*�ڵ���
En=zeros(length(t),1);
En_0=zeros(length(t),1);
%��ԪӦ���ܾ���

for i=1:length(t)
    if t(4,i)==1
        mu=mu_M;
        E=E_M;
    elseif t(4,i)==2
        mu=mu_I;
        E=E_I(3);
    end
%���ݵ�Ԫ����ȷ�����ϲ�����t(4,i)Ϊ1�����壬t(4,i)Ϊ2������
    K_i=element_stiffness_matrix(i,p,t,E,mu,t_0);
%����i��Ԫ�ĵ�Ԫ�նȾ���
    for j=1:3
        for k=1:3
        K(t(j,i)*2-1:t(j,i)*2,t(k,i)*2-1:t(k,i)*2)=K(t(j,i)*2-1:t(j,i)*2,t(k,i)*2-1:t(k,i)*2)+K_i(j*2-1:j*2,k*2-1:k*2);
        %�����հ��ڵ�˳����װ��ȥ��j��ʾ�У�k��ʾ��
        end
    end
%��Ԫ�ն�����װ������ն���
end
%F=equivalent_nodal_load(p,e,sigma_0,t_0);  %��Բ��Ч������
F=equivalent_nodal_load_rectangular(p,e,sigma_0,t_0);  %���ε�Ч�ڵ����
for i=1:length(e)
    %if (e(5,i)==5)||(e(5,i)==6)  %��Բ����
        %��������һ����Ϊy�Գ�Լ���߽�,u_x=0
    if (e(5,i)==8)||(e(5,i)==7)    %���μ���
        K((2*e(1,i)-1),(2*e(1,i)-1))=1e5*K((2*e(1,i)-1),(2*e(1,i)-1));
        K((2*e(2,i)-1),(2*e(2,i)-1))=1e5*K((2*e(2,i)-1),(2*e(2,i)-1));
        F((2*e(1,i)-1),1)=0;
        F((2*e(2,i)-1),1)=0;
        %�˴�����
    %elseif (e(5,i)==3)||(e(5,i)==4) %��Բ����
    elseif (e(5,i)==5)||(e(5,i)==6)   %���μ���
        %��������һ����Ϊx�Գ�Լ���߽磬u_y=0
        K((2*e(1,i)),(2*e(1,i)))=1e5*K((2*e(1,i)),(2*e(1,i)));
        K((2*e(2,i)),(2*e(2,i)))=1e5*K((2*e(2,i)),(2*e(2,i)));
        F(2*e(1,i),1)=0;
        F(2*e(2,i),1)=0;
    end
end
%Լ�����������룬����������
u=K\F;
%���λ������
%% ����Ԫ����
u_node=zeros(size(p,2),1);
for i=1:length(u_node)
    u_node(i)=sqrt(u(2*i-1)^2+u(2*i)^2);
end
%�ڵ�λ��magnitude
fill([p(1,t(1,:));p(1,t(2,:));p(1,t(3,:))],[p(2,t(1,:)); ...
    p(2,t(2,:));p(2,t(3,:))],[u_node(t(1,:))'; ...
    u_node(t(2,:))';u_node(t(3,:))'])
colorbar
colormap(jet)
axis equal

%dE_total_ellipse=strain_energy_ellipse(mu_M,E_M,mu_I,E_I(3),p,t,u,e);  %������Բ���������Ӧ���ܸı���
dE_total_rectangular=strain_energy_rectangular(mu_M,E_M,mu_I,E_I(3),p,t,u,e);  %������μ��������Ӧ���ܸı���




    



