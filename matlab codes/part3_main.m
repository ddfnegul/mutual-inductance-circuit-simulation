function part3_main
global UL1_tab M_tab
UL1_tab=[20 50 100 150 200 250 280 300];
M_tab=[0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];
excitations={@(t)ones(size(t)),'e=1V';@(t)240*sin(t),'e=240sin';
    @(t)210*sin(2*pi*5*t),'e=210sin5Hz'};
h_vals=[0.001 0.5];
M_types={'linear','nonlinear'};
fig_idx=1;
for ei=1:3
    e_func=excitations{ei,1};e_name=excitations{ei,2};
    for mi=1:2
        for hi=1:2
            h=h_vals(hi);
            [t,i1,i2]=simulate(e_func,h,M_types{mi});
            p=0.1*i1.^2+10*i2.^2;
            P_rect=sum(p(1:end-1))*h;
            P_trap=h*(p(1)/2+sum(p(2:end-1))+p(end)/2);
            fprintf('%s %s h=%g: rect=%.4f trap=%.4f\n',e_name,M_types{mi},h,P_rect,P_trap);
            figure(fig_idx);fig_idx=fig_idx+1;
            plot(t,p,'b');xlabel('t [s]');ylabel('R_1i_1^2+R_2i_2^2 [W]');grid on;
            title(sprintf('%s %s h=%g',e_name,M_types{mi},h));
        end
    end
end
end
function [t,i1,i2]=simulate(e_func,h,M_type)
global UL1_tab M_tab
t=0:h:30;e=e_func(t);y=zeros(3,length(t));
for i=1:length(t)-1
    emid=(e(i)+e(i+1))/2;
    k1=F(y(:,i),e(i),M_type);
    k2=F(y(:,i)+h/2*k1,emid,M_type);
    k3=F(y(:,i)+h/2*k2,emid,M_type);
    k4=F(y(:,i)+h*k3,e(i+1),M_type);
    y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
end
i1=y(1,:);i2=y(2,:);
end
function dy=F(y,e,M_type)
global UL1_tab M_tab
R1=0.1;R2=10;C=0.5;L1=3;L2=5;
i1=y(1);i2=y(2);uc=y(3);
if strcmp(M_type,'linear')
    M=0.8;
else
    M=fM(abs(e-uc-R1*i1));
end
D1=L1/M-M/L2;D2=M/L1-L2/M;
dy=[1/D1*(-R1/M*i1+R2/L2*i2-1/M*uc+1/M*e);
    1/D2*(-R1/L1*i1+R2/M*i2-1/L1*uc+1/L1*e);
    1/C*i1];
end
function M=fM(u)
global UL1_tab M_tab
if u>300;M=0.18;return;end
M=mylagrange(UL1_tab,M_tab,u);
end
function y=mylagrange(x_in,y_in,x)
y=0;N=length(x_in);
for i=1:N
    p=1;
    for j=1:N
        if i~=j;p=p*(x-x_in(j))/(x_in(i)-x_in(j));end
    end
    y=y+y_in(i)*p;
end
end
