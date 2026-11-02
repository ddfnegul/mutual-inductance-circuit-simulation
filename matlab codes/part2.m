function part2
close all;
UL1_data=[20 50 100 150 200 250 280 300];
Mdata=[0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];
h=0.001;
t=0:h:30;
excitations={ones(1,length(t)),240*sin(t),120*sin(2*pi*50*t)};
titles={'e=1 V','e=240sin(t)','e=120sin(2\pift),f=50Hz'};
anames={'Lagrange','Spline','Poly n=3','Poly n=5'};
for ex=1:3
e=excitations{ex};
for am=1:4
y=zeros(3,length(t));
Marr=zeros(1,length(t));
Marr(1)=fM(0,am,UL1_data,Mdata);
for i=1:length(t)-1
em=(e(i)+e(i+1))/2;
[k1,Mv]=F(t(i),y(:,i),e(i),am,UL1_data,Mdata);
Marr(i+1)=Mv;
k2=F(t(i)+h/2,y(:,i)+h/2*k1,em,am,UL1_data,Mdata);
k3=F(t(i)+h/2,y(:,i)+h/2*k2,em,am,UL1_data,Mdata);
k4=F(t(i+1),y(:,i)+h*k3,e(i+1),am,UL1_data,Mdata);
y(:,i+1)=y(:,i)+h*(k1+2*k2+2*k3+k4)/6;
end
i1=y(1,:); i2=y(2,:); uc=y(3,:);
R2=10; uR2=R2*i2;
figure(5*(ex-1)+1); hold on; grid on; plot(t,i1); xlabel('t [s]'); ylabel('i_1 [A]'); title(['i_1 - ' titles{ex}]);
figure(5*(ex-1)+2); hold on; grid on; plot(t,i2); xlabel('t [s]'); ylabel('i_2 [A]'); title(['i_2 - ' titles{ex}]);
figure(5*(ex-1)+3); hold on; grid on; plot(t,uR2); xlabel('t [s]'); ylabel('u_{R2} [V]'); title(['u_R2 - ' titles{ex}]);
figure(5*(ex-1)+4); hold on; grid on; plot(t,uc); xlabel('t [s]'); ylabel('u_C [V]'); title(['u_C - ' titles{ex}]);
figure(5*(ex-1)+5); hold on; grid on; plot(t,Marr); xlabel('t [s]'); ylabel('M_n [H]'); title(['M_n - ' titles{ex}]);
end
for fig_off=1:5
figure(5*(ex-1)+fig_off); legend(anames);
end
end
end
function [dy,M]=F(t,y,e,am,UL1_data,Mdata)
R1=0.1; R2=10; C=0.5; L1=3; L2=5;
i1=y(1); i2=y(2); uc=y(3);
uL1=e-uc-R1*i1;
M=fM(uL1,am,UL1_data,Mdata);
dy=[1/(L1/M-M/L2)*(-R1/M*i1+R2/L2*i2-1/M*uc+1/M*e)
1/(M/L1-L2/M)*(-R1/L1*i1+R2/M*i2-1/L1*uc+1/L1*e)
1/C*i1];
end
function M=fM(u,am,UL1_data,Mdata)
u=abs(u);
if u>300
M=0.18;
return;
end
if am==1
M=mylagrange(UL1_data,Mdata,u);
elseif am==2
M=myspline(UL1_data,Mdata,u);
elseif am==3
M=mypolyapprox(UL1_data,Mdata,u,3);
else
M=mypolyapprox(UL1_data,Mdata,u,5);
end
end
function y=mylagrange(x_in,y_in,x)
y=0;
N=length(x_in);
for i=1:N
p=1;
for j=1:N
if i~=j
p=p*(x-x_in(j))/(x_in(i)-x_in(j));
end
end
y=y+y_in(i)*p;
end
end
function y=myspline(x_in,y_in,x)
n=length(x_in);
hv=diff(x_in);
dv=diff(y_in)./hv;
A=zeros(n-2,n-2);
bv=zeros(n-2,1);
for i=1:n-2
bv(i)=6*(dv(i+1)-dv(i));
A(i,i)=2*(hv(i)+hv(i+1));
if i>1; A(i,i-1)=hv(i); end
if i<n-2; A(i,i+1)=hv(i+1); end
end
Mv=[0;A\bv;0];
for k=1:n-1
if x>=x_in(k) && x<=x_in(k+1)
hk=hv(k);
a=(Mv(k+1)-Mv(k))/(6*hk);
b2=Mv(k)/2;
c2=(y_in(k+1)-y_in(k))/hk-hk*(2*Mv(k)+Mv(k+1))/6;
d2=y_in(k);
dx=x-x_in(k);
y=a*dx^3+b2*dx^2+c2*dx+d2;
return;
end
end
y=y_in(end);
end
function y=mypolyapprox(x_in,y_in,x,deg)
n=length(x_in);
xmin=min(x_in); xmax=max(x_in);
xn=(x_in-xmin)/(xmax-xmin);
A=zeros(n,deg+1);
for i=1:n
for j=1:deg+1
A(i,j)=xn(i)^(deg+1-j);
end
end
c=(A'*A)\(A'*y_in');
xq=(x-xmin)/(xmax-xmin);
y=0;
for j=1:deg+1
y=y+c(j)*xq^(deg+1-j);
end
end
