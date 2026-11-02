function project_2026L_part_4
freq_v=0.05:0.05:5;
P_v=zeros(1,length(freq_v));
for k=1:length(freq_v)
P_v(k)=sim(freq_v(k));
end
figure(1)
plot(freq_v,P_v,'b',freq_v,406*ones(size(freq_v)),'r--')
xlabel('f [Hz]');ylabel('P [J]');title('P(f) for e(t)=100sin(2\pift)');grid on
legend('P(f)','P = 406 J')
G=@(f)(sim(f)-406);
ft=0.60;Gt=G(ft);
dfs=[0.10,0.05,0.01,0.005,0.001,0.0005];dp=NaN;
fprintf('df\t\tdP/df\t\tChange[%%]\n');
for d=1:length(dfs)
df=dfs(d);deriv=(G(ft+df)-Gt)/df;
if ~isnan(dp)
fprintf('%.4f\t\t%.2f\t\t%.2f\n',df,deriv,abs(deriv-dp)/abs(dp)*100);
else
fprintf('%.4f\t\t%.2f\t\t-\n',df,deriv);
end
dp=deriv;
end
df_ch=0.001;
fprintf('\nChosen df = %.4f\n',df_ch);
[fb,ib,cb]=bisect(G,0.5,0.7);
fprintf('Bisection: f=%.5f P=%.2f iters=%d calls=%d\n',fb,sim(fb),ib,cb);
[fs,is,cs]=secant(G,0.5,0.7);
fprintf('Secant: f=%.5f P=%.2f iters=%d calls=%d\n',fs,sim(fs),is,cs);
[fq,iq,cq]=quasi_newton(G,0.60,df_ch);
fprintf('QuasiNewton: f=%.5f P=%.2f iters=%d calls=%d\n',fq,sim(fq),iq,cq);
end

function [fs,it,nc]=bisect(G,a,b)
tol=5e-1;nc=1;Ga=G(a);
for i=1:200
c=(a+b)/2;Gc=G(c);nc=nc+1;
if abs(Gc)<tol||(b-a)/2<1e-4;it=i;fs=c;return;end
if Ga*Gc<0;b=c;else;a=c;Ga=Gc;end
end
it=200;fs=(a+b)/2;
end

function [fs,it,nc]=secant(G,f0,f1)
tol=5e-1;G0=G(f0);G1=G(f1);nc=2;
for i=1:200
if abs(G1-G0)<1e-14;break;end
f2=f1-G1*(f1-f0)/(G1-G0);G2=G(f2);nc=nc+1;
if abs(G2)<tol;it=i;fs=f2;return;end
f0=f1;G0=G1;f1=f2;G1=G2;
end
it=i;fs=f1;
end

function [fs,it,nc]=quasi_newton(G,f0,df)
tol=5e-1;nc=0;
for i=1:100
Gf=G(f0);Gfd=G(f0+df);nc=nc+2;
f1=f0-Gf/((Gfd-Gf)/df);
Gf1=G(f1);nc=nc+1;
if abs(Gf1)<tol;it=i;fs=f1;return;end
f0=f1;
end
it=100;fs=f0;
end

function P=sim(freq)
UL1=[20 50 100 150 200 250 280 300];
Md=[0.46 0.64 0.78 0.68 0.44 0.23 0.18 0.18];
h=0.001;t=0:h:30;N=length(t);
e=100.*sin(2.*pi.*freq.*t);
y=[0;0;0];i1v=zeros(1,N);i2v=zeros(1,N);
for i=1:N-1
em=(e(i)+e(i+1))/2;
k1=F(y,e(i),UL1,Md);
k2=F(y+(h/2)*k1,em,UL1,Md);
k3=F(y+(h/2)*k2,em,UL1,Md);
k4=F(y+h*k3,e(i+1),UL1,Md);
y=y+h*(k1+2*k2+2*k3+k4)/6;
i1v(i+1)=y(1);i2v(i+1)=y(2);
end
pw=0.1*i1v.^2+10.*i2v.^2;
P=h*(pw(1)/2+sum(pw(2:end-1))+pw(end)/2);
end

function dy=F(y,e,UL1,Md)
R1=0.1;R2=10;C=0.5;L1=3;L2=5;
i1=y(1);i2=y(2);uc=y(3);
M=fM(e-uc-R1*i1,UL1,Md);
D1=L1/M-M/L2;D2=M/L1-L2/M;
dy=[1/D1*(-R1/M*i1+R2/L2*i2-uc/M+e/M);
    1/D2*(-R1/L1*i1+R2/M*i2-uc/L1+e/L1);
    i1/C];
end

function M=fM(u,UL1,Md)
u=abs(u);
if u>300;M=0.18;return;end
M=mylagrange(UL1,Md,u);
end

function y=mylagrange(xin,yin,x)
y=0;N=length(xin);
for i=1:N
p=1;
for j=1:N
if i~=j;p=p*(x-xin(j))/(xin(i)-xin(j));end
end
y=y+yin(i)*p;
end
end
