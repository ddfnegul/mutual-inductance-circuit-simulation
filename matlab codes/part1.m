function part1
h = 0.001;
t = 0:h:30;
T = 3;
excitations = {ones(1,length(t)), 240*sin(t), 210*sin(2*pi*5*t), ...
               120*sin(2*pi*50*t), 120*(mod(t,T)<T/2)};
mnames = {'Euler','Heun','RK4'};
titles = {'e=1 V','e=240sin(t)','e=210sin(2\pift), f=5 Hz', ...
          'e=120sin(2\pift), f=50 Hz','e=120 square, T=3 s'};
for ex = 1:5
    e = excitations{ex};
    lbl_i = {};
    lbl_v = {};
    figure(2*ex-1); hold on; grid on;
    xlabel('t [s]'); ylabel('i [A]');
    title(['Currents - ' titles{ex}]);
    figure(2*ex);   hold on; grid on;
    xlabel('t [s]'); ylabel('u [V]');
    title(['Voltages - ' titles{ex}]);
    for method = 1:3
        y = zeros(3,length(t));
        for i = 1:length(t)-1
            ei = e(i); en = e(i+1); em = (ei+en)/2;
            k1 = F(t(i),y(:,i),ei);
            if method == 1
                y(:,i+1) = y(:,i) + h*k1;
            elseif method == 2
                k2 = F(t(i+1),y(:,i)+h*k1,en);
                y(:,i+1) = y(:,i) + h*(k1+k2)/2;
            else
                k2 = F(t(i)+h/2,y(:,i)+h/2*k1,em);
                k3 = F(t(i)+h/2,y(:,i)+h/2*k2,em);
                k4 = F(t(i+1),y(:,i)+h*k3,en);
                y(:,i+1) = y(:,i) + h*(k1+2*k2+2*k3+k4)/6;
            end
        end
        figure(2*ex-1);
        plot(t,y(1,:)); plot(t,y(2,:));
        lbl_i{end+1} = ['i_1 ' mnames{method}];
        lbl_i{end+1} = ['i_2 ' mnames{method}];
        figure(2*ex);
        plot(t,y(3,:)); plot(t,e);
        lbl_v{end+1} = ['u_C ' mnames{method}];
        lbl_v{end+1} = ['e ' mnames{method}];
    end
    figure(2*ex-1); legend(lbl_i);
    figure(2*ex);   legend(lbl_v);
end
end

function dy = F(t,y,e)
R1=0.1; R2=10; C=0.5; L1=3; L2=5; M=0.8;
i1=y(1); i2=y(2); uc=y(3);
dy = [1/(L1/M-M/L2)*(-R1/M*i1+R2/L2*i2-1/M*uc+1/M*e)
      1/(M/L1-L2/M)*(-R1/L1*i1+R2/M*i2-1/L1*uc+1/L1*e)
      1/C*i1];
end
