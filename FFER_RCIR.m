clear;
data = readtable('FFER.xlsx');
X=data.rate;
N = length(X);T = N;dt = T/N;
dR=zeros(N-1,1);dL=zeros(N-1,1);
Y = zeros(N-1,1);Z = zeros(N-1,1);V = zeros(N-1,1);U = zeros(N-1,1);
for i=1:N-1
    Y(i)=X(i+1)/X(i);
    Z(i)=1/X(i);
    U(i)=X(i+1)-X(i);
end
A = sum(X(1:N-1));B = sum(Y);C = sum(Z);D = sum(U);
a=((N-1)*sum(X(2:N))-B*A)/(dt*(N-1)^2-dt*A*C);
b=((N-1)^2-(N-1)*B+D*C)/(dt*(N-1)^2-dt*A*C);  
for i=1:N-1
    V(i)=(X(i+1)-X(i)-(a-b*X(i))*dt)^2/X(i);
end
sigma = sqrt(1/(N-1)/dt*sum(V));
m=0.05;n=0.25;
Xt=zeros(1,N);
rng(4291489);
dw=sqrt(dt)*randn(1,N);
W=[0,cumsum(dw)];
Xt(1)=X(1);

for i=1:N-1
    Xt(i+1)=Xt(i)+(a-b.*Xt(i))*dt+sigma.*sqrt(Xt(i)).*(W(i+1)-W(i));
    dR(i+1)=max(max(dR),Xt(i+1)-n);
    dL(i+1)=max(max(dL),m-Xt(i+1));
    if Xt(i+1)>n
        Xt(i+1)=Xt(i+1)-dR(i+1);
        if Xt(i+1)<m
            Xt(i+1)=m;
        else
            Xt(i+1)=Xt(i+1);
        end
    elseif Xt(i+1)<m
        Xt(i+1)=Xt(i+1)+dL(i+1);
        if Xt(i+1)>n
            Xt(i+1)=n;
        else
            Xt(i+1)=Xt(i+1);
        end
    else
        Xt(i+1)=Xt(i+1);
    end
end
time=0:dt:T-dt;
plot(time,X,'o','LineWidth', 1, 'DisplayName', 'observations');
xlabel('Time');
ylabel('Interest Rate');
xlim([0 T]);
ylim([0.02 0.30]);
yline(0.25, 'k','y=0.25' ,'LineWidth', 1.25,'HandleVisibility','off');
yline(0.05, 'k', 'y=0.05','LineWidth', 1.25, 'HandleVisibility','off');
hold on;
plot(time,Xt,'r','LineWidth', 1,'DisplayName', 'simulation curve');
hold off;
legend('show');

