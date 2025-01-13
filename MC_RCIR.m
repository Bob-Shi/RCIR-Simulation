clear;
T=250;N=1000;dt=T/N;r=1000;
a=0.25;b=0.5;sigma=0.12;m=0.4;n=0.6;
rng(100);    
dw=sqrt(dt)*randn(r,N-1);
W=[zeros(r,1),cumsum(dw,2)];
Xt=zeros(r,N);dR=zeros(r,N);dL=zeros(r,N);
Xt(:,1)=0.5;dR(:,1)=0;dL(:,1)=0;
errtop_a=zeros(r,1);
errtop_b=zeros(r,1);
errbottom=zeros(r,1);
err_a=zeros(r,1);
err_b=zeros(r,1);
At=zeros(r,N);
Bt=zeros(r,N);
Mt=zeros(r,N);
Nt=zeros(r,N);

for z=1:r

    for i=1:N-1
        time(i+1)=i*dt;
        Xt(z,i+1)=Xt(z,i)+(a-b.*Xt(z,i))*dt+sigma.*sqrt(Xt(z,i)).*(W(z,i+1)-W(z,i));
        dR(z,i+1)=max(max(dR(z,:)),Xt(z,i+1)-n);
        dL(z,i+1)=max(max(dL(z,:)),m-Xt(z,i+1));
        if Xt(z,i+1)>n
            Xt(z,i+1)=Xt(z,i+1)-dR(z,i+1);
            if Xt(z,i+1)<m
                Xt(z,i+1)=m;
            else
                Xt(z,i+1)=Xt(z,i+1);
            end
        elseif Xt(z,i+1)<m
            Xt(z,i+1)=Xt(z,i+1)+dL(z,i+1);
            if Xt(z,i+1)>n
                Xt(z,i+1)=n;
            else
                Xt(z,i+1)=Xt(z,i+1);
            end
        else
            Xt(z,i+1)=Xt(z,i+1);
        end
        At(z,i+1)=At(z,i)+Xt(z,i)*dt;
        Bt(z,i+1)=Bt(z,i)+1/Xt(z,i)*dt;
        Mt(z,i+1)=Mt(z,i)+1/sqrt(Xt(z,i))*dw(z,i);
        Nt(z,i+1)=Nt(z,i)+sqrt(Xt(z,i))*dw(z,i);
    end
    errtop_a(z,1)=At(z,N)*Mt(z,N)-T*Nt(z,N);
    errtop_b(z,1)=T*Mt(z,N)-Bt(z,N)*Nt(z,N);
    errbottom(z,1)=At(z,N)*Bt(z,N)-T*T;
    err_a(z,1)=sigma*errtop_a(z,1)/errbottom(z,1);
    err_b(z,1)=sigma*errtop_b(z,1)/errbottom(z,1);
end
err_ahat=sum(err_a)/r;
err_bhat=sum(err_b)/r;
var_a=sum(err_a.^2)/r;
var_b=sum(err_b.^2)/r;
plot(time,Xt(1,:));