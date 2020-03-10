clear; clc;

mh = 125;
Mp = 2.435*10^18;
GF = 1.166378*10^-5;
v = 1/sqrt(sqrt(2)*GF);
mt = 160;
yt = sqrt(2)*(mt/v);
mW = 80.385;
hbar = 6.58*10^(-25);
kB = 8.6173303*10^-5*10^-9;
T0 = kB*2.7255;
h = 0.673;
H0 = (0.678*hbar)/(9.777752*10^9*365.25*24*3600);
s0 = (4/3)*(pi^2/30)*(2+21/11)*T0^3;
alphaem = 1/128;
sw2 = 0.23129;
g = sqrt((8/sqrt(2))*GF*mW^2);
gprime = sqrt((4 *pi *alphaem)/(1 - sw2));

gchi = 2;
gstar = 100;

n=0;
tol=0.1;

for i=1:10^9

lambda1 = mh^2/v^2;
lambda2 = rand;
lambda3 = 3*rand;
lambda4 = rand;
lambda5 = -rand;
l345 = lambda3 + lambda4 + lambda5;

c1 = (3*lambda1 + 2*lambda3 + lambda4)/12 + (3*g^2 + gprime^2)/32 + yt^2/8;
c2 = (3*lambda2 + 2*lambda3 + lambda4)/12 + (3*g^2 + gprime^2)/32;

mchi = 100*(1+10*rand);
Mspm = mchi*(1+rand);

muH2 = lambda1*v^2;
mueta2 = lambda3*v^2 - 2*Mspm^2;

xa = (sqrt(c2)*mchi)/sqrt(mueta2);
xb = sqrt((sqrt(lambda2)*c1*mchi^2 - sqrt(lambda1)*c2*mchi^2)/(sqrt(lambda2)*muH2 - sqrt(lambda1)*mueta2));

if xa > 1 && imag(xa) == 0 && imag(xb) == 0 && xa < xb && Mspm > mchi
    
    Omegachih2 = 0.112;
    Yinf = (3*H0^2*Mp^2*Omegachih2)/(mchi*s0*h^2);

    fun = @(x) (45/(2*pi^4))*(pi/8)^(1/2)*(gchi/gstar)*x.^(3/2).*exp(-x)/(1+0.006*(100/gstar)*(c1*mueta2 - c2*muH2)/(sqrt(lambda1*lambda2)*(mchi/x)^2))-Yinf;
    
    xbactual = fzero(fun,30);
    
    sfsi = 1+0.006*(100/gstar)*(c1*mueta2 - c2*muH2)/(sqrt(lambda1*lambda2)*(mchi/xbactual)^2);
    
    if xbactual > 1 && xb > xbactual-tol && xb < xbactual+tol
    n=n+1;
    
    x(n,:)=[xa xb xbactual lambda1 lambda2 lambda3 lambda4 lambda5 mchi Mspm sfsi]
    end
end

end

sort(x)
