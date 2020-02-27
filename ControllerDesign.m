
%%
load newtankdata.mat  
t  = TankData.Time;
pV = TankData.Data(:,1); % input voltage pump
rh1 = TankData.Data(:,2); % Height of Tank 1
rh2 = TankData.Data(:,3); % Height of Tank 2
% vT = TankData.Data(:,4); % valve 
% v0 = TankData.Data(:,5); % valve

%%
%Parameters
Alpha=0.0515;
Beta=0.0673;
Gamma=0.0911;
Phi=0.0069;

%Constants around Operating Point
Vt = 1;
V01 = 0.5;
V02 = 0;
h1 = 0.602;
h2 = 0.313;
u = 1.2;
g = 9.8;


%% State Space Setup
% A Matrix Terms
A00 = (-Phi*(Vt*g)/sqrt(2*g*(h1-h2)))-((Beta*V01*g)/sqrt(2*g*h1));

A01 = (Phi*(Vt*g)/sqrt(2*g*(h1-h2)));

A10 = (Phi*(Vt*g)/sqrt(2*g*(h1-h2)));

A01 = 1.2*A01;
A10 = A01; 

A11 = (-Phi*(Vt*g)/sqrt(2*g*(h1-h2)))-((Gamma*V02*g)/sqrt(2*g*h2));

A11 =  2.2 * A11;

% B Matrix Terms
B00 = 1.15*Alpha;
B01 = Phi*(sqrt(2*g*h1));
B02 = 0;
B03 = -Gamma*(sqrt(2*g*(h1-h2)));
B10 = 0;
B11 = 0;
B12 = -Beta*(sqrt(2*g*h2));
B13 = Gamma*(sqrt(2*g*(h1-h2)));

% State Matrices
A = [A00 A01;A10 A11]

B = [B00 B01 B02 B03;B10 B11 B12 B13]

C = eye(2)

D = [0 0 0 0; 0 0 0 0]
sys = ss(A,B,C,D);
systf = tf(sys)

%% Graphing Linear & Non Linear 
% Subplot
% subplot(2,1,1); 
% stepplot(systf(1,1),150); hold on;
% plot(t,rh1,'r');
% legend('h_1 linear','h_1 nonlinear')
% axis([0 150 0 1]) 
% 
% subplot(2,1,2); 
% stepplot(systf(2,1),150); hold on;
% plot(t,rh2,'r');
% legend('h_2 linear','h_2 nonlinear')
% axis([0 150 0 1]) 


%% State-Feedback Controller
% LQR Integrator Design
C = [0 1];  %% height of the second tank


%% Check controllability 
Ap = [ zeros(1,1) C;
       zeros(2,1) A];
   
Bp = [zeros(1,1); B(:,1)];

Ap = [ zeros(1,1) C;
       zeros(2,1) A];
   
Bp = [zeros(1,1); B(:,1)];


rank(ctrb(Ap,Bp))
rank(Ap)

%% Q & R Matrices 
Q = diag([0.0001 10000 1000000000]); % Penalyty Matrice Setup
Q = diag([100 1 100]);
R =10*eye(1);                        % Cost Matrice

K = lqr(Ap,Bp,Q,R)                   % Gain vector

%% New Closed Loop System
Ar = Ap - Bp * K
Br = [-1; 0; 0]
Cr = [0 0 1]

syscl = ss(Ar, Br, Cr, []);
step(syscl);