close all

load newtankdata.mat 

%%
t  = newtankdata.Time;
V   = newtankdata.Data(:,1); % input voltage pump
rh1 = newtankdata.Data(:,2); % Height of Tank 1
rh2 = newtankdata.Data(:,3); % Height of Tank 2

%%
%Parameters
Alpha=0.0515;
Beta=0.0673;
Gamma=0.0911;
Phi=0.0069;

%Constants around Operating Point
Vt = 1;
V01 = 1;
V02 = 1;
h1 = 0.602;
h2 = 0.27;
u = 1.2;
g = 9.8;

%% State Space Setup
% A Matrix Terms
A00 = (-Phi*(Vt*g)/sqrt(2*g*(h1-h2)))-((Beta*V01*g)/sqrt(2*g*h1));

A01 = (Phi*(Vt*g)/sqrt(2*g*(h1-h2)));

A10 = (Phi*(Vt*g)/sqrt(2*g*(h1-h2)));

A11 = (-Phi*(Vt*g)/sqrt(2*g*(h1-h2)))-((Gamma*V02*g)/sqrt(2*g*h2));

%% scale parameters
A00 = 6 * A00;
A01 = 10 * A01;
A11 =  A11;

% B Matrix Terms
B00 = 15.5*Alpha;
B01 = Phi*(sqrt(2*g*h1));
B02 = 0;
B03 = -Gamma*(sqrt(2*g*(h1-h2)));
B10 = 0;
B11 = 0.11;
B12 = -Beta*(sqrt(2*g*h2));
B13 = Gamma*(sqrt(2*g*(h1-h2)));

B12 = 6 * B12; 

B12 = 10 * B12; 

% State Matrices
A = [A00 A01; A01 A11]

B = [B00 B01 B02 B03;B10 B11 B12 B13]

C = eye(2)

D = [0 0 0 0; 0 0 0 0]
sys = ss(A,B,C,D);
systf = tf(sys)

opt = stepDataOptions('StepAmplitude',1.2);

subplot(2,1,1)
plot(t,rh1,'r'); hold on; stepplot(systf(1,1),t(end));
subplot(2,1,2) 
plot(t,rh2,'r'); hold on; stepplot(systf(2,2),t(end));


%% lq integrator design
C = [0 1];  %% height of the second tank


%% check controllability 
Ap = [ zeros(1,1) C;
       zeros(2,1) A];
   
Bp = [zeros(1,1); B(:,1)];

%% penalty matrices 
Q = 100 * eye(3) %% on states
% Q = diag([1 0 1]);
R = 1*eye(1);     
K = lqr(Ap,Bp,Q,R)

%% new closed loop
Ar = Ap - Bp * K;
Br = [-1; 0; 0];
Cr = [0 0 1];

syscl = ss(Ar, Br, Cr, []);
step(syscl);