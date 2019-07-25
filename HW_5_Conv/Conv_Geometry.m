%% Geometry
%%addpath( [pwd '\my_femm']);
clear all;

%%
%%%%%%%%%%%%%%%%%%%%%%% Machine dimensions
P_d = 10000;                       % Developed Power in Watts
I_max =  20;                       % Maximum amperes allowed from the supply or the source
W_m = 523.59;                      % Base speed in rad/s 
N_rev = (60/pi)*(W_m)                % Speed in RPM (Speed is controlled by frequency of switching )
m = 1;                             % Number of phases conducting simultaniously
 A_s = 50000 ;                      % Specific electric loading A_s = (2*T_ph*i*m)/(pi*D)
                                    % 25000 < A_s < 90000
                                    
 K = 0.5;                              % Non Servo applications 0.25 < K < 0.7
                                    % Servo applications 1 < K < 3
 L = k*D ;                          % Axial length of the machine or stack length taken as multiple or sub multiple of 'D'
 N = 3  ;                           % Number of stator phases(P_s /2)




Ns = 6;                         % Stator poles  
Nr = 4;                         %Rotor Poles =  rotor teeth(360deg)
T  = 250;                       %Number of turns 
d_betar = 1/3;                  %Ratio of d_betar = Beta_r/(2pi/Nr)
Ts = 3;                         %Stator teeth per pole face
Thetea_sp =  ((Ts - 1)+d_betar)*(2*pi/Nr) % stator pole angle
Thetea_rp = (2*pi/Nr)*d_betar   %Rotor pole angle
g = 0.25e-3;                     %Air Gap
b_s = 10e-3;                      %height of stator Pole
R1  = 36.10e-3;%35.6e-3;         %Radius of rotor pole from centre
h_teeth = 0.5*b_s;               %teeth slot depth inside stator 
R2  = R1+g+b_s+h_teeth;          %Stator inner radius      
R3  = R1+g+h_teeth;              %Outer radius of stator teeth 
h_r = 2*R1*sin( (pi/Nr) - (Thetea_rp/2) )*0.4;%1.5e-2;%2.4e-2;  %Rotor pole height from rotor inner circle
y_s = sin(Thetea_sp/2)*(R1 +g+ h_teeth)*1.4;%1.25e-2;
y_r = 2.6e-2;%sin(Thetea_rp/2)*R1*1.4; %2.6e-2; // this conditon gets voilated if the 
                               % number of rotor poles gets increased 
%t_r = 2*sin(Thetea_rp/2)*R1; % Rotro pole thickness
%t_s = 2*sin(Thetea_sp/2)*(R1+g+h_teeth);% stator pole thickness 
S_r = 0.4*h_r;
L   = 5e-3;                     %Stack length
mmf_sel = (Imax*T)/2;

%% Electrical parameters 
R = 0.05;
Vdc =  100;     % Dc voltage 
T_sim = 0.3e-3;    % Simulation time
Wm = (45*pi)/(180*T_sim);

%% Calculate number of points
% Calculate points where magnetic field has to be measured in...
% in the air gap
r_middle = R1+(g/2);
n = 20; % number of points to be computed 
theta_point = linspace((-45)*(pi/180),0*(pi/180),n);
for idx =1 :1:n
x(idx,:) = r_middle*cos(theta_point(idx)); % x(Theta)column vector
y(idx,:) = r_middle*sin(theta_point(idx)); % y(Theta)column vecotr
end

%% Compute stator tooth points

if (mod(Ts,2) == 1 )
    % Odd angles
    % ((2*pi*i)/Nr)+((-1)^i)*(Thetea_rp/2); %i =0,1,2,3,4.....Ts-1
    
    for i = 0:1:(Ts-1)
    % Positive angles
    if (mod(i,2)==0)
    Theta_teeth_p(i+1) = ((pi*i)/Nr)+((-1)^i)*(Thetea_rp/2);
    
    else
         Theta_teeth_p(i+1) = ((pi*(i+1))/Nr)+((-1)^i)*(Thetea_rp/2);
    end
    
    end
    
  
    for idx = 1:1:Ts
    %Co-ordinates of inner radius of stator teeth
    x_in_teeth_p(idx,:) = (R1+g)*cos(Theta_teeth_p(idx));
    y_in_teeth_p(idx,:) = (R1+g)*sin(Theta_teeth_p(idx));
    x_in_teeth_n(idx,:) = (R1+g)*cos(-Theta_teeth_p(idx));
    y_in_teeth_n(idx,:) = (R1+g)*sin(-Theta_teeth_p(idx));
    
    %Co-ordinates of outer radius of stator teeth
    x_out_teeth_p(idx,:) = (R1+g+h_teeth)*cos(Theta_teeth_p(idx));
    y_out_teeth_p(idx,:) = (R1+g+h_teeth)*sin(Theta_teeth_p(idx));
    x_out_teeth_n(idx,:) = (R1+g+h_teeth)*cos(-Theta_teeth_p(idx));
    y_out_teeth_n(idx,:) = (R1+g+h_teeth)*sin(-Theta_teeth_p(idx));
    
    end
    
else 
    % check this formulae
    % Even angles
    % ((pi*(2*i+1))/Nr)+((-1)^(i+1))*(Thetea_rp/2); %i =0,1,2,3,4.....Ts-1
    
    for i = 0:1:(Ts-1)
    % Positive angles
    Theta_teeth_p(i+1) = ((pi*(2*i+1))/Nr)+((-1)^(i+1))*(Thetea_rp/2);
    
    end
    
    for idx = 1:1:Ts
    x_in_teeth_p(idx,:) = (R1+g)*cos(Theta_teeth_p(idx));
    y_in_teeth_p(idx,:) = (R1+g)*sin(Theta_teeth_p(idx));
    x_in_teeth_n(idx,:) = (R1+g)*cos(-Theta_teeth_p(idx));
    y_in_teeth_n(idx,:) = (R1+g)*sin(-Theta_teeth_p(idx));
    
    
    x_out_teeth_p(idx,:) = (R1+g+h_teeth)*cos(Theta_teeth_p(idx));
    y_out_teeth_p(idx,:) = (R1+g+h_teeth)*sin(Theta_teeth_p(idx));
    x_out_teeth_n(idx,:) = (R1+g+h_teeth)*cos(-Theta_teeth_p(idx));
    y_out_teeth_n(idx,:) = (R1+g+h_teeth)*sin(-Theta_teeth_p(idx));
    
    end
    
    
    
end
%% Create structure for Geometry
init_geo.Ns = Ns;
init_geo.Nr = Nr;
init_geo.N = N;
init_geo.Th_sp = Thetea_sp;
init_geo.Th_rp = Thetea_rp;
init_geo.Th = 0*pi/180;
%init_geo.Th_ss = Thetea_ss;
%init_geo.Th_rs = Thetea_rs;
init_geo.Th= 0*pi/180;
init_geo.g = g;
%init_geo.ts = t_s;
%init_geo.tr = t_r;
init_geo.ys = y_s;
init_geo.yr = y_r;
init_geo.R1 = R1;
init_geo.R2 = R2;
init_geo.hr = h_r;
init_geo.L_stk = L;

%init_geo.hu = h_u;
%init_geo.R11 = R11;
%init_geo.tm = t_m;
%init_geo.ts = t_s;
%init_geo.tb = t_b;
%init_geo.beta_ri = beta_ri;
%init_geo.beta_ro = beta_ro;
%init_geo.h_r1 = h_r1;
%init_geo.h_r2 = h_r2;
%init_geo.R_By_o = R_By_o;
%init_geo.R_By_i = R_By_i;
init_geo.x = x;
init_geo.y = y;
init_geo.theta_point = theta_point;

%%%% Teeth changes added here
init_geo.h_teeth = h_teeth;
init_geo.Theta_teeth_p = Theta_teeth_p;

init_geo.x_out_teeth_p = x_out_teeth_p;
init_geo.y_out_teeth_p = y_out_teeth_p;

init_geo.x_out_teeth_n = x_out_teeth_n;
init_geo.y_out_teeth_n = y_out_teeth_n;


init_geo.x_in_teeth_p = x_in_teeth_p;
init_geo.y_in_teeth_p = y_in_teeth_p;

init_geo.x_in_teeth_n = x_in_teeth_n;
init_geo.y_in_teeth_n = y_in_teeth_n;


save('init_geo.mat','init_geo','mmf_sel');

 