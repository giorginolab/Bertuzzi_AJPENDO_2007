% Bertuzzi's beta cell model, as described in  Alessandro Bertuzzi,
% Serenella Salinari and Geltrude Mingrone "Insulin granule trafficking in
% β-cells: mathematical model of glucose-induced insulin secretion" Am J
% Physiol Endocrinol Metab 293:E396-E409, 2007.
% doi:10.1152/ajpendo.00647.2006.


% Warning: gamma overrides a built-in function 

% Naming convention: not using underscores after single-letter variables


% Additional Variables (unused)
% ---------------------
% C Pool of unbound L-type Ca2+ channels
% G Extracellular glucose concentration (input variable) (mmol/l) (Eq. 7)
% ISR Insulin secretion rate (pmol/min) (Eq. 11)
% IS Insulin amount secreted in a time T during the first-phase response (pmol) (Eq. A4)
% Gm Measured plasma glucose concentration (mmol/l) (Eq. C1)
% Ip Measured plasma insulin concentration (pmol/l) (Eq. C1)
% Ga Estimated arterial glucose concentration (mmol/l) (Eq. C1)
% 


function dx=BertuzziModel(t,x,Z)

%% State Variables (Equation Number Where Quantities First Appear)
% ---------------------
% I Pool of proinsulin aggregates (Eq. 1)
% V Pool of granule membrane material (Eq. 1)
% R Reserve pool (Eq. 3)
% D Pool of docked granules (Eq. 4)
% DIR Pool of immediately releasable granules (Eq. 4)
% F Pool of granules fused with plasma membrane (Eq. 6)
% gamma Rate coefficient of granule externalization and priming, related to the ATP-to-ADP ratio (min^-1) (Eq. 3)
% rho  rho Rate coefficient of granule fusion with cell membrane, related to [Ca2+]i (min^-1) (Eq. 5)
    tmp=num2cell(x);
    [I V R D DIR F gamma rho]=tmp{:};
 
    
%% Parameters in the Equations of Granule Dynamics
% ---------------------
% k Rate constant of formation of proinsulin-containing granules in trans-Golgi network (min^-1) (Eq. 1)
k=1e-2;
% bI Biosynthesis rate of proinsulin aggregates at a given glucose concentration (min^-1) (Eq. 1)
bI=4;
% alpha_I Rate constant of degradation of proinsulin aggregates (min^-1) (Eq. 1)
alpha_I=0.3;
% bV Rate of biosynthesis of granule membrane material (min^-1) (Eq. 2)
bV=6;
% alpha_V Rate constant of degradation of granule membrane material (min^-1) (Eq. 2)
alpha_V=2*alpha_I;
% tau_V Time delay related to recycling of granule membrane material (min) (Eq. 2)
tau_V=5;
% CT Pool of total Ca2+ channels (Eq. 4)
CT=500;
% k1p + Rate constant of association for the binding between granule and Ca2+ channel (min^-1) (Eq. 4)
k1p=1.447e-5;
% k1m - Rate constant of dissociation for the binding between granule and Ca2+ channel (min^-1) (Eq. 4)
k1m=0.10375;
% sigma Rate constant of insulin release from granules fused with cell membrane (min^-1) (Eq. 2)
sigma=30;


% Additional Parameters and Functions
% ---------------------
% I0 Insulin amount contained in a granule (amol) (Eq. 11)
I0=1.6;
% Nc Average number of β-cells in an islet (Eq. 11)
Nc=1000;
% Ni Number of islets in pancreas (Eq. 11)
Ni=1;
% N Total number of β-cells in pancreas (Eq. 17)
N=1;
% f Fraction of β-cells responding to glucose stimulus (function of glucose concentration G) (Eq. 11)

% fb Basal value of the fraction f (Eq. A5)
fb=0.05;
% Kf Parameter in the function f(G) (mmol/l) (Eq. A5)
Kf=3.43;

% VG Glucose distribution volume (liters) (Eq. C1)
% VC Volume of C-peptide accessible compartment (liters) (Eq. 17)
% Q Cardiac output (l/min) (Eq. C1)
% SI Insulin sensitivity index [l min^-1 (pmol/l)^-1] (Eq. C1)
% Sd Dynamic responsivity index (Eq. 17)
% Ss Static responsivity index (min^-1) (Eq. 18)
% a, b Parameters in the approximate expression of the first-phase ISR (min^-1) (Eq. A1)
% xb Denotes the basal value of the generic variable x
% x_bar Denotes the steady-state value of the variable x

% 
% Parameters and Functions in the Equations Representing Stimulus-secretion Coupling
% ---------------------
% eta Rate constant in the equation for gamma (min^-1) (Eq. 7)
eta=4;
% gamma_b Basal value of gamma (min^-1) (Eq. 7)
gamma_b=1e-4;
% psi Oscillatory function that represents events inducing oscillations in gamma (min^-1) (Eq. 7)
% BELOW
% h_gamma Function of G representing the activatory action of glucose on gamma (min^-1) (Eq. 7)
% BELOW
% tau_G Time delay related to time required by glucose metabolism for activation of gamma (min) (Eq. 7)
tau_G=1;
% G_star Glucose concentration threshold for the activation of gamma (mmol/l) (Eq. 8)
G_star=4.58;
% h_hat Maximal value of hgamma (min^-1) (Eq. 8)
h_hat=3.93e-3;
% G_hat Glucose concentration over which hgamma remains constant and equal to h  (mmol/l) (Eq. 8)
G_hat=10;
% zeta Rate constant in the equation for rho (min^-1) (Eq. 9)
zeta=4;
% rho_b Basal value of rho (min^-1) (Eq. 9)
rho_b=0.02;
% h_rho Function of gamma representing the activatory action of rho (min^-1) (Eq. 9)
% BELOW
% k_rho Parameter representing the sensitivity of rho on the activatory action of gamma (Eq. 10)
k_rho=350;



    %% ---------------------
    dI=-k*I*V-alpha_I*I+bI;                     % (1)
    F_delayed=Z(6,1);  % F(t-tau_V) as per DDE
    dV=-k*I*V-alpha_V*V+bV+sigma*F_delayed;     % (2)
    dR=k*I*V-gamma*R;                           % (3)
    dD=gamma*R-k1p*(CT-DIR)*D+k1m*DIR;          % (4)
    dDIR=k1p*(CT-DIR)*D-k1m*DIR-rho*DIR;        % (5)
    dF=rho*DIR-sigma*F;                         % (6)
    dgamma=eta*(-gamma+gamma_b+psi(t)+h_gamma(G(t-tau_G))); % (7)
    drho=zeta*(-rho+rho_b+h_rho(gamma));        % (9)
    
    dx= [dI dV dR dD dDIR dF dgamma drho]';



%% f Fraction of beta-cells responding to glucose stimulus (function of glucose concentration G) (Eq. 11)
function out=f(G)
    if G<G_star
        out=fb;
    else
        out=fb+(1-fb)*(G-G_star)/(Kf+G-G_star);
    end
end

%% h_gamma Function of G representing the activatory action of glucose on gamma (min^-1) (Eq. 8)
function out=h_gamma(G)
    if G<=G_star
        out=0;
    elseif G_star<G && G<=G_hat
        out=h_hat*(G-G_star)/(G_hat-G_star);
    else
        out=h_hat;
    end
end

%% h_rho Function of gamma representing the activatory action of rho (min^-1) (Eq. 10)
function out=h_rho(gamma)
    if gamma<gamma_b
        out=0;
    else
        out=k_rho*(gamma-gamma_b);
    end
end

%% psi Oscillatory function that represents events inducing oscillations in gamma (min^-1) (Eq. 7)
function out=psi(t)
    out=0;
end

% Glucose (input) - depends on the experiment
function out=G(t)
    out=30;
end



end

