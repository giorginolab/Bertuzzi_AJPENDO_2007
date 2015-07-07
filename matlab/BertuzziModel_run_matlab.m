
% Still needs work, esp. for the various experiments

tau_V=5;
x0=BertuzziModelInit;   % initial conditions
xH=BertuzziModelInit;   % history

% Not explicit, but it is a DDE
sol=dde23('BertuzziModel',[tau_V],xH,[0 60]);

%     [I V R D DIR F gamma rho]=tmp{:};

I0=1.6;
sigma=30;

ISR=sol.y(6,:)*I0*sigma; % omitting other factors

%% Figure 5? Forgot it.
plot(sol.x,ISR);
ylabel 'ISR (units?)'
xlabel 'Time (min)'



