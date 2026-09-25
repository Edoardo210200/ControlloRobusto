%% H2_FROM_LQG
% Formulazione H2 equivalente al problema LQG del progetto.

close all;
clc;

%% 1. Eseguo la sintesi LQG gia' validata
run('LQG_2DOF_Synthesis.m');

requiredVars = { ...
    'A','B','Cmeas','Q','R','W','V','Gnoise', ...
    'Nbar','K_LQG'};

for k = 1:numel(requiredVars)
    assert(exist(requiredVars{k},'var') == 1, ...
        'Variabile mancante: %s',requiredVars{k});
end

A      = double(A);
B      = double(B);
Cmeas  = double(Cmeas);
Q      = double(Q);
R      = double(R);
W      = double(W);
V      = double(V);
Gnoise = double(Gnoise);

n  = size(A,1);
nu = size(B,2);
ny = size(Cmeas,1);
nw = size(Gnoise,2);
nr = size(Nbar,2);

%% 2. Fattori dei pesi e delle covarianze

Q = 0.5*(Q+Q');
R = 0.5*(R+R');
W = 0.5*(W+W');
V = 0.5*(V+V');

% Qhalf'*Qhalf = Q
Qhalf = chol(Q);

% Rhalf'*Rhalf = R
Rhalf = chol(R);

% Whalf*Whalf' = W
Whalf = chol(W,'lower');

% Vhalf*Vhalf' = V
Vhalf = chol(V,'lower');

%% 3. Generalized plant H2

% Ingressi esogeni:
% eta = [eta_w ; eta_v]
%
% xdot = A*x + Gnoise*Whalf*eta_w + B*u
%
% z = [Qhalf*x ; Rhalf*u]
%
% y = Cmeas*x + Vhalf*eta_v

B1 = [ ...
    Gnoise*Whalf, ...
    zeros(n,ny) ...
    ];

B2 = B;

C1 = [ ...
    Qhalf;
    zeros(nu,n) ...
    ];

D11 = zeros(n+nu,nw+ny);

D12 = [ ...
    zeros(n,nu);
    Rhalf ...
    ];

C2 = Cmeas;

D21 = [ ...
    zeros(ny,nw), ...
    Vhalf ...
    ];

D22 = zeros(ny,nu);

P_H2 = ss( ...
    A, ...
    [B1 B2], ...
    [C1; C2], ...
    [D11 D12; D21 D22]);

%% 4. Sintesi H2

[K_H2_y,CL_H2,gamma_H2] = ...
    h2syn(P_H2,ny,nu);

K_H2_y = minreal(ss(K_H2_y),1e-9);

fprintf('\n=============================================\n');
fprintf('H2 EQUIVALENTE AL PROBLEMA LQG\n');
fprintf('=============================================\n');

fprintf('Ordine K_H2_y = %d\n',order(K_H2_y));
fprintf('||Tzw||_2 H2  = %.10g\n',gamma_H2);
fprintf('Stabile        = %d\n',isstable(CL_H2));

%% 5. Valuto LQG sullo stesso identico generalized plant

% K_LQG ha ingressi [r ; y_m].
% Nel problema di regolazione estraggo solamente y_m -> u.

K_LQG_y = ...
    K_LQG(:,nr+(1:ny));

CL_LQG_as_H2 = ...
    minreal( ...
        lft(P_H2,K_LQG_y), ...
        1e-9);

gamma_LQG_as_H2 = ...
    norm(CL_LQG_as_H2,2);

relGap = ...
    abs(gamma_H2-gamma_LQG_as_H2) / ...
    max(gamma_H2,eps);

fprintf('\nVERIFICA EQUIVALENZA\n');

fprintf('||Tzw||_2 H2  = %.10g\n', ...
    gamma_H2);

fprintf('||Tzw||_2 LQG = %.10g\n', ...
    gamma_LQG_as_H2);

fprintf('Scarto relativo = %.3e\n', ...
    relGap);

if relGap < 1e-4
    fprintf('Equivalenza LQG-H2 verificata.\n');
else
    warning(['Scarto non trascurabile: controllare ', ...
             'generalized plant e convenzioni di segno.']);
end

%% 6. Stessa interfaccia del controllore LQG per il tracking
%
% u = Nbar*r + K_H2_y*y_m

Kref = ss([],[],[],Nbar);

K_H2 = ...
    minreal( ...
        [Kref K_H2_y], ...
        1e-9);

fprintf('\nController H2 finale:\n');
fprintf('Ingressi = %d\n',size(K_H2,2));
fprintf('Uscite   = %d\n',size(K_H2,1));
fprintf('Ordine   = %d\n',order(K_H2));

%% 7. Confronto frequenziale sulla stessa mappa H2

omega = logspace(-3,3,700);

figure;
sigma( ...
    CL_H2, ...
    CL_LQG_as_H2, ...
    omega);

grid on;

legend( ...
    'H_2 via h2syn', ...
    'LQG', ...
    'Location','best');

title('Equivalenza LQG - H_2');

%% 8. Salvataggio

save('H2_LQG_controller.mat', ...
    'K_H2', ...
    'K_H2_y', ...
    'K_LQG_y', ...
    'P_H2', ...
    'CL_H2', ...
    'CL_LQG_as_H2', ...
    'gamma_H2', ...
    'gamma_LQG_as_H2', ...
    'relGap');