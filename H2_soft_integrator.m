%% H2_SOFT_INTEGRATOR
%
% Sintesi H2 con disturbance shaping a bassa frequenza.
%
% Il controllore usa:
%   - lo stesso plant nominale aumentato dell'LQG;
%   - gli stessi Q e R;
%   - le stesse covarianze W e V;
%   - due dinamiche lente associate alle coppie d_alpha e d_beta.
%
% Il soft integrator NON integra l'errore di tracking.
% Introduce nel problema di sintesi un modello di disturbo quasi-integrale.

close all;
clc;

%% ========================================================================
% 0. INIZIALIZZAZIONE
% ========================================================================

projectRoot = fileparts(mfilename('fullpath'));

if isempty(projectRoot)
    projectRoot = pwd;
end

cd(projectRoot);
addpath(projectRoot);

% Stessa sequenza usata negli altri script del progetto.
build_uncertain_linear_model;

projectRoot = pwd;
addpath(projectRoot);

init;

LQG_2DOF_Synthesis;

%% ========================================================================
% 1. VERIFICA VARIABILI
% ========================================================================

requiredVariables = {
    'A'
    'B'
    'Cmeas'
    'Ctrack'
    'Q'
    'R'
    'W'
    'V'
    'Gnoise'
    'P_nominal_ext'
};

for k = 1:numel(requiredVariables)

    assert( ...
        exist(requiredVariables{k},'var') == 1, ...
        'Variabile mancante: %s', ...
        requiredVariables{k});

end

A      = double(A);
B      = double(B);
Cmeas  = double(Cmeas);
Ctrack = double(Ctrack);

Q = double(Q);
R = double(R);
W = double(W);
V = double(V);

Gnoise = double(Gnoise);

n  = size(A,1);
nu = size(B,2);
ny = size(Cmeas,1);
nr = size(Ctrack,1);
nw = size(Gnoise,2);

nd = 2;

%% ========================================================================
% 2. MATRICE DELLE COPPIE AERODINAMICHE SUL MODELLO AUMENTATO
% ========================================================================
%
% P_nominal_ext ha gli ingressi:
%
%   [delta_F1;
%    delta_F2;
%    d_alpha;
%    d_beta]
%
% La parte meccanica ha quattro stati.
% Nel modello LQG gli stati degli attuatori precedono quelli meccanici.

Bd_helicopter = ...
    double(P_nominal_ext.B(:,3:4));

nh = size(Bd_helicopter,1);
na = n - nh;

assert(na >= 0, ...
    'Dimensioni incompatibili tra modello LQG e P_nominal_ext.');

Bd_augmented = [
    zeros(na,nd);
    Bd_helicopter
];

%% ========================================================================
% 3. PARAMETRI DEL SOFT INTEGRATOR
% ========================================================================
%
% Wd_i(s) = wd_i / (s + epsilon)
%
% epsilon molto piccolo -> comportamento quasi-integrale nella regione
% di frequenza di interesse.
%
% Come punto di partenza usiamo come numeratori le ampiezze delle coppie
% aerodinamiche gia' adottate nel progetto:
%
%   d_alpha = 5e-3 Nm
%   d_beta  = 2e-3 Nm
%
% soft.gain permette il tuning senza cambiare il resto della sintesi.

soft.epsilon = 1e-3;       % [rad/s]
soft.gain    = 1.0;

soft.wd = ...
    soft.gain * [
        5e-3;
        2e-3
    ];

WdSoft = diag(soft.wd);

fprintf('\n============================================================\n');
fprintf('H2 SOFT INTEGRATOR\n');
fprintf('============================================================\n');

fprintf('epsilon = %.3e rad/s\n',soft.epsilon);
fprintf('wd_alpha = %.3e\n',soft.wd(1));
fprintf('wd_beta  = %.3e\n',soft.wd(2));

%% ========================================================================
% 4. FATTORI DI Q, R, W, V
% ========================================================================

Q = 0.5*(Q + Q.');
R = 0.5*(R + R.');
W = 0.5*(W + W.');
V = 0.5*(V + V.');

Qhalf = chol(Q);
Rhalf = chol(R);

Whalf = chol(W,'lower');
Vhalf = chol(V,'lower');

%% ========================================================================
% 5. MODELLO AUMENTATO CON GLI STATI DEL DISTURBO
% ========================================================================
%
% Stato complessivo:
%
%   x_soft = [x;
%             x_d]
%
% con
%
%   xdot   = A*x + B*u + Bd_augmented*x_d + Gnoise*w
%
%   xdot_d = -epsilon*x_d + WdSoft*eta_d

Asoft = [
    A, ...
    Bd_augmented;

    zeros(nd,n), ...
    -soft.epsilon*eye(nd)
];

%% ========================================================================
% 6. INGRESSI ESOGENI
% ========================================================================
%
% eta =
%
%   [eta_w;
%    eta_v;
%    eta_d]
%
% dove:
%   eta_w = rumore di processo normalizzato
%   eta_v = rumore di misura normalizzato
%   eta_d = eccitazione normalizzata del soft integrator

B1soft = [
    Gnoise*Whalf, ...
    zeros(n,ny), ...
    zeros(n,nd);

    zeros(nd,nw), ...
    zeros(nd,ny), ...
    WdSoft
];

B2soft = [
    B;
    zeros(nd,nu)
];

%% ========================================================================
% 7. USCITE DI PRESTAZIONE
% ========================================================================
%
% Si continua a penalizzare esclusivamente lo stato fisico x e il comando u.
% Gli stati x_d appartengono al modello del disturbo e non sono direttamente
% inclusi nella funzione costo.

C1soft = [
    Qhalf, zeros(n,nd);

    zeros(nu,n), zeros(nu,nd)
];

D11soft = ...
    zeros(n + nu, nw + ny + nd);

D12soft = [
    zeros(n,nu);
    Rhalf
];

%% ========================================================================
% 8. MISURE
% ========================================================================

C2soft = [
    Cmeas, ...
    zeros(ny,nd)
];

D21soft = [
    zeros(ny,nw), ...
    Vhalf, ...
    zeros(ny,nd)
];

D22soft = ...
    zeros(ny,nu);

%% ========================================================================
% 9. GENERALIZED PLANT
% ========================================================================

P_H2_soft = ss( ...
    Asoft, ...
    [B1soft B2soft], ...
    [C1soft; C2soft], ...
    [D11soft D12soft; ...
     D21soft D22soft]);

P_H2_soft = ...
    minreal(P_H2_soft,1e-10);

%% ========================================================================
% 10. SINTESI H2
% ========================================================================

[K_H2_soft_y, ...
 CL_H2_soft, ...
 gamma_H2_soft] = ...
    h2syn( ...
        P_H2_soft, ...
        ny, ...
        nu);

K_H2_soft_y = ...
    minreal(ss(K_H2_soft_y),1e-9);

CL_H2_soft = ...
    minreal(ss(CL_H2_soft),1e-9);

fprintf('\nSINTESI\n');

fprintf( ...
    'Ordine K_H2_soft_y = %d\n', ...
    order(K_H2_soft_y));

fprintf( ...
    'gamma H2 soft = %.10g\n', ...
    gamma_H2_soft);

fprintf( ...
    'Ciclo chiuso generalized plant stabile = %d\n', ...
    isstable(CL_H2_soft));

%% ========================================================================
% 11. PREFILTRO STATICO PER IL TRACKING
% ========================================================================
%
% Il controllore restituito da h2syn chiude il problema di regolazione:
%
%       u_fb = K_H2_soft_y * y
%
% Per mantenere la stessa interfaccia di LQG introduciamo:
%
%       u = u_fb + Nbar_H2_soft*r
%
% Nbar viene calcolato sul guadagno statico DEL NUOVO closed loop,
% non viene riciclato quello dell'LQG.

Gq = ss( ...
    A, ...
    B, ...
    Ctrack, ...
    zeros(nr,nu));

Gym = ss( ...
    A, ...
    B, ...
    Cmeas, ...
    zeros(ny,nu));

LinputSoft = ...
    minreal( ...
        K_H2_soft_y * Gym, ...
        1e-9);

SuSoft = ...
    minreal( ...
        feedback( ...
            ss(eye(nu)), ...
            LinputSoft, ...
            +1), ...
        1e-9);

GvirtualToTrack = ...
    minreal( ...
        Gq * SuSoft, ...
        1e-9);

GdcSoft = ...
    dcgain(GvirtualToTrack);

if rcond(GdcSoft) < 1e-10

    warning( ...
        'Guadagno statico quasi singolare: uso pinv.');

    Nbar_H2_soft = ...
        pinv(GdcSoft);

else

    Nbar_H2_soft = ...
        inv(GdcSoft);

end

fprintf('\nNbar H2 soft:\n');
disp(Nbar_H2_soft);

%% ========================================================================
% 12. CONTROLLERE FINALE
% ========================================================================
%
% Ingressi:
%
%   [delta_r_alpha;
%    delta_r_beta;
%    delta_y_acc;
%    delta_y_mx;
%    delta_y_my]
%
% Uscite:
%
%   [delta_F1_cmd;
%    delta_F2_cmd]

Kref_H2_soft = ...
    ss([],[],[],Nbar_H2_soft);

K_H2_soft = ...
    minreal( ...
        [Kref_H2_soft K_H2_soft_y], ...
        1e-9);

K_H2_soft.InputName = {
    'delta_r_alpha'
    'delta_r_beta'
    'delta_y_acc'
    'delta_y_mx'
    'delta_y_my'
};

K_H2_soft.OutputName = {
    'delta_F1_cmd'
    'delta_F2_cmd'
};

fprintf('\nController H2 soft finale:\n');

fprintf( ...
    'Numero ingressi = %d\n', ...
    size(K_H2_soft,2));

fprintf( ...
    'Numero uscite   = %d\n', ...
    size(K_H2_soft,1));

fprintf( ...
    'Ordine          = %d\n', ...
    order(K_H2_soft));

%% ========================================================================
% 13. MATRICI PER SIMULINK
% ========================================================================

[A_H2_soft, ...
 B_H2_soft, ...
 C_H2_soft, ...
 D_H2_soft] = ...
    ssdata(K_H2_soft);

%% ========================================================================
% 14. FILTRI DI DISTURBO
% ========================================================================

s = tf('s');

Wdist_alpha = ...
    soft.wd(1) / ...
    (s + soft.epsilon);

Wdist_beta = ...
    soft.wd(2) / ...
    (s + soft.epsilon);

Wdist_soft = ...
    blkdiag( ...
        Wdist_alpha, ...
        Wdist_beta);

figure('Name','Soft integrator disturbance shaping');

sigma(Wdist_soft);

grid on;

title( ...
    'Disturbance shaping del controllore H_2 soft');

%% ========================================================================
% 15. SALVATAGGIO
% ========================================================================

save( ...
    'H2_soft_controller.mat', ...
    'K_H2_soft', ...
    'K_H2_soft_y', ...
    'P_H2_soft', ...
    'CL_H2_soft', ...
    'gamma_H2_soft', ...
    'Nbar_H2_soft', ...
    'soft', ...
    'Wdist_soft', ...
    'Bd_augmented', ...
    'A_H2_soft', ...
    'B_H2_soft', ...
    'C_H2_soft', ...
    'D_H2_soft');

fprintf('\nH2_soft_controller.mat salvato correttamente.\n');