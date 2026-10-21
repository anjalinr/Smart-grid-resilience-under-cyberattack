function smart_grid_FDIA_IEEE14()
% SMART GRID RESILIENCE UNDER FALSE DATA INJECTION ATTACK - IEEE 14-BUS
% Network data : MATPOWER case14 (IEEE 14-bus test case, UW Power Systems Test Case Archive)
%                bus loads Pd and branch reactance x / tap ratio only (DC model).
% Model        : x(k+1) = x(k) + B*u(k) + w(k)      (A = I, quasi-steady-state DC angles)
%                y(k)   = H*x(k) + v(k) + a(k)      (attack a(k) added to received data)
%   x : 13 bus voltage angles (rad), bus 1 = reference
%   u : scheduled generation change at bus 2 (known input), B = Binv(:,bus 2)
%   w : random slow load variation mapped to angles, w = Binv*dP, Q = cov(w)
%   y : 34 meters = 20 line active-power flows + 14 bus power injections
% Detection    : Kalman filter, a-priori residual r = y - H*xpred,
%                J = r'*inv(S)*r vs chi-square threshold; static WLS bad-data test for comparison
% Mitigation   : iterative largest-normalised-residual isolation inside the filter
% ASSUMED (not from a source): sensor noise, load-variation level, attack sizes, false-alarm rate.

rng(1);

%% ---------- MATPOWER case14 data (DC model) ----------
base = 100;                                   % MVA
Pd = [0 21.7 94.2 47.8 7.6 11.2 0 0 29.5 9 3.5 6.1 13.5 14.9]';   % MW, buses 1..14
br = [ 1  2 0.05917 1;   1  5 0.22304 1;   2  3 0.19797 1;   2  4 0.17632 1;
       2  5 0.17388 1;   3  4 0.17103 1;   4  5 0.04211 1;   4  7 0.20912 0.978;
       4  9 0.55618 0.969; 5 6 0.25202 0.932; 6 11 0.19890 1; 6 12 0.25581 1;
       6 13 0.13027 1;   7  8 0.17615 1;   7  9 0.11001 1;   9 10 0.08450 1;
       9 14 0.27038 1;  10 11 0.19207 1;  12 13 0.19988 1;  13 14 0.34802 1 ];
nb = 14; nl = size(br,1);
Ainc = zeros(nl,nb);
for l = 1:nl, Ainc(l,br(l,1)) = 1; Ainc(l,br(l,2)) = -1; end
bl   = 1./(br(:,3).*br(:,4));                 % branch susceptance incl. tap (as in MATPOWER DC)
Bf   = diag(bl)*Ainc;   Bbus = Ainc'*Bf;
Hf = Bf(:,2:nb);  Hp = Bbus(:,2:nb);          % reduced (slack angle = 0)
H  = [Hf; Hp];  nm = size(H,1);               % 34 x 13
Binv = inv(Bbus(2:nb,2:nb));

%% ---------- Simulation settings (assumed) ----------
N = 200; k_on = 80; k_off = 140; k = 1:N;
sf = 0.01; sp = 0.01;                         % meter noise (p.u.) = 1 MW
spmu = 0.0005;                                % trusted PMU angle noise (rad)
R  = diag([sf*ones(nl,1); sp*ones(nb,1)].^2);
sigP = max(0.004*Pd(2:nb)/base, 1e-5);        % per-step load variation (p.u.)
Q  = Binv*diag(sigP.^2)*Binv';
u0 = 0.0005;                                  % gen-2 scheduled ramp per step (p.u.)
Bu = Binv(:,1);                               % bus 2 is first reduced bus
pfa = 1e-4;                                   % designed false-alarm probability per sample
chi = zeros(1,40);
for d = 1:40, chi(d) = 2*gammaincinv(pfa, d/2, 'upper'); end

P0inj = -Pd/base;  P0inj(2) = P0inj(2) + 0.40;  P0inj(1) = -sum(P0inj(2:nb));
p.x0 = Binv*P0inj(2:nb);  p.P0 = 1e-4*eye(13);
p.N = N; p.Q = Q; p.Bu = Bu; p.u0 = u0; p.chi = chi; p.Hp = Hp; p.base = base;
% Noise is stored in sim_noise.csv so MATLAB and Octave reproduce the SAME run.
% Delete sim_noise.csv to generate a fresh random realisation.
if exist('sim_noise.csv','file')
    Nz = dlmread('sim_noise.csv');  p.w = Nz(1:13,:);  p.vn = Nz(14:end,:);
else
    p.w  = Binv*(sigP.*randn(13,N));
    p.vn = randn(nm+1,N);
    dlmwrite('sim_noise.csv', [p.w; p.vn], 'precision', '%.12g');
end
Lb = [3 4 5 6 9 10 11 12 13 14];              % load buses used for demand estimate
p.Lb = Lb;

%% ---------- Attack definitions ----------
win = (k>=k_on & k<=k_off);
% (B) coordinated attack: 4 demand meters under-report by 30 %
atkbus = [3 4 9 14];
aB = zeros(nm,N);
for b = atkbus, aB(nl+b, win) = 0.30*Pd(b)/base; end
% (C) stealthy attack a = H*c : offset angle of bus 3 so that apparent demand falls by 20 MW
c = zeros(13,1);  c(3-1) = 0.20/bl(3);        % bl(3) = susceptance of line 2-3
aStep = zeros(nm,N);  aRamp = zeros(nm,N);
for kk = find(win)
    aStep(:,kk) = H*c;
    aRamp(:,kk) = H*c*(kk-k_on+1)/(k_off-k_on+1);
end
stealth_meters = find(abs(H*c) > 1e-9);

%% ---------- Runs ----------
Z0 = zeros(nm,N);
s0  = run_kf(p,H,R,Z0,false);                 % normal
sBs = run_kf(p,H,R,aB,false);                 % coordinated FDIA, standard KF
sBm = run_kf(p,H,R,aB,true);                  % coordinated FDIA, with mitigation
sCs = run_kf(p,H,R,aStep,false);              % stealth step, standard
sCm = run_kf(p,H,R,aStep,true);               % stealth step, mitigation
sRs = run_kf(p,H,R,aRamp,false);              % stealth ramp, standard
sRm = run_kf(p,H,R,aRamp,true);               % stealth ramp, mitigation
% (D) stealth ramp with one TRUSTED PMU at bus 3 (attacker cannot alter it)
Haug = [H; zeros(1,13)];  Haug(nm+1,3-1) = 1;
Raug = blkdiag(R, spmu^2);
aRampD = [aRamp; zeros(1,N)];
sRp = run_kf(p,Haug,Raug,aRampD,true,nm+1);   % trusted PMU row is never isolated

%% ---------- Report ----------
fprintf('\n=== IEEE 14-bus FDIA simulation (attack samples %d-%d, designed false-alarm prob. %.0e) ===\n',k_on,k_off,pfa);
fprintf('Total demand on load buses = %.1f MW\n\n', -base*sum(Hp(Lb,:)*p.x0));
fprintf('Normal operation: KF false alarms %d/%d, static-WLS false alarms %d/%d\n', sum(s0.alarm),N,sum(s0.alarmW),N);
fprintf('\n(B) Coordinated FDIA: demand meters at buses [%s] under-report by 30%% (%.1f MW total)\n', num2str(atkbus), 0.30*sum(Pd(atkbus)));
rep('  standard KF    ',sBs,win,k_on,k);
rep('  KF + mitigation',sBm,win,k_on,k);
fprintf('\n(C) Stealthy attack a=H*c on %d meters (apparent demand -20 MW)\n', numel(stealth_meters));
rep('  step, standard  ',sCs,win,k_on,k);
rep('  step, mitigation',sCm,win,k_on,k);
rep('  ramp, standard  ',sRs,win,k_on,k);
rep('  ramp, mitigation',sRm,win,k_on,k);
fprintf('\n(D) Ramp stealth attack + trusted PMU at bus 3\n');
rep('  ramp, PMU + mitigation',sRp,win,k_on,k);

fprintf('\nMinimum meters to compromise for a stealthy single-bus attack (nonzero rows of H column):\n');
card = zeros(1,13);
for j = 1:13, card(j) = sum(abs(H(:,j)) > 1e-9); fprintf('  bus %2d: %2d meters\n', j+1, card(j)); end

%% ---------- Single-meter residual test  |r(k)| > 3*sqrt(S)  (coordinated attack, standard KF) ----------
im = nl + 3;                                    % demand (injection) meter at bus 3 - attacked
io = 1;                                         % line-flow meter 1-2 - not attacked
tau_a = 3*sBs.Sd(im,:);  tau_o = 3*sBs.Sd(io,:);
al_a = abs(sBs.r(im,:)) > tau_a;  al_o = abs(sBs.r(io,:)) > tau_o;
pre = k < k_on;  post = k > k_off;
fprintf('\nSingle-meter test |r| > 3*sqrt(S), standard KF, coordinated attack:\n');
fprintf('  attacked demand meter (bus 3): alarms before %d/%d | during %d/%d | after %d/%d\n', sum(al_a(pre)),sum(pre),sum(al_a(win)),sum(win),sum(al_a(post)),sum(post));
fprintf('  unattacked flow meter (1-2)  : alarms before %d/%d | during %d/%d | after %d/%d\n', sum(al_o(pre)),sum(pre),sum(al_o(win)),sum(win),sum(al_o(post)),sum(post));

%% ---------- Consequence for control: dispatch follows the estimated demand ----------
% ASSUMED illustrative parameters (typical orders of magnitude, not from a specific source):
f0 = 50;  Hin = 5;  Rdroop = 0.05;  Ts_s = 1;           % Hz, inertia constant (s), 5 % droop, 1 s per sample
Ssys = sum(Pd);                                          % system size = total load (MW)
Minert = 2*Hin*Ssys/f0;                                  % MW.s/Hz
beta   = Ssys/(Rdroop*f0);                               % MW/Hz primary-response strength (damping neglected)
dfN = freq_sim(s0.D  - s0.Dh,  Minert, beta, Ts_s);      % normal
dfS = freq_sim(sBs.D - sBs.Dh, Minert, beta, Ts_s);      % attack, standard KF
dfM = freq_sim(sBm.D - sBm.Dh, Minert, beta, Ts_s);      % attack, with mitigation
fprintf('\nGeneration shortfall (true minus estimated demand) and frequency deviation, coordinated attack:\n');
fprintf('  worst shortfall (MW): normal %.1f | standard KF %.1f | mitigated %.1f\n', max(s0.D-s0.Dh), max(sBs.D-sBs.Dh), max(sBm.D-sBm.Dh));
fprintf('  lowest frequency deviation (Hz): normal %.3f | standard KF %.3f | mitigated %.3f   (M=%.1f MW.s/Hz, beta=%.1f MW/Hz, assumed)\n', min(dfN), min(dfS), min(dfM), Minert, beta);

fprintf('ramp mit, meters isolated after attack: %.1f\n', mean(sRm.nIso(k>k_off)));
%% ---------- Plots ----------
c1=[0 0.35 0.75]; c2=[0.85 0.15 0.1]; c3=[0 0.6 0.2]; c4=[0.5 0.2 0.7];
f1 = figure('Name','Demand'); 
plot(k,-base*sum(Hp(Lb,:)*s0.x,1),'k','LineWidth',1.8); hold on;
plot(k,sBs.Dh,'Color',c2,'LineWidth',1.3);
plot(k,sBm.Dh,'Color',c3,'LineWidth',1.3);
mark(k_on,k_off);
xlabel('Sample k'); ylabel('Total demand on load buses (MW)');
title('Coordinated FDIA: actual vs estimated demand (IEEE 14-bus)');
legend('Actual demand','Estimate: standard KF','Estimate: with mitigation','Location','southwest'); grid on;

f2 = figure('Name','Detection');
semilogy(k,sBs.J+1e-3,'Color',c2,'LineWidth',1.2); hold on;
semilogy(k,sBm.J+1e-3,'Color',c3,'LineWidth',1.2);
semilogy([1 N],[chi(nm) chi(nm)],'k--','LineWidth',1.2);
ylim([10 1e4]);
mark(k_on,k_off);
xlabel('Sample k'); ylabel('Detection statistic J(k) (log scale)');
title('Coordinated FDIA: J(k) = r^T S^{-1} r vs threshold \gamma');
legend('Standard KF','KF with mitigation',sprintf('\\gamma = %.1f (P_{fa}=10^{-4})',chi(nm)),'Location','northeast'); grid on;

f3 = figure('Name','StateError');
plot(k,s0.eRMS,'Color',c1,'LineWidth',1.2); hold on;
plot(k,sBs.eRMS,'Color',c2,'LineWidth',1.2);
plot(k,sBm.eRMS,'Color',c3,'LineWidth',1.2);
mark(k_on,k_off);
xlabel('Sample k'); ylabel('Angle estimation error, RMS over buses (deg)');
title('State-estimation accuracy under coordinated FDIA');
legend('Normal','Attack, standard KF','Attack, with mitigation','Location','northwest'); grid on;

f4 = figure('Name','Stealth');
semilogy(k,sCs.J/chi(nm)+1e-3,'Color',c2,'LineWidth',1.2); hold on;
semilogy(k,sRs.J/chi(nm)+1e-3,'Color',c1,'LineWidth',1.2);
semilogy(k,sRs.JW/sRs.chiW+1e-3,'Color',[0.5 0.5 0.5],'LineWidth',1.2);
semilogy(k,sRp.J/chi(nm+1)+1e-3,'Color',c3,'LineWidth',1.2);
semilogy([1 N],[1 1],'k--','LineWidth',1.2);
ylim([0.1 1e3]);
mark(k_on,k_off);
xlabel('Sample k'); ylabel('J(k) / threshold  (alarm if > 1)');
title('Stealthy attack a = Hc: detection statistic relative to threshold');
legend('Step attack, KF','Ramp attack, KF','Ramp attack, static WLS','Ramp attack, KF + trusted PMU','Threshold','Location','northwest'); grid on;

f5 = figure('Name','StealthImpact');
plot(k,sRs.Dh-sRs.D,'Color',c1,'LineWidth',1.3); hold on;
plot(k,sRm.Dh-sRm.D,'Color',c4,'LineWidth',1.3);
plot(k,sRp.Dh-sRp.D,'Color',c3,'LineWidth',1.3);
mark(k_on,k_off);
xlabel('Sample k'); ylabel('Demand estimation error (MW)  [negative = under-estimate]');
title('Stealth ramp attack: effect on estimated demand');
legend('Standard KF','KF + mitigation','KF + mitigation + trusted PMU','Location','southwest'); grid on;

f6 = figure('Name','Vulnerability');
bar(2:14,card,'FaceColor',c1);
xlabel('Bus whose estimated angle is offset'); ylabel('Meters to compromise (stealthy attack)');
title('Simplified measurement-coverage indicator for stealthy single-bus attacks'); grid on;

%% ---------- Export numerical results (used for the report figures) ----------
M = [k(:), s0.D(:), sBs.Dh(:), sBm.Dh(:), sBs.J(:), sBm.J(:), s0.eRMS(:), sBs.eRMS(:), sBm.eRMS(:), ...
     sCs.J(:), sRs.J(:), sRs.JW(:), sRp.J(:), sRs.Dh(:)-sRs.D(:), sRm.Dh(:)-sRm.D(:), sRp.Dh(:)-sRp.D(:), ...
     sCm.Dh(:)-sCm.D(:), sCs.Dh(:)-sCs.D(:)];
dlmwrite('sim_results.csv', M, 'precision', '%.6g');
dlmwrite('sim_constants.csv', [chi(nm) chi(nm+1) sRs.chiW card], 'precision', '%.6g');

f7 = figure('Name','Residual');
subplot(2,1,1);
kp = 4:N;                                      % skip first samples (initial-covariance transient)
plot(k(kp), base*sBs.r(im,kp), 'Color',c2, 'LineWidth',1.2); hold on;
plot(k(kp), base*tau_a(kp), 'k--'); plot(k(kp), -base*tau_a(kp), 'k--');
mark(k_on,k_off);
xlabel('Sample k'); ylabel('Residual r(k) (MW)');
title('Attacked demand meter (bus 3): r(k) = y(k) - yhat(k) with \pm3\surd S band');
legend('Residual r(k)','\pm3\surd S threshold','Location','southwest'); grid on;
subplot(2,1,2);
plot(k(kp), base*sBs.r(io,kp), 'Color',c1, 'LineWidth',1.2); hold on;
plot(k(kp), base*tau_o(kp), 'k--'); plot(k(kp), -base*tau_o(kp), 'k--');
mark(k_on,k_off);
xlabel('Sample k'); ylabel('Residual r(k) (MW)');
title('Unattacked flow meter 1-2: drifts above the band during the attack');
legend('Residual r(k)','\pm3\surd S threshold','Location','southwest'); grid on;

f8 = figure('Name','Frequency');
plot(k, dfN, 'Color',c1, 'LineWidth',1.3); hold on;
plot(k, dfS, 'Color',c2, 'LineWidth',1.3);
plot(k, dfM, 'Color',c3, 'LineWidth',1.3);
mark(k_on,k_off);
xlabel('Sample k'); ylabel('Frequency deviation \Delta f (Hz)');
title('Frequency deviation from the estimation error (assumed M, droop)');
legend('Normal','Attack, standard KF','Attack, with mitigation','Location','southwest'); grid on;

names = {'fig1_demand','fig2_detection','fig3_state_error','fig4_stealth_detection','fig5_stealth_impact','fig6_vulnerability','fig7_residual','fig8_frequency'};
fh = [f1 f2 f3 f4 f5 f6 f7 f8];
for i = 1:8
    try, print(fh(i),'-dpng','-r130',[names{i} '.png']); catch, end
end
end

%% =====================================================================
function s = run_kf(p,Hm,Rm,atk,mitig,trusted)
if nargin < 6, trusted = []; end   % indices of trusted meters that may never be isolated
N = p.N; nm = size(Hm,1); nx = 13;
sd  = sqrt(diag(Rm));
x   = p.x0;  xh = p.x0;  P = p.P0;
Ri  = inv(Rm);  Gw = (Hm'*Ri*Hm);
s.J=zeros(1,N); s.JW=zeros(1,N); s.alarm=false(1,N); s.alarmW=false(1,N);
s.r=zeros(nm,N); s.Sd=zeros(nm,N);
s.nIso=zeros(1,N); s.eRMS=zeros(1,N); s.D=zeros(1,N); s.Dh=zeros(1,N); s.x=zeros(nx,N);
dofW = nm - nx; s.chiW = 2*gammaincinv(1e-4, dofW/2, 'upper');
for k = 1:N
    y = Hm*x + sd.*p.vn(1:nm,k) + atk(:,k);
    % static WLS bad-data test on this snapshot only
    xw = Gw\(Hm'*Ri*y);  rw = y - Hm*xw;  s.JW(k) = rw'*Ri*rw;
    % Kalman residual from the a-priori prediction
    S = Hm*P*Hm' + Rm;  r = y - Hm*xh;
    U = true(nm,1);
    s.J(k) = r'*(S\r);
    if mitig
        for it = 1:20
            idx = find(U); SU = S(idx,idx); rU = r(idx);
            if rU'*(SU\rU) <= p.chi(numel(idx)), break; end
            z = abs(rU)./sqrt(diag(SU));  z(ismember(idx,trusted)) = 0;
            if max(z) == 0, break; end
            [~,j] = max(z);  U(idx(j)) = false;
        end
    end
    Hs = Hm(U,:);  Rs = Rm(U,U);
    K  = P*Hs'/(Hs*P*Hs' + Rs);
    xf = xh + K*(y(U) - Hs*xh);
    Pf = (eye(nx) - K*Hs)*P;  Pf = (Pf+Pf')/2;
    s.nIso(k) = sum(~U);
    s.r(:,k) = r;  s.Sd(:,k) = sqrt(diag(S));   % residual and its standard deviation
    s.x(:,k) = x;
    s.eRMS(k) = sqrt(mean((xf-x).^2))*180/pi;
    s.D(k)  = -p.base*sum(p.Hp(p.Lb,:)*x);
    s.Dh(k) = -p.base*sum(p.Hp(p.Lb,:)*xf);
    x  = x + p.Bu*p.u0 + p.w(:,k);
    xh = xf + p.Bu*p.u0;
    P  = Pf + p.Q;
end
s.alarm  = s.J  > p.chi(nm);
s.alarmW = s.JW > s.chiW;
end

function rep(name,s,win,k_on,k)
idx = find(s.alarm & k>=k_on,1);
if isempty(idx), dly = NaN; else dly = k(idx)-k_on; end
e = s.Dh(win)-s.D(win);
fprintf('%s: KF alarms in attack %d/%d (before: %d) delay %g | WLS alarms %d/%d | angle RMS err %.3f deg | demand err mean %.1f, worst %.1f MW | meters isolated (mean) %.1f | demand err at attack end %.1f MW\n', ...
    name,sum(s.alarm(win)),sum(win),sum(s.alarm(~win & k<k_on)),dly,sum(s.alarmW(win)),sum(win), ...
    sqrt(mean(s.eRMS(win).^2)), mean(e), e(find(abs(e)==max(abs(e)),1)), mean(s.nIso(win)), s.Dh(find(win,1,'last'))-s.D(find(win,1,'last')));
end

function df = freq_sim(imb, M, beta, Ts)
% First-order frequency response to a generation shortfall imb (MW, positive = short):
%   M*d(df)/dt = -imb - beta*df     (imb held constant over each sample)
% Exact discretisation (stable for any Ts): df(n+1) = a*df(n) - (1-a)*imb(n)/beta,  a = exp(-beta*Ts/M)
a = exp(-beta*Ts/M);
df = zeros(size(imb));
for n = 1:numel(imb)-1
    df(n+1) = a*df(n) - (1-a)*imb(n)/beta;
end
end

function mark(k_on,k_off)
yl = ylim; plot([k_on k_on],yl,'k:'); plot([k_off k_off],yl,'k:'); ylim(yl);
end
