function [ P_next, V_next] = computeNext( P, V)
%  The numerical iteration function
%   INPUT: Bubbles' last Velocity & Position
%   OUTPUT: Bubbles' next Velocity & Position

%% DECLARE PRE-SET GLOBAL PARAMETERS
global F_g_set; % gravity (1 x 3)
global k_v;
global k_of;
global k_r;
global k_a;
global k_air;
global m; % number of bubbles
global R; % radii (m x 1)
global PO;
global RO;
global mo;
global k_ro;
global k_ao;

%% MATRIX VARIABLE DESCRIPTION
% NB ------- Number of bubbles in contact(m x 1)
% NO ------- Number of objecs in contact (m x 1)
% DB ------- Distances with bubbles      (m x m)
% DO ------- Distances with objecs       (m x m)
% R1R2 ----- Sums of two radii           (m x m)
% diffB ---- Difference of DB & R1R2     (m x m)
% diffO ---- Difference of DO & R        (m x m)
% FLAGB ---- Flags that shows two bubbles are overlapped(1) or not(0)             (m x m)
% FLAGO ---- Flags that shows a bubble and an object are overlapped(1) or not(0)  (m x m)
% V_meanB -- Mean value of the overlapped bubbles' velocities (m x 3)
% V_meanO -- Mean value of the overlapped objects' velocities (m x 3)
% F_ra ----- All repulsive and attractive forces    (m x 3)
% F_rB ----- Repulsive force by overlapped bubbles  (3m x m)
% F_aB ----- Attractive force by overlapped bubbles (3m x m)
% F_rO ----- Repulsive force by overlapped objects  (3m x m)
% F_aO ----- Attractive force by overlapped objects (3m x m)

%% start to compute
% bubbles with bubbles
DB = pdist2(P,P); % compute distances: DB
R1R2 = bsxfun(@plus,R,R.');% compute sums of two radii
diffB = DB - R1R2;% compute difference of DB & R1R2 (m x m)
FLAGB = diffB<0; % compute overlap flags for bubbles
FLAGB(1:m+1:end)=0;
NB = sum(FLAGB, 2); % compute overlap nubmbers for bubbles

% bubbles with objects
DO = pdist2(P,PO); % compute distances: DB
R1RO = bsxfun(@plus,R,RO.');% compute sums of two radii
diffO = DO - R1RO;% compute difference of DB & R1R2 (m x m)
FLAGO = diffO<0; % compute overlap flags for bubbles
NO = sum(FLAGO, 2); % compute overlap nubmbers for bubbles

%% compute V_meanB
NB_3D = kron(NB, ones(3,1)); % extend velocity matrix for manipulation (3m x 1);
V_3D = repmat(V.', m, 1); % extend velocity matrix for manipulation (3m x m);
FLAGB_3D = kron(FLAGB, ones(3,1)); % extend FLAGB matrix for manipulation (3m x m);
V_overlappedB = V_3D .* FLAGB_3D; % keep the overlapped bubbles' velocities only (3m x m);

V_meanB = sum(V_overlappedB,2)./NB_3D; % (3m x 1)
V_meanB = reshape(V_meanB, 3, m)'; % (m x 3)
V_meanB(isnan(V_meanB)) = 0; % fix NAN

%% compute V_meanO
V_meanO = [0,0,0];
V_meanO = repmat(V_meanO, m, 1); % (m x 3)

FLAGO_3D = kron(FLAGO, ones(3,1)); % extend FLAGB matrix for manipulation (3m x m);


%% compute F_ra
DB_3D = kron(DB, ones(3,1)); % extend distance matrix for manipulation (3m x m);
R1R2_3D = kron(R1R2, ones(3,1)); % extend R1R2 matrix for manipulation (3m x m);
DO_3D = kron(DO, ones(3,1)); % extend distance matrix for manipulation (3m x m);
R1RO_3D = kron(R1RO, ones(3,1)); % extend R1RO matrix for manipulation (3m x m);

% compute F_rB
P_trans = P.';
PP_trans = repmat(P_trans,m,1);
Disp = bsxfun(@minus,P_trans(:),PP_trans);% compute "displacements" (3m x m)
F_rB = k_r .* (1./DB_3D - 1./R1R2_3D) .* Disp; % compute F_rB  (3m x m)
F_rB = F_rB .* FLAGB_3D; % keep the overlapped bubbles' forces only (3m x m);
F_rB(isnan(F_rB)) = 0; % fix NAN

% compute F_rO
PO_trans = PO.';
PO_trans = repmat(PO_trans,m,1);
DispO = bsxfun(@minus,P_trans(:),PO_trans);% compute "displacements" (3m x mo)

F_rO = k_ro .* (1./DO_3D - 1./R1RO_3D) .* DispO; % compute F_rB  (3m x mo)
% F_rO = k_ro .* (R1RO_3D - DispO); % compute F_rB

F_rO = F_rO .* FLAGO_3D; % keep the overlapped bubbles & objects forces only (3m x mo);
F_rO(isnan(F_rO)) = 0; % fix NAN

% compute F_aB
% (1)C_nb
NB_inv = 1./NB; % (m x 1)
NB_inv(isinf(NB_inv)) = 0;
NB_temp = bsxfun(@plus,NB_inv,NB_inv.');
NB_temp_3D = kron(NB_temp, ones(3,1)); % (3m x m)
C_nb = 1/2 .* NB_temp_3D;
% (2)C_dist
maxR = max(repmat(R,1,m), repmat(R',m,1));
minR = min(repmat(R,1,m), repmat(R',m,1));
maxR_3D = kron(maxR, ones(3,1)); % (3m x m)
minR_3D = kron(minR, ones(3,1)); % (3m x m)
C_dist = (DB_3D - maxR_3D)./minR_3D;
% (3)F_aB
F_aB = k_a .* C_nb .* C_dist .* (-1.*Disp) ./ DB_3D;% compute F_aB  (3m x m)
F_aB(isnan(F_aB)) = 0; % fix NAN
F_aB = F_aB .* FLAGB_3D; % keep the overlapped bubbles' forces only (3m x m);

% compute F_aO
% (1)C_nbO
NO_inv = 1./NO;
NO_inv(isinf(NO_inv)) = 0;
NO_temp = bsxfun(@plus,NB_inv,NO_inv);
NO_temp_3D = kron(NO_temp, ones(3,mo)); % (3m x mo)
C_nbO = 1/2 .* NO_temp_3D;
% (2)C_distO
C_distO = ones(3*m,mo);
% (3)F_aO
% F_aO = k_ao .* C_nbO .* C_distO .* (-1.*DispO) ./ DO_3D;% compute F_aO  (3m x mo)
% ABOVE is the RIGHT FORMULA, however, strangely, the following one WORKS BETTER!
F_aO = k_ao .* C_nbO .* C_distO .* (DispO) ./ DO_3D;
F_aO(isnan(F_aO)) = 0; % fix NAN
F_aO = F_aO .* FLAGO_3D; % keep the overlapped bubbles & objects forces only (3m x mo);

% Final: compute F_ra
F_ra = sum(F_rB,2) + sum(F_aB,2) + sum(F_rO,2) + sum(F_aO,2); % (3m x 1)
F_ra = reshape(F_ra, 3, m)'; % (m x 3)

%% compute F_g
F_g = repmat(F_g_set, m, 1); % (m x 3)

%% compute next Velocity & Position
V_next = 1/(k_v + k_of +k_air) .* (k_v.*V_meanB + k_of.*V_meanO + F_ra + F_g);
P_next = P + V_next;

end
