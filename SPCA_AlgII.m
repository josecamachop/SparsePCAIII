% Alg II from Guerra-Urzola et al. psychometrika 86.4 (2021): 893-919. 
% Results for SPCA of Zou et al. consisten with the matching sparsity,
% Fig 1. We chose the situation in which SPCA showed the worst performance
% in comparison to GPower in the paper, namely:
%
% Number of PCs: 2
% Number of observations (I): 500
% Number of variables (J): 1000
% Sparsity: 80%
% VAF: 80%
%
% Requierements: MEDA Toolbox v1.3
%
% coded by: Torfinn Støve Madssen (torfinn.s.madssen@ntnu.no)
%       Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/March/2025
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Data simulation

clear
close all
clc

sre = @(A,B) (norm(A-B,'fro')^2)/(norm(A,'fro')^2);
pev = @(A,B) 1 - (norm(A - B,'fro')^2)/(norm(A,'fro')^2);

pcs = 1:2
obs = 500
vars = 1000
sparsity = 80
s2n = 80
    
ridge = Inf

nonzero = round(vars*(1-sparsity/100));

X = randn(obs,vars);
[U,D,V] = svd(X,0);
U = U(:,pcs);
D = D(pcs,pcs);
V = V(:,pcs);

for i = pcs % sparsize
    kk = sort(abs(V(:,i)),'descend');
    V(find(abs(V(:,i))<kk(nonzero)),i) = 0;
    V(:,i) = V(:,i)/norm(V(:,i));
end

T = X*V;
P = ((T'*T)\T'*X)';

X = T*P'; 
E = randn(size(X));
f = sqrt((1-s2n/100)*(norm(X,'fro')^2)/((s2n/100)*(norm(E,'fro')^2)));
X = X + f*E;


%% SPCA simultaneous uncorrected (SPCA I)

vec = nonzero*ones(size(pcs));
p7 = spca_zouhastie(X, [], max(pcs), ridge, -vec);

t7 = X*p7;

sig = sign(diag(P'*p7));
sig(find(sig==0)) = 1;
p7u = repmat(sig',vars,1).*p7; % sign correction
t7u = repmat(sig',obs,1).*t7;

SRE_LW_u = sre(V,p7u)
SRE_S_u = sre(T,t7u)
PEV_u = pev(X,t7u*p7u')

%% SPCA simultaneous corrected only with P (SPCA II)

vec = nonzero*ones(size(pcs));
p7 = spca_zouhastie(X, [], max(pcs), ridge, -vec);

t7 = X*p7*pinv(p7'*p7);

sig = sign(diag(P'*p7));
sig(find(sig==0)) = 1;
p7c = repmat(sig',vars,1).*p7; % sign correction
t7c = repmat(sig',obs,1).*t7;

SRE_LW_c = sre(V,p7c)
SRE_S_c = sre(T,t7c)
PEV_c = pev(X,t7c*p7c') 


%% SPCA simultaneous properly corrected (SPCA III)

vec = nonzero*ones(size(pcs));
p7 = spca_zouhastie(X, [], max(pcs), ridge, -vec);
[u,s,v]=svd(X'*X*p7,0);
q7=u*v';

t7 = X*p7*pinv(q7'*p7);

sig = sign(diag(P'*p7));
sig(find(sig==0)) = 1;
q7 = repmat(sig',vars,1).*q7; % sign correction
t7q = repmat(sig',obs,1).*t7;

SRE_LW_u = sre(V,p7u)
SRE_S_u = sre(T,t7q)
PEV_u = pev(X,t7q*q7')


%% PLOTs

f = figure; subplot(2,1,1), hold on, plot(T(:,pcs),'o-'), axis tight
a=get(f,'Children'); set(a,'FontSize',20)
plot([1 obs],[0 0],'k--'),
title('Simulated','FontSize',20), subplot(2,1,2), hold on, bar(P(:,pcs)),axis tight 
a=get(f,'Children'); set(a(1),'FontSize',20), axis tight

f = figure; subplot(2,1,1), hold on, plot(t7u(:,pcs),'o-'), axis tight 
a=get(f,'Children'); set(a(1),'FontSize',20)
plot([1 obs],[0 0],'k--'),
title('SPCA uncorrected','FontSize',20), subplot(2,1,2), hold on, bar(p7u(:,pcs)),axis tight 
a=get(f,'Children'); set(a(1),'FontSize',20), axis tight

f = figure; subplot(2,1,1), hold on, plot(t7c(:,pcs),'o-'), axis tight 
a=get(f,'Children'); set(a(1),'FontSize',20)
plot([1 obs],[0 0],'k--'),
title('SPCA corrected','FontSize',20), subplot(2,1,2), hold on, bar(p7c(:,pcs)),axis tight 
a=get(f,'Children'); set(a(1),'FontSize',20), axis tight

