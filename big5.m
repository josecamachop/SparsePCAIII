%% Big5 from R package qgraph "All Sparse PCA Models Are Wrong, But Some Are Useful. Part III: model interpretation" 
% Submitted to Chemometrics and Intelligent Laboratory Systems. 2025
%
% Dependencies: MEDA Toolbox v1.8 at https://github.com/codaslab/MEDA-Toolbox
%               SPASM Toolbox at https://www2.imm.dtu.dk/projects/spasm
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 10/Apr/2025
%
% Copyright (C) 2025  University of Granada, Granada
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

%% Big5 data from R package qgraph

data = importdata('big5.csv');

var_l = data.textdata;
X = data.data;

M = size(X,2);
reord = [1:5:M 2:5:M 3:5:M 4:5:M 5:5:M];
X = X(:,reord);
var_l = var_l(reord);

var_class(1:48) = repmat("Neuroticism",1,48); % Variable classes
var_class(48 + (1:48)) = repmat("Extraversion",1,48);
var_class(2*48 + (1:48)) = repmat("Openness",1,48);
var_class(3*48 + (1:48)) = repmat("Agreeableness",1,48);
var_class(4*48 + (1:48)) = repmat("Conscientiousness",1,48);

pcs = 1:5;


%% Compute variance estimates for PCA with 5 PCs

X = preprocess2D(X); % autoscale
model = pcaEig(X,'PCs',pcs);

model.lvs = 1:2;
biplot(model)

% Compute variance estimates in PCA
model.lvs = 1:5;
Xest = model.scores*model.loads';
totVPCA = 100*sum(sum(Xest.^2))/trace(X'*X)

%% PEV vs sparsity: SPCA-Z multi-component, truncated search

clc
PEVpq = [];
fp = [];
D = min(size(X));
nze = [1:24:240 240];
ridge = [0 1 10 100 10000 Inf];
thres = 0.05
flag = 0;
tic
for j1=1:length(nze)
    for j2= 1:j1
        for j3= 1:j2
            for j4= 1:j3
                for j5= 1:j4
                    vec = nze([j1 j2 j3 j4 j5])
                    for j=1:length(ridge)
                        for i=1:length(pcs)
                            p = spca_zouhastie(X, [], pcs(i), ridge(j), -vec(1:pcs(i)));
                            [u,s,v]=svd(X'*X*p,0);
                            q=u*v';

                            fp(i,j,j1,j2,j3,j4,j5) = length(find(p.^2)) - length(find(sum(p.^2,1)));
                            PEVpq(i,j,j1,j2,j3,j4,j5) = 1 - sum(sum((X - X*p*inv(q'*p)*q').^2))/sum(sum(X.^2));

                            if (100*PEVpq(i,j,j1,j2,j3,j4,j5)/totVPCA) > (1 - thres)
                                flag = 1;
                                break
                            end
                        end
                    end
                    if flag, break; end
                end
                if flag, break; end
            end
            if flag, break; end
        end
        if flag, break; end
    end
    if flag, break; end
end
total_time=toc

close all


%% Plot the truncated razor plot

ufp = unique(fp);
PEVfp = [];
for i=1:length(ufp)
    ind = find(fp == ufp(i));
    mind = find(PEVpq(ind)==max(PEVpq(ind)),1);
    PEVfp(i) = PEVpq(ind(mind));
end

val = num2cell(ufp);
val{end+1} = 'Ref';
f = plotVec([PEVfp totVPCA/100],'ObsClass',[2*ones(1,length(PEVfp)) 1]);
legend('off')
ylabel('PEV')
xlabel('f')
a=get(f,'Children');
set(a,'XTickLabel',val);
set(a,'XTick',1:length(val));
set(a,'XTickLabelRotation',45);

saveas(gcf,'Figures/razorBig5');
saveas(gcf,'Figures/razorBig5.eps','epsc');

    
%% Plot the truncated razor plot per component

ufp = unique(fp);
PEVfp2D = [];
for i=1:length(ufp)
    for j = pcs
        ind = find(fp(j,:) == ufp(i));
        mind = find(PEVpq(j,ind)==max(PEVpq(j,ind)),1);
        if ~isempty(mind)
            PEVfp2D(j,i) = PEVpq(j,ind(mind));
        end
    end
end

figure
surf((((ones(length(pcs),1)*ufp')))',(pcs'*ones(1,length(ufp)))',(PEVfp2D)')
hold on
pcolor((((ones(length(pcs),1)*ufp')))',(pcs'*ones(1,length(ufp)))',(PEVfp2D)')
axis([ufp(1) ufp(end) pcs(1) pcs(end)])
colorbar
ylabel('# Components')
xlabel('f')
zlabel('PEV')
saveas(gcf,'Figures/surfaceBig5');
saveas(gcf,'Figures/surfaceBig5.eps','epsc');


%% Select ridge penalty: 10

f=plotVec(PEVpq(5,:,3,3,3,3,2));
legend('off')
ylabel('PEV')
xlabel('Ridge')
a=get(f,'Children');
set(a,'XTickLabel',ridge);
saveas(gcf,'Figures/ridgeBig5');
saveas(gcf,'Figures/ridgeBig5.eps','epsc');


%% Visualize multi-model selected with 5 components with two non-zero weights each.

p = spca_zouhastie(X, [], max(pcs), 10, -vec);

[u,s,v]=svd(X'*X*p,0);
q=u*v';
r = q*inv(p'*q);

for i=pcs
    plotVec(p(:,i),'XYLabel',{'Variables',sprintf('Sparse weights (p_%d)',i)},'ObsClass',var_class);
    legend('Location','southeast');
    if i>1, legend('off'); end
    saveas(gcf,sprintf('Figures/pBig5%d',i));
    saveas(gcf,sprintf('Figures/pBig5%d.eps',i),'epsc');

    for j=1:5
        restab(j,i) = sum(p((j-1)*48+1:j*48,i)~=0);
    end

    plotVec(q(:,i),'XYLabel',{'Variables',sprintf('Auxiliary loadings (q_%d)',i)},'ObsClass',var_class);
    if i>1, legend('off'); end
    saveas(gcf,sprintf('Figures/qBig5%d',i));
    saveas(gcf,sprintf('Figures/qBig5%d.eps',i),'epsc');

    f = plotMap(1/(size(X,1)-1)*[r(:,i)*p(:,i)'*X'*X*p(:,i)*r(:,i)']);
    ylabel('Variables','Fontsize',16)
    xlabel('Variables','Fontsize',16)
    saveas(gcf,sprintf('Figures/mapBig5_%d',i));
    saveas(gcf,sprintf('Figures/mapBig5_%d.eps',i),'epsc');
end
restab(6,:) = sum(restab);
rows = unique(var_class);
rows(end+1) = 'Total nonzero';
T = table(restab(:,1),restab(:,4),restab(:,3),restab(:,4),restab(:,5),'VariableNames',{'w_1','w_2','w_3','w_4','w_5'},'RowNames',rows)

f = plotMap(1/(size(X,1)-1)*[r*p'*X'*X*p*r'],'VarsLabel',var_l);
a = get(f,'Children');
set(a(2),'XTickLabelRotation',45);
saveas(gcf,'Figures/mapBig5_total');
saveas(gcf,'Figures/mapBig5_total.eps','epsc'); 

