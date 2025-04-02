%% Pitprops from R package elsticnet "All Sparse PCA Models Are Wrong, But Some Are Useful. Part III: model interpretation" 
% Submitted to Chemometrics and Intelligent Laboratory Systems. 2025
%
% Dependencies: MEDA Toolbox v1.8 at https://github.com/codaslab/MEDA-Toolbox
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 31/Mar/2025
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

%% Pitprops data from R package elsticnet

var_l = {'topdiam' 'length'  'moist' 'testsg' 'ovensg' 'ringtop' 'ringbut' 'bowmax' 'bowdist' 'whorls'  'clear'  'knots' 'diaknot'};

XX=[ 1.000  0.954  0.364  0.342 -0.129   0.313   0.496  0.424   0.592  0.545  0.084 -0.019   0.134
     0.954  1.000  0.297  0.284 -0.118   0.291   0.503  0.419   0.648  0.569  0.076 -0.036   0.144
     0.364  0.297  1.000  0.882 -0.148   0.153  -0.029 -0.054   0.125 -0.081  0.162  0.220   0.126
     0.342  0.284  0.882  1.000  0.220   0.381   0.174 -0.059   0.137 -0.014  0.097  0.169   0.015
    -0.129 -0.118 -0.148  0.220  1.000   0.364   0.296  0.004  -0.039  0.037 -0.091 -0.145  -0.208
     0.313  0.291  0.153  0.381  0.364   1.000   0.813  0.090   0.211  0.274 -0.036  0.024  -0.329
     0.496  0.503 -0.029  0.174  0.296   0.813   1.000  0.372   0.465  0.679 -0.113 -0.232  -0.424
     0.424  0.419 -0.054 -0.059  0.004   0.090   0.372  1.000   0.482  0.557  0.061 -0.357  -0.202
     0.592  0.648  0.125  0.137 -0.039   0.211   0.465  0.482   1.000  0.526  0.085 -0.127  -0.076
     0.545  0.569 -0.081 -0.014  0.037   0.274   0.679  0.557   0.526  1.000 -0.319 -0.368  -0.291
     0.084  0.076  0.162  0.097 -0.091  -0.036  -0.113  0.061   0.085 -0.319  1.000  0.029   0.007
    -0.019 -0.036  0.220  0.169 -0.145   0.024  -0.232 -0.357  -0.127 -0.368  0.029  1.000   0.184
     0.134  0.144  0.126  0.015 -0.208  -0.329  -0.424 -0.202  -0.076 -0.291  0.007  0.184   1.000];

clc


%% PEV vs sparsity: SPCA-Z multi-component, truncated search

[PEVpq, fp] = razorPlot([], XX, 6);


%% Visualize multi-model selected with 6 components with two non-zero weights each.

model = spcaZou([], XX, 6, -[3,2,2,2,1,1]);
p = model.weights;
q = model.loads;
r = model.altweights;

for i=pcs
    plotVec(-p(:,i),'XYLabel',{'Variables','Sparse weights (p)'});
    axis([.5 13.5 -1 1])
    
    plotVec(-q(:,i),'XYLabel',{'Variables','Sparse weights (p)'});
    axis([.5 13.5 -1 1])
    
    f = plotMap([r(:,i)*p(:,i)'*X'*X*p(:,i)*r(:,i)']);
    ylabel('Variables','Fontsize',16)
    xlabel('Variables','Fontsize',16);
end

f = plotMap([r*p'*X'*X*p*r'],'VarsLabel',var_l);
a = get(f,'Children');
set(a(2),'XTickLabelRotation',45); 

