%% Create data matrix
%
% Nishimura, Yuhei, Christa L. Martin, Araceli Vazquez-Lopez, Sarah J. 
% Spence, Ana Isabel Alvarez-Retuerto, Marian Sigman, Corinna Steindler et 
% al. "Genome-wide expression profiling of lymphoblastoid cell lines 
% distinguishes different forms of autism and reveals shared pathways." 
% Human molecular genetics 16, no. 14 (2007): 1682-1698.
%
% The microarray expression data from this study have been submitted to GEO 
% under accession number GSE7329.
%
% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7329
%
%
% LogRatio extracted from variable 17
%
% To Run this script, download first https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7329
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

files = gunzip('*.gz');

annotations = readtable('GPL1708_old_annotations.txt');
ind = find(annotations{:,6}=="FALSE");
ind2 = annotations{ind,1};

X = zeros(length(files)-1,length(ind2)); 
for i = 2:length(files)
    data = readtable(files{i});
    [~,ia,ib] = intersect(data{:,2},ind2);
    X(i-1,ib) = data{ia,17};
    obs_l{i-1} = files{i}(1:end-4);
end

var_l = annotations(ind2,4);
class = repmat("dup(15q)",7,1);
class(8:15) = repmat("FMR1-FM",8,1);
class(16:30) = repmat("control",15,1);

save GSE7329_RAW X obs_l var_l class

delete *.txt


