function [hess, massMat] = sbNMA_PSF(xyz, bonds, angles, torsional, improper, atomTypes, charmmBonds, charmmAngles, charmmPhis, charmmImps, charmmNb, charmmMass)
% sbNMA normal modes
% Author: Guang song
% date: 08/03/17
% last update: 01/16/2018. set an upperbound for K_vdW to be 10.
% xyz: coordinates of the protein
% bonds: nb-by-4 matrix, each row of which has 2 ints that represent indices of 2 interacting atoms
% angles: na-by-4 matrix, each row of which has 3 ints that represent indices of 3 interacting atoms
% torsional: nt=by-4 matrix, each row of which has 4 ints that represent indices of 4 interacting atoms
% improper: ni-by-4 matrix, each row of which has 4 ints that represent indices of 4 interacting atoms
% atomTypes: a n-by-1 matrix where n is the number of atoms in this protein, representing the charmm atom type of each atom in this protein.  
% charmm*: charmm parameters
%
% What to cite: 
% Hyuntae Na and Guang Song, “Bridging between NMA and Elastic Network Models”, Proteins, 2014, 82(9):2157-68, 2014
%
% Copyright (c) 2020 Guang Song. All Rights Reserved. 
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% pairwise distance among all the atoms
dist = squareform(pdist(xyz));
n = size(xyz,1);

ub = angles(:, [1,3]);
ub =[min(ub,[], 2), max(ub,[],2)];

c14 = []; % 1-4 contacts
bondMat = sparse(bonds(:,1),bonds(:,2),1,n,n);
bondMat = bondMat + bondMat';
for i=1:size(bonds,1)
	b1 = bonds(i,1);
	b2 = bonds(i,2);
	c2b1 = setdiff(find(bondMat(:,b1)>0), b2); % c2b1: atoms connect 2 b1, excluding b2
	c2b2 = setdiff(find(bondMat(:,b2)>0), b1); % c2b2: atoms connect 2 b2, excluding b1
	if (~isempty(c2b1) && ~isempty(c2b2)) 
		c14 = [c14; combine(c2b1(:), c2b2(:))];
	end
end
% sort the order of each 1-4 contact
c14 =[min(c14,[], 2), max(c14,[],2)];
c14 = setdiff(c14,ub,'rows'); % remove 1-3 contacts from it.

% check if c14 is computed correctly
for i=1:size(c14, 1)
   if bondMat(c14(i,1), c14(i,2))~=0
	  display('find a bond');
   elseif length(find(bondMat(:, c14(i,1)) & bondMat(:, c14(i,2))))~=0
	  display('find an angle');
   end
end

bonds = sort(bonds,2);
% vc: vdW contacts. 12: cutoff dist.
vc = dist <= 12;
[vi, vj] = find(triu(vc,1));
vc = [vi, vj];
vc = setdiff(vc, bonds, 'rows'); % remove 1-2 contacts
vc = setdiff(vc, ub, 'rows'); % remove 1-3 contacts
vc = setdiff(vc, c14, 'rows');% remove 1-4 contacts


% determine the vdW contacts with spring constants
% k_vdW = 12epsilon/r^2(13(r_0/r)^12 - 7(r0/r)^6)
vAll = [vc; c14];  % combine vc and c14 
[epsilon, r0] = computeKvdW(vc, c14, atomTypes, charmmNb);
% r0: vdW radii sums. r_0's for vc and c14 are different
r = dist(sub2ind(size(dist), vAll(:,1), vAll(:,2)));
r6 = (r0./r).^6;
k_vdW = 12*epsilon.*(13*r6.^2-7*r6)./(r.^2);
k_vdW(find(k_vdW<0)) = 0; % zero all the negatives
k_vdW(find(k_vdW>10)) = 10; % 01/16/2018 set an upper limit on k_vdW



[kb] = computeKb(bonds, atomTypes, charmmBonds);
[ktheta] = computeKtheta(angles,atomTypes, charmmAngles);
[kphi] = computeKphi(torsional, atomTypes, charmmPhis, [1,4]);
[kimp] = computeKphi(improper,  atomTypes, charmmImps, [2,3]);

% contact matrix for UB and vdW terms
k_UB = ktheta(:,2);
cxUB = sparse(ub(:,1),ub(:,2),k_UB*2,n,n); % *2 because non-standard k_UB*(s-s0)^2
cx_vdW = sparse(vAll(:,1), vAll(:,2), k_vdW, n, n);
% now determine the contact matrix of all 2-body interactions
cx_bond = sparse(bonds(:,1),bonds(:,2),kb*2,n,n);% *2 because non-standard kb*(b-b0)^2 form
cx12 = (cx_bond + cx_bond') + (cxUB + cxUB') + (cx_vdW + cx_vdW');



hess=baseHess(xyz, cx12);
hess=hess+STeMtheta(xyz,angles,ktheta(:,1));
n = kphi(:,6); % mulplicity
hess=hess+STeMphi(xyz,kphi(:,1:4),abs(kphi(:,5)).*n.^2/2); % apply abs since some kphi may be negative in charmm
hess=hess+STeMphi(xyz,kimp(:,1:4),kimp(:,5));

%hess = (hess+hess')/2; % make it fully symmetric
id3 = mat2vec([atomTypes, atomTypes, atomTypes]);
massMat = charmmMass(id3+1); % mass matrix, +1: make it 1-based
%[V, D] = eig(hess, massMat); % generalize eigenvalues
%dia = diag(D);


function [epsilon, r0] = computeKvdW(vc, c14, atomTypes, charmmNb)
charmmNb = charmmNb(find(charmmNb(:,1)>-1),:);
charmmNb = sortrows(charmmNb,1);
vc = atomTypes(vc) + 1; % make it 1-based
epsilon = sqrt(charmmNb(vc(:,1), 2).*charmmNb(vc(:,2), 2));
r0 = charmmNb(vc(:,1), 3) + charmmNb(vc(:,2), 3);
c14 = atomTypes(c14) + 1;
epsilon14 = sqrt(charmmNb(c14(:,1), 4).*charmmNb(c14(:,2), 4));
r0_14 = charmmNb(c14(:,1), 5) + charmmNb(c14(:,2), 5);
epsilon = [epsilon; epsilon14];
r0 = [r0; r0_14];


function [kb] = computeKb(bonds, atomTypes, charmmBonds)
charmmBonds = rearrage(charmmBonds, 1,2, 1:2);
bonds = atomTypes(bonds);
bonds = rearrage(bonds, 1,2, 1:2);
kb = zeros(size(bonds,1), 1);
for i=1:size(bonds,1)
  row = findlines(charmmBonds(:,1:2),bonds(i,:));
  if length(row)~=1
	display(i);
	display(row);
	error('matching row is not unique.');
  end
  kb(i) = charmmBonds(row,3);
end

function [ktheta] = computeKtheta(angles, atomTypes, charmmAngles)
charmmAngles = rearrage(charmmAngles, 1,3, 1:3);
angles = atomTypes(angles);
angles = rearrage(angles, 1, 3, 1:3);
ktheta = zeros(size(angles,1), 2);
for i=1:size(angles,1)
  row = findlines(charmmAngles(:,1:3),angles(i,:));
  if length(row)~=1
	display(i);
	display(row);
	error('matching row is not unique.');
  end
  ktheta(i,:) = charmmAngles(row,[4,6]);
end

function [kphi] = computeKphi(torsional, atomTypes, charmmPhis, wildcardCol)
charmmPhis = rearrage(charmmPhis, 2,3, 1:4);
phis = atomTypes(torsional);
phis = rearrage(phis, 2,3, 1:4);
col = setdiff(1:4, wildcardCol);
kphi = []; % kphi may have more rows than phis since some torsional 
		   % have two k_phi's with different mulplicity

for i=1:size(phis,1) 
  row = findlines(charmmPhis(:,1:4),phis(i,:));
  if length(row) == 0
     idx = find(charmmPhis(:,wildcardCol(1))==-1);
	 row = findlines(charmmPhis(idx,col), sort(phis(i, col),2));
	 row = idx(row);
	 if length(row)==0
		display(phis(i,:));
		error('no matching row is found.');
	 end
  end
  kphi = [kphi; torsional(i*ones(length(row),1),:), charmmPhis(row,5:end)];
end


function [A] = rearrage(A, i, j, range)
idx = find(A(:,i)>A(:,j));
if length(range)==4 % in this case, i=2, j=3
   idx2 = find(A(:,i)==A(:,j) & A(:,1)>A(:,4));
   idx = [idx; idx2];
end
A(idx, range) = A(idx, range(end:-1:1)); % rearrange of the colomns of these rows by reversing.


function [set] = combine(v1, v2)
n1 = length(v1);
n2 = length(v2);
set = [mat2vec(repmat(v1,1,n2))', repmat(v2,n1,1)];