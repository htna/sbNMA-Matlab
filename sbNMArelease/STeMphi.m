function [hessian] =STeMphi(x,torsional,kphi)
% a new implmentation of SbNMA torsional/improper terms. 
% Date: 03/07/2018
% Guang Song
% derive the hessian of the second term
% reference: Blondel and Karplus, NEW FORMULATION FOR ELIMINATION OF StNGULARITIES, Journal of Computational Chemistry, Vol. 17, No. 9, 1132-1141 (1996)
% x: xyz coordinates of protein
% torsional: a n-by-4 matrix with each row of 4 ints representing 0-based indices of the four atoms that have torsional or improper interactions. 
% kphi: the spring constant of each torsional/improper interaction.
% hessian:  hessian matrix due to torsional/improper interactions.
%
% Note: all eq. numbers are from the reference. 
% Note: at present, we are NOT consider d2theta/dr2 term but only outer product of (dphi/dr) with itself. The derivation of d2theta/dr2 term is included before but commented out. Its correctness has NOT been checked.
% Copyright © 2020 Guang Song. All Rights Reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
if length(kphi)==1
   kphi = kphi(ones(size(torsional,1),1)); % duplicate it for each dihedral 09/27/17
end
		
hessian=zeros(size(x,1)*3, size(x,1)*3);
for m=1:size(torsional,1)
	K_phi = kphi(m); % added on 09/27/17 for sbNMA
    r=x(torsional(m,:),1:3);
    F = r(1,:)-r(2,:);
	G = r(2,:)-r(3,:);
	H = r(4,:)-r(3,:);
	A = cross(F,G);
	B = cross(H,G);
	A2 = sum(A.^2);
	B2 = sum(B.^2);
	G2 = sum(G.^2);
	dphidFGH = [ -norm(G)/A2*A;	
			dot(F,G)/A2/norm(G)*A - dot(H,G)/B2/norm(G)*B;
			norm(G)/B2*B]'; % eq. 21, 25, 22
	dFGHdr = [1 -1 0 0; 
			  0 1 -1 0; 
			  0 0 -1 1]; % eq. (26)
	dphiDr = dphidFGH*dFGHdr;
	dphi2 = dphiDr(:)*dphiDr(:)';
	
	% d2dFG =  1/norm(G)/A2^2*(G2*outer(cross(A,F),A)+dot(F,G)*outer(A,cross(A,G)));
	% d2dGH = -1/norm(G)/B2^2*(G2*outer(cross(B,H),B)+ dot(H,G)*outer(B,cross(B,G)));
	% d2dFH = 0;
	% d2dF2 = norm(G)/A2^2*outer(A,cross(G,A));
	% d2dF2 = d2dF2 + d2dF2';
	% d2dG2 = 1/2/norm(G)^3/A2*outer(cross(G,A),A) + dot(F,G)/norm(G)/A2^2*outer(A,cross(F,A)) - 1/2/norm(G)^3/B2*outer(cross(G,B),B) - dot(H,G)/norm(G)/B2^2*outer(B,cross(H,B));
	% d2dG2 = d2dG2 + d2dG2'; % eq. (44)
	% d2dH2 = -norm(G)/B2^2*outer(B,cross(G,B));
	% d2dH2 = d2dH2 + d2dH2';
	% d2d = [d2dF2; d2dG2; d2dH2; d2dFG; d2dGH];
	% hess = zeros(12,12);
	% for i=1:4
		% for j=1:4
			% wt = dFGHdr(:,i)*dFGHdr(:,j)';
			% wt = [diag(wt); diag(wt,-1) + diag(wt,1)];
			% hess((i-1)*3+(1:3),(j-1)*3+(1:3)) = blockProduct(wt, d2d);
		% end
	% end
	% dE2dr2 = dphi2 + hess; % eq. (28)
	idx = [3*(torsional(m,:)-1)+1; 3*(torsional(m,:)-1)+2; 3*(torsional(m,:)-1)+3];
	hessian(idx,idx) = hessian(idx,idx) + (2*K_phi)*dphi2; 
end
%hessian = 2*kphi*dE2dr2;	
     
% function [prod] = outer(x, y)
% prod = x(:)*y(:)';