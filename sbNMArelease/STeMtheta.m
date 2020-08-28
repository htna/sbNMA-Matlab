function hessian=STeMtheta(x,angles, ktheta)
% a new implmentation of sbNMA bond angle term. 
% Guang Song
% derive the hessian of the second term
% reference: Blondel and Karplus, NEW FORMULATION FOR ELIMINATION OF StNGULARITIES, Journal of Computational Chemistry, Vol. 17, No. 9, 1132-1141 (1996)
% x: xyz coordinates of protein
% angles: a n-by-3 matrix with each row of 3 ints representing 0-based indices of the three atoms that have bond angle interactions. 
% ktheta: the spring constant of each interaction.
% hessian:  hessian matrix due to bond angle interactions.
%
%
% Note: all eq. numbers are from the reference. 
% Note: at present, we are NOT consider d2theta/dr2 term but only outer product of (dtheta/dr) with itself.
%
% Copyright © 2020 Guang Song. All Rights Reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

if length(ktheta)==1
   ktheta = ktheta(ones(size(angles,1),1)); % duplicate it for each angle 09/27/17
end
hessian=zeros(size(x,1)*3);

for m=1:size(angles,1)
	K_theta = ktheta(m); % added on 09/27/17 for sbNMA
    r=x(angles(m,:),1:3);
    A = r(1,:)-r(2,:);
	B = r(3,:)-r(2,:);
	G = cross(A,B);
	G = G/norm(G);
	dphidAB = [cross(G,A)/sum(A.^2); cross(B,G)/sum(B.^2);]'; % eqs. (16) and (17)
	dABdr = [1 -1 0; 
			 0 -1 1;]; % eq. (26)
	dphiDr = dphidAB*dABdr;
	dphi2 = dphiDr(:)*dphiDr(:)';
	
	idx = [3*(angles(m,:)-1)+1; 3*(angles(m,:)-1)+2; 3*(angles(m,:)-1)+3];
	hessian(idx,idx) = hessian(idx,idx) + (2*K_theta)*dphi2; 
end
%hessian = 2*kphi*dE2dr2;	
     
% function [prod] = outer(x, y)
% prod = x(:)*y(:)';