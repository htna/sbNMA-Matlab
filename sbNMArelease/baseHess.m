function anm = baseHess(x,cx)
% Basic Hessian Matrix
% Author: Guang Song
% Created: Feb 23, 2005
% Updated: 09/28/2019
%
% x: a 3-by-N matrix that contains the coordinates of N atoms.
% cx: the contact matrix. Also with gama info 
% anm: ANM-like Hessian matrix 
%
% Copyright © 2020 Guang Song. All Rights Reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

cx = tril(cx); % contact matrix is symmetric. Only half of the matrix is needed to find all the contacts. 
dim = size(cx,1);
[I,J, K] = find(cx); % pairs of contacts
%nct = length(I); % # of contacts

dR  = x(J,:) - x(I,:);
nrm = normMat(dR,2);
dR = dR./repmat(nrm, 1, 3);
dRa = [repmat(dR(:,1),1,3), repmat(dR(:,2),1,3), repmat(dR(:,3),1,3)];
dRb = repmat(dR, 1, 3);
dR2 = -dRa.*dRb;
Is = [repmat(3*I-2, 1, 3), repmat(3*I-1, 1, 3), repmat(3*I, 1, 3)];
Js = repmat([3*J-2, 3*J-1, 3*J], 1, 3);
Ks = repmat(K, 1, 9);
anm = sparse(Is(:), Js(:), Ks(:).*dR2(:), 3*dim, 3*dim);
anm = anm + anm';
anm = anm - blockSumDiag(anm);
anm = full(anm); % change it from a sparse matrix to a full matrix