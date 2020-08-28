function [atomTypes, bonds, angles, phis, imps, charges]=psfread(filename, charmmAtoms)
% Read a protein structure file and saved the results into the following arrays:
% atomTypes: 
% bonds:
% angles: bond angles
% phis: dihedral angles
% imps: inpromper angles
%
% Copyright (c) 2020 Guang Song. All Rights Reserved. 
%
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

fid=fopen(filename,'r');
A = fread(fid);
lines = splitlines(string(char(A')));
atomLine = find(endsWith(lines, "!NATOM"));
temp = split(strip(lines(atomLine)));
natom = str2num(temp(1));
atomRows = lines(atomLine+(1:natom));
atomRows = split(strip(atomRows));
% The columns of atomRows in a psf file looks like:
%1 AP1  1    GLU  N    NH3   -0.300000       14.0070           0
charges = str2double(atomRows(:,7)); % the 7th column is charge
%atomTypes = atomRows(:, 6);
atomTypes = zeros(size(atomRows,1),1);
for i=1:size(atomRows,1)
	atomTypes(i) = find(matches(charmmAtoms, atomRows(i, 6)))-1; % column 6 is charmm atom type of atom i; =1: make it 0-based
end


blocks = ["!NBOND: bonds", "!NTHETA: angles", "!NPHI: dihedrals", "!NIMPHI: impropers", "!NDON: donors"];
for i=1:length(blocks)
  row(i) = find(endsWith(lines, blocks(i)));
end
% retrieve the data in each block
for i=1:length(blocks)-1
  data{i} = str2double(split(strip(join(lines(row(i)+1:row(i+1)-1)))));
end
bonds = vec2mat(data{1}, 2);
angles = vec2mat(data{2}, 3);
phis = vec2mat(data{3}, 4);
imps = vec2mat(data{4}, 4);


