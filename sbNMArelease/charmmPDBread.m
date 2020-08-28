function [xyz]=charmmPDBread(charmmPDBfilename)
% Read a charmm pdb file and saved the results into the following array:
% xyz: a N by 7 arrays, where N is the number of atoms, whose columns are: 
% column 1: atom sequence number (1 to N)
% column 2: residue number
% column 3-5: coordinates
% column 6 (char): atom types, e.g., 'C', 'N', 'O', 'H', 'S'.  
% column 7 (char): backbone atoms, 'A' for C_alpha, 'N' for N, 'C' for C, 'O' for O, 0 for everything else. 

% Copyright © 2020 Guang Song. All Rights Reserved.
% 
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
fid=fopen(charmmPDBfilename,'r');
A = fread(fid);
lines = splitlines(string(char(A')));
lines = lines(find(startsWith(lines, "ATOM")),:); % key only ATOM rows
lines = char(lines); % convert a matrix of chars
xyz = zeros(size(lines,1), 7);
xyz(:,1:5) = [str2num(lines(:,7:11)), str2num(lines(:,23:26)), str2num(lines(:,31:38)), str2num(lines(:,39:46)), str2num(lines(:,47:54))];
atomType = strip(string(lines(:,13:16)));
%ATOM      1  N   GLU A   1      25.081 -10.699 -18.389  1.00  0.00      AP1  N
atomtypes = ['C', 'N', 'O', 'H', 'S'];
for i=1:length(atomtypes)
	xyz(find(startsWith(atomType, atomtypes(i))),6)=atomtypes(i);
end
bbatoms = ["CA", "C", "N", "O", "OT1"]; % backbone atoms. OT1: C-terminal O
bbatomsNames = ['A', 'C', 'N', 'O', 'O'];
for i=1:length(bbatoms)
	xyz(find(matches(atomType, bbatoms(i))),7)=bbatomsNames(i);
end