To use sbNMA to compute Hessian matrix and normal modes, take the following steps:

1. Download the pdb file (e.g., 1ubq.pdb)

2. Use VMD automatic PSF builder to generate pdb/psf files (e.g., 1ubq_autopsf.pdb, 1ubq_autopsf.psf)

4. run charmm36prm.m in MATLAB to load the charmm parameters
>> charmm36prm

4. run charmmPDBread.m to extract xyz coordinates from the generated pdb file (e.g., 1ubq_autopsf.pdb). 
>> [ubq] = charmmPDBread('1ubq_autopsf.pdb');

5. run psfread.m to extreact bonds, angles, torsional, improper interaction terms as well as the charmm atom type of each atom in the protein from the generated psf file (e.g. 1ubq_autopsf.psf). 
>> [ubqAtomTypes, ubqBonds, ubqAngles, ubqPhis, ubqImps] = psfread('1ubq_autopsf.psf', charmmAtomTypes);

(Note: each atom in the protein is listed in the psf file under the !NATOM section and is given a charmm atom type (the 6th column, such as NH3 for the N-terminal Nitrogen; see 1ubq_autopsf.psf). A comprehensive list of all charmm atom types based on CHARMM parameter files (par_all36_prot.prm and par_all36_carb.prm) is given as a String array called charmmAtomTypes inside charmm36prm.m. There are 110 charmm atom types in the list and each charmm atom type is given a 0-based index based on its position in the list. For example, the first charmm atom type is "H" and its index is 0; the last charmm atom type in the list is "NC311" and its index is 109. By using the indices, any charmm atom type can be referred to by a number.)

6. Lastly, to generate the Hessian and mass matrices, run in MATLAB:
>> [hess, massMat] = sbNMA_PSF(ubq(:,3:5), ubqBonds, ubqAngles, ubqPhis, ubqImps, ubqAtomTypes, charmmBonds, charmmAngles, charmmPhis, charmmImps, charmmNb, charmmMass);

7. (optional) The following sequence of commands compute the mass-weighted Hessian matrix, solve for normal modes and the vibrational frequencies, and plot the vibrational frequency spectrum.  
>> massSqrt = 1./sqrt(massMat(:));
>> hessM = massSqrt*massSqrt'.*hess;
>> [Vm, Dm] = eig(hessM);
>> lambdas  = abs(diag(Dm)); % eigenvalues
>> omegas = sqrt(lambdas)*108.52; % vibrational frequency, in cm-1.
>> hist(omegas, 100);
(Note: sample run results are included in sampleRunResults.mat.)

Lastly, other matlab scripts in the directory but not mentioned above are needed directly or indirectly from the main script sbNMA_PSF.m. 


