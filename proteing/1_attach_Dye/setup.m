%% read pdb of ww-domain
[pdb_pro, crd_pro] = readpdb('K10C_Alexa488.pdb');

%% attachh Alexa 594 at C-terminal
index = selectid(pdb_pro.resseq, 56) & selectname(pdb_pro.name, 'C');
crd = crd_pro(to3(index));
crd_pro(1:3:end) = crd_pro(1:3:end) - crd(1);
crd_pro(2:3:end) = crd_pro(2:3:end) - crd(2);
crd_pro(3:3:end) = crd_pro(3:3:end) - crd(3);

[pdb_dye, crd_dye] = readpdb('A59_C1C.pdb');
index = selectname(pdb_dye.resname, 'C1C') & selectname(pdb_dye.name, 'N');
crd = crd_dye(to3(index));
crd_dye(1:3:end) = crd_dye(1:3:end) - crd(1);
crd_dye(2:3:end) = crd_dye(2:3:end) - crd(2);
crd_dye(3:3:end) = crd_dye(3:3:end) - crd(3);

crd_dye(1:3:end) = crd_dye(1:3:end) - 1.335;

[pair, dist] = searchrange(crd_pro, crd_dye, 8.0);
dist_min = min(dist);
dist_mean = mean(dist);

for i = 1:1000
  [R, ~] = qr(randn(3));
  crd = reshape(crd_dye, 3, []);
  crd = R*crd;
  crd = crd(:)';
  [pair, dist] = searchrange(crd_pro, crd, 8.0);

  if (min(dist) > dist_min) & (mean(dist) > dist_mean)
    crd_dye = crd;
    dist_min = min(dist);
    dist_mean = mean(dist);
  end
end

%% cleaning of residue id, etc.
pdb_dye.xyz = reshape(crd_dye, 3, [])';
pdb_pro.xyz = reshape(crd_pro, 3, [])';
pdb_pro = addstruct(pdb_pro, pdb_dye);

%pdb_pro.resseq = pdb_pro.resseq + 3;

% index = selectname(pdb_pro.resname, 'A48');
% pdb_pro.resseq(index) = 1;n
% index = selectname(pdb_pro.resname, 'C1R');
% pdb_pro.resseq(index) = 2;
% index = selectname(pdb_pro.resname, 'A64');
% pdb_pro.resseq(index) = 42;
% index = selectname(pdb_pro.resname, 'C1R');
% pdb_pro.resseq(index) = 43;

natom = size(pdb_pro.name, 1);
pdb_pro.serial = (1:natom)';

for iatom = 1:natom
  pdb_pro.chainid(iatom) = 'A';
end

%% write pdb and save the workspace
writepdb('setup.pdb', pdb_pro);
save -v7.3 setup.mat;

%% delete OXT manually
%% renumber residue index manually
