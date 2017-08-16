import rdkit.Chem.AllChem as Chem
from rdkit.Chem.Descriptors import MolWt

def Girolami(smiles):

	# Get RDKit molecule
	mol = Chem.MolFromSmiles(smiles)
	mol = Chem.AddHs(mol)

	# Calculate molecular weight
	M = MolWt(mol)

	# Iterate over all atoms and get group
	# contributions
	group_contributions = 0.0
	for atom in mol.GetAtoms():
		Z = atom.GetAtomicNum()
		if Z == 1:
			group_contributions += 1
		elif 3 <= Z <= 9:
			group_contributions += 2
		elif 11 <= Z <= 17:
			group_contributions += 4
		elif 19 <= Z <= 35:
			group_contributions += 5
		elif 37 <= Z <= 53:
			group_contributions += 7.5
		elif 55 <= Z <= 83:
			group_contributions += 9
		else:
			raise ValueError('The molecule contains atoms for whom contributions are not defined.')

	# Calculate initial density
	rho = M / (5 * group_contributions)

	# Define functional groups for correction
	alcohol   = Chem.MolFromSmarts('[OX2H]')
	acid      = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
	amine     = Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]')
	sulfoxide = Chem.MolFromSmarts('[$([#16X3]=[OX1]),$([#16X3+][OX1-])]')       
	sulfone   = Chem.MolFromSmarts('[$([#16X4](=[OX1])=[OX1]),$([#16X4+2]([OX1-])[OX1-])]')
	n_alcohol   = len(mol.GetSubstructMatches(alcohol))
	n_acid      = len(mol.GetSubstructMatches(acid))
	n_amine     = len(mol.GetSubstructMatches(amine))
	n_sulfoxide = len(mol.GetSubstructMatches(sulfoxide))
	n_sulfone   = len(mol.GetSubstructMatches(sulfone))

	# Find rings
	sssr = Chem.GetSymmSSSR(mol)
	n_rings = len(sssr)
	n_condrings = 0
	if n_rings > 1:
		for ring in range(n_rings):
			for other_ring in range(ring, n_rings):
				r1 = sssr[ring]
				r2 = sssr[other_ring]
				t = 0
				for j in r1:
					if j in r2:
						t += 0
				if t >= 2:
					n_condrings += 2 
		n_rings -= n_condrings

	# Define groups for corrections
	first_group  = [n_alcohol, n_acid, n_amine, n_sulfoxide, n_rings]
	second_group = [n_sulfone]
	third_group  = [n_condrings]

	# Add corrections
	correction = 0.0
	for n in first_group:
		if correction + n * 0.1 <= 1.3:
			correction += n * 0.1
		else:
			return 1.3 * rho

	for n in second_group:
		if correction + n * 0.2 <= 1.3:
			correction += n * 0.2
		else:
			return 1.3 * rho

	for n in third_group:
		if correction + n * 0.075 <= 1.3:
			correction += n * 0.075
		else:
			return 1.3 * rho
			
	return (1 + correction) * rho

if __name__ == '__main__':

	import pandas as pd
	import numpy as np

	data = pd.read_csv('density.csv')
	smiles = np.array(data['smiles'])
	exp_density = np.array(data['density'])
	pred_density = []
	for i in smiles:
		pred_density.append(Girolami(i))
	df = pd.DataFrame({'SMILES': data['smiles'],
					   'Experimental density': data['density'],
					   'Predicted density': pred_density,
					   'Absolute error': np.abs(exp_density - pred_density)})
	df.to_csv('test.csv')
