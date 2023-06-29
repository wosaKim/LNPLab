def plogp(smi):
    return None

def qed(smi):
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    return Chem.QED.qed(mol)


def lo5(smi):
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    mol.UpdatePropertyCache(strict=False)
    mol = Chem.AddHs(mol, addCoords=True)

    mw = Descriptors.ExactMolWt(mol)
    hba = Descriptors.NOCount(mol)
    hbd = Descriptors.NHOHCount(mol)
    logp = Descriptors.MolLogP(mol)

    conditions = [mw <= 500, hba <= 10, hbd <= 5, logp <= 5]
    return conditions.count(True)


if __name__ == "__main__":
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    print(lo5('CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5'))
    print(qed('CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5'))


