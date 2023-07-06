import os
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


class Druglikeness:
    def __init__(self, smi):
        self.smi = smi
        self.mol = self.plain_mol()
        self.mol_h = self.h_mol()

    def plain_mol(self):
        return Chem.MolFromSmiles(self.smi, sanitize=True)

    def h_mol(self):
        result = self.plain_mol()
        result.UpdatePropertyCache(strict=False)
        return Chem.AddHs(result, addCoords=True)

    def sascore(self):
        return sascorer.calculateScore(self.mol)

    def plogp(self):
        return None

    def qed(self):
        return Chem.QED.qed(self.mol)

    def lo5(self):
        mw = Descriptors.ExactMolWt(self.mol_h)
        hba = Descriptors.NOCount(self.mol_h)
        hbd = Descriptors.NHOHCount(self.mol_h)
        logp = Descriptors.MolLogP(self.mol_h)
        conditions = [mw <= 500, hba <= 10, hbd <= 5, logp <= 5]
        return conditions.count(True)


if __name__ == "__main__":
    dln = Druglikeness('CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5')
    print(dln.qed(), dln.lo5(), dln.sascore())


