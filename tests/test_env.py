"""Test environment and core cheminformatics functionality."""

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


def test_aspirin_fingerprint_bits():
    """Verify that Morgan fingerprint generation for aspirin produces valid on-bits."""
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = Chem.MolFromSmiles(aspirin_smiles)
    assert mol is not None, "Failed to parse aspirin SMILES string."

    mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = mfgen.GetFingerprint(mol)
    on_bits = fp.GetNumOnBits()

    assert 0 < on_bits < 2048, f"Expected 0 < on_bits < 2048, got {on_bits}"

