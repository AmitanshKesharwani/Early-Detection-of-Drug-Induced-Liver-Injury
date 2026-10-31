"""Environment check script for DILI project.

Verifies:
1. Python version and installed package versions from requirements.txt.
2. RDKit installation and functionality by computing Morgan fingerprint,
   molecular weight, and logP for aspirin.
3. Saves all diagnostic output to reports/env_check.txt.
"""

from pathlib import Path
import sys
import importlib.metadata
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdFingerprintGenerator


def run_env_check() -> str:
    lines = []
    lines.append("=== Python & Environment Check ===")
    lines.append(f"Python Version: {sys.version}")
    lines.append("")
    lines.append("=== Installed Package Versions ===")

    requirements_file = Path(__file__).resolve().parent.parent / "requirements.txt"
    if requirements_file.exists():
        with open(requirements_file, "r", encoding="utf-8") as f:
            packages = [line.strip() for line in f if line.strip() and not line.startswith("#")]
    else:
        packages = [
            "pandas",
            "numpy",
            "scikit-learn",
            "rdkit",
            "lightgbm",
            "shap",
            "requests",
            "matplotlib",
            "pytest",
            "jupyter",
        ]

    for pkg in packages:
        try:
            ver = importlib.metadata.version(pkg)
            lines.append(f"{pkg}: {ver}")
        except importlib.metadata.PackageNotFoundError:
            lines.append(f"{pkg}: NOT INSTALLED")

    lines.append("")
    lines.append("=== RDKit Aspirin Smoke Test ===")
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    mol = Chem.MolFromSmiles(aspirin_smiles)
    if mol is None:
        raise ValueError(f"Failed to parse SMILES: {aspirin_smiles}")

    # Generate Morgan fingerprint (radius 2, 2048 bits)
    mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = mfgen.GetFingerprint(mol)
    on_bits = fp.GetNumOnBits()

    # Molecular properties
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)

    lines.append(f"Aspirin SMILES: {aspirin_smiles}")
    lines.append(f"Morgan fingerprint bits on: {on_bits} (out of 2048)")
    lines.append(f"Molecular weight: {mw:.4f}")
    lines.append(f"LogP: {logp:.4f}")

    output_text = "\n".join(lines)
    return output_text


def main():
    output_text = run_env_check()
    print(output_text)

    # Save to reports/env_check.txt
    reports_dir = Path(__file__).resolve().parent.parent / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    report_file = reports_dir / "env_check.txt"
    with open(report_file, "w", encoding="utf-8") as f:
        f.write(output_text + "\n")
    print(f"\nReport written to: {report_file}")


if __name__ == "__main__":
    main()

