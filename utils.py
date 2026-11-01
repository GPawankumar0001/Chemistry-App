"""Helper functions for SMILES in Chemistry (SRM IST, Tiruchirappalli)."""
import requests
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Draw, rdMolDescriptors

PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
TIMEOUT = 12

PRESETS = {
    "Aspirin": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Paracetamol": "CC(=O)NC1=CC=C(O)C=C1",
    "Benzene": "C1=CC=CC=C1",
    "Glucose": "C(C1C(C(C(C(O1)O)O)O)O)O",
    "Ibuprofen": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
}

INSIGHTS = {
    "Aspirin": {
        "use": "Pain, fever, inflammation; low-dose antiplatelet therapy to prevent heart attack and stroke.",
        "moa": "Irreversibly acetylates COX-1 and COX-2, blocking prostaglandin and thromboxane synthesis.",
        "pros": ["Inexpensive and widely available", "Proven cardiovascular prevention", "Anti-inflammatory and antipyretic"],
        "cons": ["Gastric irritation and bleeding risk", "Reye's syndrome risk in children", "Not suitable with some anticoagulants"],
    },
    "Caffeine": {
        "use": "CNS stimulant for alertness; used in apnea of prematurity and in some analgesic combinations.",
        "moa": "Antagonises adenosine A1/A2A receptors and inhibits phosphodiesterases.",
        "pros": ["Improves alertness and reaction time", "Rapidly absorbed", "Well studied at normal intake"],
        "cons": ["Anxiety, tremor, insomnia", "Dependence and withdrawal headache", "Raises heart rate at high doses"],
    },
    "Paracetamol": {
        "use": "Mild to moderate pain and fever.",
        "moa": "Thought to inhibit central COX pathways and modulate the endocannabinoid system; weak peripheral anti-inflammatory effect.",
        "pros": ["Gentle on the stomach", "Safe in pregnancy at normal doses", "Fast onset"],
        "cons": ["Liver toxicity in overdose", "Little anti-inflammatory effect", "Hidden in many combination products"],
    },
    "Benzene": {
        "use": "Industrial solvent and feedstock for styrene, phenol and nylon precursors. Not a drug.",
        "moa": "Metabolised to reactive epoxides and quinones that damage bone marrow and DNA.",
        "pros": ["Key building block in petrochemicals", "Excellent nonpolar solvent", "High octane contribution"],
        "cons": ["Known human carcinogen", "Flammable and volatile", "Strict exposure limits; largely replaced by safer solvents"],
    },
    "Glucose": {
        "use": "Primary energy substrate; IV dextrose for hypoglycaemia and fluid therapy.",
        "moa": "Enters glycolysis and oxidative phosphorylation to produce ATP.",
        "pros": ["Immediate energy source", "Very safe and well tolerated", "Cheap and stable"],
        "cons": ["Hyperglycaemia in diabetes", "Contributes to dental decay", "Excess intake promotes weight gain"],
    },
    "Ibuprofen": {
        "use": "Pain, fever and inflammation, including arthritis and dysmenorrhoea.",
        "moa": "Reversibly inhibits COX-1 and COX-2, reducing prostaglandin synthesis.",
        "pros": ["Effective anti-inflammatory", "Reversible COX inhibition", "Available over the counter"],
        "cons": ["Gastric ulcer and kidney risk", "Raises blood pressure with long use", "Avoid in late pregnancy"],
    },
}


@st.cache_data(show_spinner=False, ttl=3600)
def _pubchem(kind: str, value: str):
    """Return a PubChem property dict for a compound name or SMILES, or None."""
    props = "SMILES,IUPACName,MolecularFormula,MolecularWeight,Title"
    try:
        if kind == "name":
            r = requests.get(f"{PUG}/compound/name/{requests.utils.quote(value)}/property/{props}/JSON", timeout=TIMEOUT)
        else:
            r = requests.post(f"{PUG}/compound/smiles/property/{props}/JSON", data={"smiles": value}, timeout=TIMEOUT)
        if r.status_code != 200:
            return None
        return r.json()["PropertyTable"]["Properties"][0]
    except (requests.RequestException, KeyError, IndexError, ValueError):
        return None


def resolve(query: str):
    """Turn a name or SMILES into {'name','smiles','iupac','source'}; None if unresolved."""
    q = (query or "").strip()
    if not q:
        return None
    for name, smi in PRESETS.items():
        if q.lower() == name.lower():
            hit = _pubchem("smiles", smi)
            return {"name": name, "smiles": smi, "iupac": (hit or {}).get("IUPACName", "n/a"), "source": "Preset"}
    if " " not in q and Chem.MolFromSmiles(q) is not None:
        hit = _pubchem("smiles", q)
        return {"name": (hit or {}).get("Title", "Custom molecule"), "smiles": q,
                "iupac": (hit or {}).get("IUPACName", "n/a"), "source": "SMILES input"}
    hit = _pubchem("name", q)
    if hit:
        smi = hit.get("SMILES") or hit.get("ConnectivitySMILES") or hit.get("IsomericSMILES")
        if smi and Chem.MolFromSmiles(smi):
            return {"name": hit.get("Title", q.title()), "smiles": smi,
                    "iupac": hit.get("IUPACName", "n/a"), "source": "PubChem"}
    return None


def properties(smiles: str) -> dict:
    """Core descriptors computed with RDKit."""
    m = Chem.MolFromSmiles(smiles)
    hbd, hba = rdMolDescriptors.CalcNumHBD(m), rdMolDescriptors.CalcNumHBA(m)
    mw, logp = Descriptors.MolWt(m), Crippen.MolLogP(m)
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return {
        "Molecular Formula": rdMolDescriptors.CalcMolFormula(m),
        "Molecular Weight (g/mol)": round(mw, 2),
        "LogP": round(logp, 2),
        "H-Bond Donors": hbd,
        "H-Bond Acceptors": hba,
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(m), 2),
        "Rotatable Bonds": rdMolDescriptors.CalcNumRotatableBonds(m),
        "Rings": rdMolDescriptors.CalcNumRings(m),
        "Heavy Atoms": m.GetNumHeavyAtoms(),
        "Lipinski Violations": violations,
    }


def element_table(smiles: str):
    """Per-element atom counts with atomic number and mass."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    pt, counts = Chem.GetPeriodicTable(), {}
    for a in m.GetAtoms():
        counts[a.GetAtomicNum()] = counts.get(a.GetAtomicNum(), 0) + 1
    return [{"Symbol": pt.GetElementSymbol(z), "Element": pt.GetElementName(z) if hasattr(pt, "GetElementName") else pt.GetElementSymbol(z),
             "Atomic Number": z, "Atomic Mass": round(pt.GetAtomicWeight(z), 3), "Atoms in Molecule": n}
            for z, n in sorted(counts.items())]


def mol_image(smiles: str, size=(420, 320)):
    return Draw.MolToImage(Chem.MolFromSmiles(smiles), size=size)


def sdf_block(smiles: str):
    """3D coordinates (ETKDG + MMFF) as an SDF block for py3Dmol; None on failure."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(m, randomSeed=42) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception:
        pass
    return Chem.MolToMolBlock(m)


def viewer_html(sdf: str, style="Ball and stick", spin=False, width=700, height=480) -> str:
    import py3Dmol
    v = py3Dmol.view(width=width, height=height)
    v.addModel(sdf, "sdf")
    styles = {
        "Ball and stick": {"stick": {"radius": 0.15}, "sphere": {"scale": 0.28}},
        "Stick": {"stick": {"radius": 0.2}},
        "Space filling": {"sphere": {}},
        "Wireframe": {"line": {}},
    }
    v.setStyle(styles[style])
    v.setBackgroundColor("#0B2545")
    v.zoomTo()
    if spin:
        v.spin(True)
    return v._make_html()
