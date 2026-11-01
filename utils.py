import requests
from rdkit import Chem
from rdkit.Chem import AllChem

# Comprehensive Molecular & Clinical Knowledge Database
MOLECULES_DB = {
    "Aspirin": {
        "smiles": "CC(=O)Oc1ccccc1C(=O)O",
        "formula": "C9H8O4",
        "iupac": "2-acetoxybenzoic acid",
        "mw": 180.16,
        "logp": 1.19,
        "h_donors": 1,
        "h_acceptors": 4,
        "tpsa": 63.6,
        "category": "NSAID / Antiplatelet",
        "uses": "Mild-to-moderate pain, fever reduction, inflammation, and secondary cardiovascular prevention (low-dose antiplatelet).",
        "mechanism": "Irreversibly inhibits cyclooxygenase (COX-1 and COX-2) enzymes, blocking prostaglandin and thromboxane synthesis.",
        "pros": ["Inexpensive & widely accessible", "Proven cardiovascular protection", "Effective antipyretic & anti-inflammatory"],
        "cons": ["Gastric mucosal irritation / ulcer risk", "Reye's syndrome risk in pediatric viral infections", "Increased bleeding propensity"],
        "dosages": {
            "Adults (Pain/Fever):": "325 mg - 650 mg every 4 hours orally (Max: 4g/day)",
            "Adults (Cardiovascular):": "75 mg - 100 mg once daily orally",
            "Pediatrics:": "Contraindicated under 16 years due to Reye's syndrome risk"
        },
        "global_context": "Widely utilized in cardiovascular disease management globally as an affordable primary/secondary prophylactic agent against myocardial infarction."
    },
    "Paracetamol": {
        "smiles": "CC(=O)Nc1ccc(O)cc1",
        "formula": "C8H9NO2",
        "iupac": "N-(4-hydroxyphenyl)acetamide",
        "mw": 151.16,
        "logp": 0.46,
        "h_donors": 1,
        "h_acceptors": 2,
        "tpsa": 29.3,
        "category": "Analgesic / Antipyretic",
        "uses": "Management of mild to moderate pain (headaches, musculoskeletal pain) and reduction of pyrexia (fever).",
        "mechanism": "Acts primarily in the central nervous system by inhibiting prostaglandin synthesis pathways and modulating cannabinoid receptors.",
        "pros": ["Gentle on gastric lining (minimal ulceration)", "Safe during pregnancy/lactation at standard doses", "Compatible with broad patient demographics"],
        "cons": ["Severe hepatotoxicity risk in overdosage", "Negligible peripheral anti-inflammatory action", "Frequent hidden ingredient in combination cold formulas"],
        "dosages": {
            "Adults:": "500 mg - 1000 mg every 4 to 6 hours orally (Max: 4g/day)",
            "Pediatrics (6-12 yrs):": "15 mg/kg body weight every 6 hours",
            "Neonates/Infants:": "Weight-based pediatric suspension per clinician instruction"
        },
        "global_context": "Essential WHO model medicine for primary healthcare pain management systems worldwide."
    },
    "Ibuprofen": {
        "smiles": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
        "formula": "C13H18O2",
        "iupac": "2-(4-isobutylphenyl)propanoic acid",
        "mw": 206.28,
        "logp": 3.5,
        "h_donors": 1,
        "h_acceptors": 2,
        "tpsa": 37.3,
        "category": "Propionic Acid Derivative NSAID",
        "uses": "Treatment of inflammatory conditions, rheumatoid arthritis, dysmenorrhea, acute mild-to-moderate pain, and fever.",
        "mechanism": "Reversibly inhibits both COX-1 and COX-2 enzymes, lowering systemic tissue concentrations of pro-inflammatory prostaglandins.",
        "pros": ["Stronger anti-inflammatory potency than paracetamol", "Rapid onset of action", "Available in diverse pediatric formulations"],
        "cons": ["Nephrotoxicity risks with chronic long-term usage", "Gastrointestinal bleeding and dyspepsia", "Exacerbation of asthma symptoms in sensitive patients"],
        "dosages": {
            "Adults:": "200 mg - 400 mg every 4 to 6 hours orally (Max: 1200 mg/day OTC)",
            "Pediatrics (>6 months):": "10 mg/kg every 6-8 hours as required"
        },
        "global_context": "Pivotal agent in musculoskeletal injury management and multi-national acute trauma care."
    },
    "Caffeine": {
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "formula": "C8H10N4O2",
        "iupac": "1,3,7-trimethylpurine-2,6-dione",
        "mw": 194.19,
        "logp": -0.07,
        "h_donors": 0,
        "h_acceptors": 3,
        "tpsa": 58.4,
        "category": "Central Nervous System Stimulant",
        "uses": "Enhancing alertness, fatigue mitigation, adjuvant in headache medication formulations, neonatal apnea treatment.",
        "mechanism": "Competitively antagonizes adenosine receptors in the central nervous system, preventing sleep-promoting adenosine signaling.",
        "pros": ["Increases cognitive vigilance and reaction speed", "Synergistic enhancement with analgesics", "Metabolic rate stimulant"],
        "cons": ["Insomnia and nervousness with excessive intake", "Physiological tolerance and mild withdrawal dependency", "Tachycardia and blood pressure elevation in sensitive hosts"],
        "dosages": {
            "Adults (Moderate):": "100 mg - 200 mg per serving (Max: 400 mg/day)",
            "Pediatrics:": "Not recommended for daily supplemental routines"
        },
        "global_context": "The world's most widely consumed psychoactive substance, embedded in economic food trade networks."
    },
    "Glucose": {
        "smiles": "C([C@@H]1[C@H]([C@@H]([C@H]([C@H](O1)O)O)O)O)O",
        "formula": "C6H12O6",
        "iupac": "(2R,3R,4S,5S,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol",
        "mw": 180.16,
        "logp": -3.24,
        "h_donors": 5,
        "h_acceptors": 6,
        "tpsa": 110.38,
        "category": "Monosaccharide / Energy Substrate",
        "uses": "Primary cellular metabolic fuel, intravenous management of severe hypoglycemia, oral rehydration solutions.",
        "mechanism": "Transported into cells via insulin-dependent GLUT pathways to enter glycolysis and the Krebs cycle for ATP generation.",
        "pros": ["Immediate biological energy source", "Crucial for cerebral and myocardial metabolic function", "Life-saving intervention in acute hypoglycemia"],
        "cons": ["Hyperglycemia and metabolic complications in uncontrolled diabetes", "Obesity risks from chronic excess simple sugar intake"],
        "dosages": {
            "Hypoglycemia Rescue:": "15 g - 20 g oral fast-acting glucose tablets",
            "Intravenous Infusion:": "25 mL to 50 mL of 50% Dextrose solution under clinical supervision"
        },
        "global_context": "Fundamental biochemical pillar underpinning cellular respiration across all aerobic life domains."
    },
    "Benzene": {
        "smiles": "c1ccccc1",
        "formula": "C6H6",
        "iupac": "benzene",
        "mw": 78.11,
        "logp": 2.13,
        "h_donors": 0,
        "h_acceptors": 0,
        "tpsa": 0.0,
        "category": "Aromatic Hydrocarbon (Industrial Chemical)",
        "uses": "Chemical feedstock for manufacturing plastics, synthetic resins, lubricants, dyes, rubbers, and pharmaceutical intermediates.",
        "mechanism": "Volatile hydrocarbon absorbed via inhalation or dermal pathways, metabolizing in the liver via cytochrome P450 into reactive epoxides.",
        "pros": ["Indispensable industrial precursor for chemical synthesis", "High octane blending component historically in petrochemical streams"],
        "cons": ["Confirmed Group 1 Human Carcinogen (Leukemia risk)", "Bone marrow depression and immunotoxicity", "Acute central nervous system depression at high vapor concentrations"],
        "dosages": {
            "Therapeutic Use:": "None (Strictly toxic industrial chemical)"
        },
        "global_context": "Subject to stringent international workplace safety regulations and strict environmental emission controls."
    }
}

GLOBAL_DISEASES_DB = {
    "Cardiovascular Diseases (CVDs)": {
        "global_burden": "Leading global cause of mortality, responsible for an estimated 17.9 million deaths annually worldwide.",
        "key_pathology": "Atherosclerosis, plaque accumulation, arterial stiffness leading to myocardial infarction and stroke.",
        "standard_medicines": ["Aspirin", "Statins (e.g., Atorvastatin)", "ACE Inhibitors (e.g., Lisinopril)"],
        "chemical_alternatives": "Omega-3 polyunsaturated fatty acids, natural polyphenol antioxidants (resveratrol, curcumin) for endothelial support.",
        "preventative_strategy": "Blood pressure management, strict glycemic control, balanced nutrition, and structured aerobic activity."
    },
    "Diabetes Mellitus (Type 2)": {
        "global_burden": "Over 500 million adults affected globally, presenting chronic metabolic strain on healthcare systems.",
        "key_pathology": "Peripheral insulin resistance combined with progressive pancreatic beta-cell dysfunction.",
        "standard_medicines": ["Metformin", "Insulin formulations", "SGLT2 Inhibitors (e.g., Empagliflozin)"],
        "chemical_alternatives": "Berberine (natural isoquinoline alkaloid targeting AMPK pathways), alpha-lipoic acid for neuropathy support.",
        "preventative_strategy": "Dietary refinement emphasizing low glycemic index intake, weight management, and physical conditioning."
    },
    "Infectious Respiratory Pandemics & Pathogens": {
        "global_burden": "Recurrent global health crises driven by mutating viral strains (Influenza, Coronaviruses).",
        "key_pathology": "Severe acute respiratory distress syndrome (ARDS), systemic cytokine storm, pulmonary epithelial damage.",
        "standard_medicines": ["Antivirals (Paxlovid/Remdesivir)", "Corticosteroids (Dexamethasone)", "Bronchodilators"],
        "chemical_alternatives": "Immuno-modulatory botanical extracts (Echinacea, Elderberry), zinc ionophores (Quercetin).",
        "preventative_strategy": "Global vaccination tracking, broad-spectrum antiviral screening, and advanced epidemiological surveillance."
    },
    "Antimicrobial Resistance (AMR)": {
        "global_burden": "Emerging silent pandemic threatening modern surgical procedures and standard infection controls.",
        "key_pathology": "Bacterial mutation and horizontal gene transfer rendering traditional beta-lactams and antibiotics ineffective.",
        "standard_medicines": ["Reserved Carbapenems", "Novel beta-lactamase inhibitor combinations"],
        "chemical_alternatives": "Bacteriophage therapy, antimicrobial peptide engineering, silver nanoparticle formulations.",
        "preventative_strategy": "Rational antibiotic stewardship, strict hygiene protocols, and rapid diagnostic biomarker assays."
    }
}

def get_molecule_data(query):
    """Fetches molecule data from local dictionary or queries PubChem PUG REST API."""
    q_clean = query.strip()
    # Check local db first
    for name, data in MOLECULES_DB.items():
        if q_clean.lower() == name.lower() or q_clean.upper() == data["smiles"].upper():
            return {
                "name": name,
                "smiles": data["smiles"],
                "formula": data["formula"],
                "iupac": data["iupac"],
                "mw": data["mw"],
                "logp": data["logp"],
                "h_donors": data["h_donors"],
                "h_acceptors": data["h_acceptors"],
                "tpsa": data["tpsa"],
                "category": data["category"],
                "uses": data["uses"],
                "mechanism": data["mechanism"],
                "pros": data["pros"],
                "cons": data["cons"],
                "dosages": data["dosages"],
                "global_context": data["global_context"],
                "source": "SRM Verified Curated Database"
            }
            
    # Fallback to PubChem API for custom queries
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q_clean}/property/MolecularFormula,MolecularWeight,IUPACName,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA,CanonicalSMILES/JSON"
        res = requests.get(url, timeout=5)
        if res.status_code == 200:
            props = res.json()["PropertyTable"]["Properties"][0]
            return {
                "name": q_clean.capitalize(),
                "smiles": props.get("CanonicalSMILES", q_clean),
                "formula": props.get("MolecularFormula", "N/A"),
                "iupac": props.get("IUPACName", "N/A"),
                "mw": props.get("MolecularWeight", 0.0),
                "logp": props.get("XLogP", 0.0),
                "h_donors": props.get("HBondDonorCount", 0),
                "h_acceptors": props.get("HBondAcceptorCount", 0),
                "tpsa": props.get("TPSA", 0.0),
                "category": "Custom Queried Compound",
                "uses": "General chemical entity queried from public repositories.",
                "mechanism": "Standard molecular conformation analyzed via RDKit engine.",
                "pros": ["Dynamic query resolution", "Standardized cheminformatics properties"],
                "cons": ["Requires clinical trial validation for human therapeutics"],
                "dosages": {"Standard Guidance:": "Consult certified pharmacological literature."},
                "global_context": "Evaluated through international open-access chemical indexing.",
                "source": "PubChem Rest API"
            }
    except Exception:
        pass

    # RDKit direct SMILES parser fallback
    mol = Chem.MolFromSmiles(q_clean)
    if mol:
        formula = Chem.MolToSmiles(mol)
        return {
            "name": "Custom SMILES Structure",
            "smiles": q_clean,
            "formula": Chem.CalcMolFormula(mol),
            "iupac": "Computed IUPAC",
            "mw": round(Chem.Descriptors.MolWt(mol), 2),
            "logp": round(Chem.Crippen.MolLogP(mol), 2),
            "h_donors": Chem.Lipinski.NumHDonors(mol),
            "h_acceptors": Chem.Lipinski.NumHAcceptors(mol),
            "tpsa": round(Chem.RdMolDescriptors.CalcTPSA(mol), 2),
            "category": "User Synthesized SMILES",
            "uses": "Custom structure parsed directly via RDKit backend.",
            "mechanism": "Evaluated via structural bond topology.",
            "pros": ["Instant computational parsing"],
            "cons": ["Requires experimental verification"],
            "dosages": {"Guidance:": "Experimental compound - No clinical dosage assigned."},
            "global_context": "Research scale chemical component.",
            "source": "RDKit Engine"
        }
    return None

def generate_3d_sdf(smiles):
    """Generates 3D coordinates and returns an SDF block using RDKit."""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        # Embed 3D coordinates
        id = AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        if id < 0:
            return None
        AllChem.UFFOptimizeMolecule(mol)
        block = Chem.MolToMolBlock(mol)
        return block
    except Exception:
        return None
