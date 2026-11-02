"""
utils.py: backend engine for "SMILES in Chemistry"
SRM Institute of Science and Technology (SRM IST), Tiruchirappalli
Department of Computer Science & Chemistry Research

Handles: PubChem REST lookups, RDKit parsing, physicochemical + ADME descriptors,
drug-likeness rules, structural alerts, 3D (ETKDG + UFF) SDF generation, 2D depiction,
a curated educational molecule library, and the Global Health knowledge base.

EDUCATIONAL USE ONLY. Nothing here is medical advice.
"""
from __future__ import annotations

import math
from functools import lru_cache
from typing import Any, Dict, List, Optional

import requests
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, QED, rdMolDescriptors, rdMolDraw2D
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

RDLogger.DisableLog("rdApp.*")

# --------------------------------------------------------------------------- #
# Constants
# --------------------------------------------------------------------------- #
PUBCHEM_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_COMPOUND_PAGE = "https://pubchem.ncbi.nlm.nih.gov/compound/"
DRUGBANK_PAGE = "https://go.drugbank.com/drugs/"
DAILYMED_SEARCH = "https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query="
HTTP_HEADERS = {"User-Agent": "SMILES-in-Chemistry/1.0 (SRM IST Tiruchirappalli; educational)"}
HTTP_TIMEOUT = 12
MAX_SMILES_LEN = 600
MAX_HEAVY_ATOMS_3D = 160

PUBCHEM_PROPERTY_SETS = [
    "Title,IUPACName,MolecularFormula,MolecularWeight,ExactMass,XLogP,TPSA,SMILES,InChIKey",
    "Title,IUPACName,MolecularFormula,MolecularWeight,ExactMass,XLogP,TPSA,IsomericSMILES,InChIKey",
]


def _d(dose: str, limit: str, safety: str) -> Dict[str, str]:
    return {"dose": dose, "limit": limit, "safety": safety}


# --------------------------------------------------------------------------- #
# Curated educational library
# --------------------------------------------------------------------------- #
LIBRARY: Dict[str, Dict[str, Any]] = {
    "aspirin": {
        "name": "Aspirin (Acetylsalicylic acid)",
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "cid": 2244, "drugbank": "DB00945", "dailymed": "aspirin",
        "aliases": ["acetylsalicylic acid", "asa", "acetylsalicylate"],
        "drug_class": "NSAID; irreversible COX inhibitor; antiplatelet",
        "summary": "Aspirin irreversibly acetylates cyclooxygenase (COX-1/COX-2), reducing prostaglandin and thromboxane A2 synthesis. It is analgesic, antipyretic, anti-inflammatory and antiplatelet.",
        "indications": ["Pain and fever", "Secondary prevention of myocardial infarction and ischaemic stroke", "Inflammatory conditions (specialist use)"],
        "pros": ["Inexpensive and widely available", "Proven secondary cardiovascular prevention", "Multiple therapeutic actions"],
        "cons": ["GI irritation and bleeding risk", "Reye syndrome risk in children with viral illness", "Hypersensitivity / NSAID-exacerbated asthma"],
        "adme": {
            "Absorption": "Rapidly absorbed in stomach and upper small intestine. Substantial pre-systemic hydrolysis to salicylate. Enteric coating delays absorption.",
            "Distribution": "Salicylate distributes widely (Vd about 0.15 L/kg) and binds albumin about 80-90% (concentration-dependent). Crosses the placenta.",
            "Metabolism": "Esterases hydrolyse aspirin to salicylic acid; then conjugation to salicyluric acid and glucuronides. Pathways saturate at higher doses (non-linear kinetics).",
            "Excretion": "Renal. Aspirin t1/2 is about 15-20 min; salicylate t1/2 is 2-3 h at low dose and up to 15-30 h at high dose. Urinary alkalinisation increases elimination.",
        },
        "pk": {"dose_mg": 500, "f": 0.8, "vd_l_kg": 0.15, "t_half_h": 3.0, "tmax_h": 1.0,
               "note": "Parameters approximate the active metabolite salicylate after a low oral dose."},
        "pk_table": {"Bioavailability": "about 80-100% (as salicylate)", "Tmax": "about 1-2 h", "Half-life": "2-3 h (low dose)", "Protein binding": "80-90% (salicylate)", "Vd": "about 0.15 L/kg"},
        "tox": {
            "ld50": [{"species": "Rat", "route": "Oral", "value": 200.0, "unit": "mg/kg", "source": "HSDB via PubChem", "url": PUBCHEM_COMPOUND_PAGE + "2244#section=Toxicity"}],
            "clinical": [
                {"marker": "GI bleeding / ulceration", "detail": "Dose-related; higher with alcohol, anticoagulants, age >60.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "aspirin"},
                {"marker": "Reye syndrome", "detail": "Avoid in children/adolescents with viral illness (influenza, varicella).", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "aspirin"},
                {"marker": "Salicylism", "detail": "Tinnitus, hyperventilation, metabolic acidosis in overdose.", "source": "PubChem toxicity section", "url": PUBCHEM_COMPOUND_PAGE + "2244#section=Toxicity"},
            ],
        },
        "dosing": {
            "Adult": _d("325-650 mg every 4-6 h for pain/fever", "Max about 4 g/day (OTC labelling)", "Low-dose (75-100 mg) cardiovascular use only on clinician advice. Avoid with active ulcer, bleeding disorders, late pregnancy."),
            "Pediatric": _d("Not for routine analgesia/antipyresis in children and adolescents", "Specialist-directed only (e.g., Kawasaki disease)", "Reye syndrome risk. Do not self-administer to children or teenagers with viral illness."),
            "Geriatric": _d("Lowest effective dose, shortest duration", "Avoid high-dose chronic use", "Higher GI/intracranial bleeding risk. Guidelines advise against routinely starting low-dose aspirin for primary prevention in adults 60 and older (USPSTF 2022)."),
        },
        "conditions": ["cardiovascular"],
    },
    "paracetamol": {
        "name": "Paracetamol (Acetaminophen)",
        "smiles": "CC(=O)NC1=CC=C(O)C=C1",
        "cid": 1983, "drugbank": "DB00316", "dailymed": "acetaminophen",
        "aliases": ["acetaminophen", "tylenol", "apap", "n-acetyl-p-aminophenol"],
        "drug_class": "Analgesic / antipyretic",
        "summary": "Paracetamol is a widely used analgesic and antipyretic with weak peripheral anti-inflammatory activity. Its mechanism is incompletely understood and includes central COX modulation and descending serotonergic pathways.",
        "indications": ["Mild-to-moderate pain", "Fever"],
        "pros": ["Good GI tolerability compared with NSAIDs", "Safe at therapeutic doses in most populations", "Suitable in pregnancy at recommended doses"],
        "cons": ["Dose-dependent hepatotoxicity in overdose", "Little anti-inflammatory effect", "Narrow margin when combined with other paracetamol-containing products"],
        "adme": {
            "Absorption": "Rapid oral absorption, mainly small intestine. Oral bioavailability about 63-89%. Food delays Tmax but not total absorption.",
            "Distribution": "Vd about 0.9 L/kg. Low plasma protein binding (10-25%). Crosses the placenta and enters breast milk.",
            "Metabolism": "Hepatic glucuronidation (about 50-60%) and sulfation (about 25-35%). A minor CYP2E1/1A2/3A4 route (about 5-10%) forms toxic NAPQI, normally detoxified by glutathione.",
            "Excretion": "About 90% renal as glucuronide/sulfate conjugates; less than 5% unchanged. t1/2 about 1.5-3 h.",
        },
        "pk": {"dose_mg": 1000, "f": 0.8, "vd_l_kg": 0.9, "t_half_h": 2.0, "tmax_h": 1.0, "note": "Typical adult single dose; first-pass effect is dose-dependent."},
        "pk_table": {"Bioavailability": "about 63-89%", "Tmax": "10-60 min (IR)", "Half-life": "1.5-3 h", "Protein binding": "10-25%", "Vd": "about 0.9 L/kg"},
        "tox": {
            "ld50": [{"species": "Rat", "route": "Oral", "value": 1944.0, "unit": "mg/kg", "source": "ChemIDplus via PubChem", "url": PUBCHEM_COMPOUND_PAGE + "1983#section=Toxicity"}],
            "clinical": [
                {"marker": "Hepatotoxicity", "detail": "Acute overdose can cause hepatic necrosis; N-acetylcysteine is the antidote. Risk rises with chronic alcohol use, fasting, malnutrition.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "acetaminophen"},
                {"marker": "Hidden duplication", "detail": "Many cold/flu and prescription combinations contain paracetamol.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "acetaminophen"},
            ],
        },
        "dosing": {
            "Adult": _d("500-1000 mg every 4-6 h", "Max 4 g/day (many authorities advise 3 g/day for routine use)", "Count all paracetamol-containing products. Seek urgent care for any suspected overdose, even if asymptomatic."),
            "Pediatric": _d("10-15 mg/kg per dose every 4-6 h (weight-based)", "Max 75 mg/kg/day, not exceeding 4 g/day", "Use the product's own weight/age chart and the supplied device. Concentrations of liquid products differ."),
            "Geriatric": _d("500-1000 mg every 6 h", "Consider 3 g/day or less in frail, underweight or hepatic impairment", "Preferred first-line analgesic in older adults, but reduce dose with liver disease, malnutrition or heavy alcohol use."),
        },
        "conditions": [],
    },
    "ibuprofen": {
        "name": "Ibuprofen",
        "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "cid": 3672, "drugbank": "DB01050", "dailymed": "ibuprofen",
        "aliases": ["advil", "motrin", "brufen"],
        "drug_class": "NSAID (propionic acid derivative); reversible non-selective COX inhibitor",
        "summary": "Ibuprofen reversibly inhibits COX-1 and COX-2 to reduce prostaglandin synthesis, giving analgesic, antipyretic and anti-inflammatory effects. The S-enantiomer is mainly active; partial R-to-S inversion occurs in vivo.",
        "indications": ["Pain and fever", "Dysmenorrhoea", "Inflammatory arthritis"],
        "pros": ["Effective anti-inflammatory analgesic", "Rapid onset", "OTC availability"],
        "cons": ["GI ulcer/bleed risk", "Renal impairment and fluid retention", "Cardiovascular thrombotic risk at high doses", "Interacts with anticoagulants, ACE inhibitors, diuretics"],
        "adme": {
            "Absorption": "Rapid oral absorption, bioavailability about 80-100%. Food slows rate, not extent.",
            "Distribution": "Highly protein-bound (about 99%), Vd about 0.1-0.2 L/kg. Appears in synovial fluid; small amounts in breast milk.",
            "Metabolism": "Hepatic oxidation by CYP2C9 (and CYP2C8) to inactive hydroxy and carboxy metabolites, followed by glucuronidation. CYP2C9 poor metabolisers have higher exposure.",
            "Excretion": "Mostly renal as metabolites (less than 1% unchanged). t1/2 about 2 h.",
        },
        "pk": {"dose_mg": 400, "f": 0.9, "vd_l_kg": 0.15, "t_half_h": 2.0, "tmax_h": 1.5, "note": "Typical adult single oral dose."},
        "pk_table": {"Bioavailability": "about 80-100%", "Tmax": "1-2 h", "Half-life": "about 2 h", "Protein binding": "about 99%", "Vd": "about 0.1-0.2 L/kg"},
        "tox": {
            "ld50": [{"species": "Rat", "route": "Oral", "value": 636.0, "unit": "mg/kg", "source": "ChemIDplus via PubChem", "url": PUBCHEM_COMPOUND_PAGE + "3672#section=Toxicity"}],
            "clinical": [
                {"marker": "GI bleeding / perforation (boxed warning)", "detail": "Can occur without warning, more common in older adults.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "ibuprofen"},
                {"marker": "Cardiovascular thrombotic events (boxed warning)", "detail": "Risk may increase with duration and dose. Contraindicated peri-CABG.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "ibuprofen"},
                {"marker": "Renal injury", "detail": "Risk higher with dehydration, heart failure, CKD, ACE-i/diuretic co-use.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "ibuprofen"},
            ],
        },
        "dosing": {
            "Adult": _d("200-400 mg every 4-6 h with food", "OTC max 1200 mg/day; prescription up to 3200 mg/day under supervision", "Use the lowest effective dose for the shortest time. Avoid in late pregnancy, active GI bleed, severe renal/heart failure."),
            "Pediatric": _d("5-10 mg/kg every 6-8 h (age 6 months and older)", "Max 40 mg/kg/day (and within OTC pack limits)", "Not for infants under 6 months without medical advice. Maintain hydration. Use the weight chart on the product."),
            "Geriatric": _d("Lowest effective dose, shortest duration", "Avoid chronic use where possible", "Beers Criteria: avoid chronic use in older adults. Monitor renal function, blood pressure, GI symptoms."),
        },
        "conditions": [],
    },
    "caffeine": {
        "name": "Caffeine",
        "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "cid": 2519, "drugbank": "DB00201", "dailymed": "caffeine",
        "aliases": ["1,3,7-trimethylxanthine", "coffee"],
        "drug_class": "Methylxanthine; CNS stimulant; adenosine receptor antagonist",
        "summary": "Caffeine antagonises adenosine A1/A2A receptors and inhibits phosphodiesterases at higher concentrations, increasing alertness and mildly raising heart rate and diuresis.",
        "indications": ["Fatigue / alertness (dietary)", "Apnoea of prematurity (as caffeine citrate, specialist use)", "Adjunct in some analgesic combinations"],
        "pros": ["Improves vigilance and reaction time", "Rapid, nearly complete absorption", "Well-characterised safety at moderate intake"],
        "cons": ["Anxiety, insomnia, tachycardia at high intake", "Dependence and withdrawal headache", "Interacts with CYP1A2 inhibitors (e.g., fluvoxamine, ciprofloxacin)"],
        "adme": {
            "Absorption": "Almost complete (about 100%) oral absorption within 30-120 min.",
            "Distribution": "Vd about 0.6 L/kg; distributes into total body water; crosses the blood-brain barrier and placenta. Protein binding about 36%.",
            "Metabolism": "Hepatic CYP1A2 N-demethylation to paraxanthine (about 80%), theobromine, theophylline. Smoking induces, pregnancy and some drugs inhibit CYP1A2.",
            "Excretion": "Renal, with 1-3% unchanged. t1/2 about 3-7 h (prolonged in pregnancy, liver disease, neonates).",
        },
        "pk": {"dose_mg": 200, "f": 1.0, "vd_l_kg": 0.6, "t_half_h": 5.0, "tmax_h": 1.0, "note": "About 2 cups of brewed coffee."},
        "pk_table": {"Bioavailability": "about 100%", "Tmax": "0.5-2 h", "Half-life": "about 3-7 h", "Protein binding": "about 36%", "Vd": "about 0.6 L/kg"},
        "tox": {
            "ld50": [{"species": "Rat", "route": "Oral", "value": 192.0, "unit": "mg/kg", "source": "ChemIDplus via PubChem", "url": PUBCHEM_COMPOUND_PAGE + "2519#section=Toxicity"}],
            "clinical": [
                {"marker": "Acute caffeine toxicity", "detail": "Tachyarrhythmia, seizures, vomiting. Concentrated powders/tablets pose high risk of accidental overdose.", "source": "PubChem toxicity section", "url": PUBCHEM_COMPOUND_PAGE + "2519#section=Toxicity"},
            ],
        },
        "dosing": {
            "Adult": _d("Up to 400 mg/day for most healthy adults", "FDA general guidance: about 400 mg/day", "Pregnancy: commonly advised 200 mg/day or less (ACOG). Reduce with anxiety disorders, arrhythmia, hypertension."),
            "Pediatric": _d("No therapeutic dose for routine use", "Children: commonly cited limit 2.5 mg/kg/day (Health Canada); adolescents: 100 mg/day or less (AAP guidance)", "Energy drinks and concentrated caffeine products are not recommended for children and adolescents."),
            "Geriatric": _d("Individualise; often lower tolerance", "Staying under 400 mg/day, usually less if sleep-disturbed", "Slower clearance with age/liver disease; may worsen insomnia, palpitations, urinary urgency."),
        },
        "conditions": [],
    },
    "metformin": {
        "name": "Metformin",
        "smiles": "CN(C)C(=N)N=C(N)N",
        "cid": 4091, "drugbank": "DB00331", "dailymed": "metformin",
        "aliases": ["glucophage", "metformin hydrochloride"],
        "drug_class": "Biguanide; first-line antihyperglycaemic for type 2 diabetes",
        "summary": "Metformin lowers hepatic glucose output (AMPK-dependent and independent effects, complex I inhibition) and improves peripheral insulin sensitivity without causing hypoglycaemia on its own.",
        "indications": ["Type 2 diabetes mellitus", "Prediabetes in selected high-risk adults (guideline-dependent)", "PCOS (off-label)"],
        "pros": ["Weight-neutral to mildly weight-reducing", "Low hypoglycaemia risk", "Inexpensive; long safety record", "Cardiovascular benefit data in overweight T2D"],
        "cons": ["GI intolerance (diarrhoea, nausea)", "Rare lactic acidosis, mainly with renal impairment", "Vitamin B12 deficiency with long-term use"],
        "adme": {
            "Absorption": "Absorbed mainly in the small intestine (PMAT/OCT1-mediated). Oral bioavailability about 50-60% (fasting). Food slightly reduces and delays absorption.",
            "Distribution": "Negligible plasma protein binding. Large apparent Vd, accumulating in gut, liver and kidney.",
            "Metabolism": "Not metabolised; no known metabolites.",
            "Excretion": "Excreted unchanged in urine by glomerular filtration and tubular secretion (OCT2/MATE1/MATE2-K). Plasma t1/2 about 6 h; accumulates if eGFR is low.",
        },
        "pk": {"dose_mg": 500, "f": 0.55, "vd_l_kg": 3.0, "t_half_h": 6.2, "tmax_h": 2.5, "note": "Immediate-release tablet; Vd is highly variable across sources."},
        "pk_table": {"Bioavailability": "about 50-60%", "Tmax": "about 2.5 h (IR)", "Half-life": "about 6 h", "Protein binding": "negligible", "Vd": "about 60-280 L (wide range)"},
        "tox": {
            "ld50": [{"species": "Rat", "route": "Oral", "value": 1000.0, "unit": "mg/kg", "source": "Manufacturer SDS (varies by salt/vendor)", "url": PUBCHEM_COMPOUND_PAGE + "4091#section=Toxicity"}],
            "clinical": [
                {"marker": "Lactic acidosis (boxed warning)", "detail": "Risk factors: renal impairment, hypoxia, sepsis, excessive alcohol, hepatic impairment, iodinated contrast.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "metformin"},
                {"marker": "Renal contraindication", "detail": "Contraindicated at eGFR below 30 mL/min/1.73 m2; initiation not recommended at 30-45.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "metformin"},
            ],
        },
        "dosing": {
            "Adult": _d("500 mg once/twice daily with meals, titrated", "IR max 2550 mg/day; ER usually 2000 mg/day", "Check renal function before and during treatment. Hold around iodinated contrast and acute illness per clinician."),
            "Pediatric": _d("Age 10 years and older (T2D): 500 mg twice daily initially", "Max 2000 mg/day (IR)", "Paediatric endocrinology supervision recommended. Monitor growth, GI tolerance, renal function."),
            "Geriatric": _d("Start low; titrate by renal function", "Do not up-titrate to maximum in frail elderly", "Avoid if eGFR below 30. Monitor B12 and renal function at least annually."),
        },
        "conditions": ["diabetes"],
    },
    "atorvastatin": {
        "name": "Atorvastatin",
        "smiles": "CC(C)C1=C(C(=C(N1CC[C@H](C[C@H](CC(=O)O)O)O)C2=CC=C(C=C2)F)C3=CC=CC=C3)C(=O)NC4=CC=CC=C4",
        "cid": 60823, "drugbank": "DB01076", "dailymed": "atorvastatin",
        "aliases": ["lipitor", "atorvastatin calcium"],
        "drug_class": "HMG-CoA reductase inhibitor (statin)",
        "summary": "Atorvastatin competitively inhibits hepatic HMG-CoA reductase, lowering LDL-cholesterol and cardiovascular events. It also has pleiotropic anti-inflammatory and endothelial effects.",
        "indications": ["Primary and secondary cardiovascular prevention", "Hyperlipidaemia", "Heterozygous familial hypercholesterolaemia (10 years and older)"],
        "pros": ["Strong LDL-C reduction", "Large outcome-trial evidence base", "No renal dose adjustment"],
        "cons": ["Myalgia / rare rhabdomyolysis", "Raised transaminases", "Small increase in new-onset diabetes risk", "CYP3A4 interactions (e.g., clarithromycin, grapefruit in large amounts)"],
        "adme": {
            "Absorption": "Rapid absorption; absolute bioavailability about 12% due to pre-systemic clearance. Food lowers rate but not LDL-C effect.",
            "Distribution": "Vd about 381 L; plasma protein binding 98% or more; hepatic uptake via OATP1B1.",
            "Metabolism": "CYP3A4 to active ortho- and para-hydroxylated metabolites (about 70% of HMG-CoA reductase activity).",
            "Excretion": "Biliary/faecal; less than 2% recovered in urine. t1/2 about 14 h; inhibitory activity lasts 20-30 h.",
        },
        "pk": {"dose_mg": 20, "f": 0.12, "vd_l_kg": 5.4, "t_half_h": 14.0, "tmax_h": 1.5, "note": "Illustrative: parent drug only; active metabolites extend pharmacodynamic effect."},
        "pk_table": {"Bioavailability": "about 12% (systemic)", "Tmax": "1-2 h", "Half-life": "about 14 h (activity 20-30 h)", "Protein binding": "98% or more", "Vd": "about 381 L"},
        "tox": {
            "ld50": [],
            "clinical": [
                {"marker": "Myopathy / rhabdomyolysis", "detail": "Risk increases with high dose, age 65 and older, hypothyroidism, renal impairment and interacting drugs.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "atorvastatin"},
                {"marker": "Hepatic enzyme elevation", "detail": "Baseline and as-clinically-indicated liver enzyme monitoring.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "atorvastatin"},
                {"marker": "Acute LD50", "detail": "No curated value is listed here. Consult the manufacturer SDS and PubChem toxicity section.", "source": "PubChem toxicity section", "url": PUBCHEM_COMPOUND_PAGE + "60823#section=Toxicity"},
            ],
        },
        "dosing": {
            "Adult": _d("10-80 mg once daily (start 10-20 mg; 40-80 mg for high-intensity)", "Max 80 mg/day", "Avoid in pregnancy and active liver disease. Report unexplained muscle pain or dark urine."),
            "Pediatric": _d("Age 10-17 years (HeFH): 10 mg once daily", "Titrate to max 20 mg/day", "Specialist-supervised only. Not studied in children under 10."),
            "Geriatric": _d("No age-based adjustment; start moderate", "Use lowest effective dose to reach target", "Higher myopathy risk. Review interacting drugs (macrolides, azole antifungals, diltiazem)."),
        },
        "conditions": ["cardiovascular"],
    },
    "amoxicillin": {
        "name": "Amoxicillin",
        "smiles": "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)[C@@H](C3=CC=C(C=C3)O)N)C(=O)O)C",
        "cid": 33613, "drugbank": "DB01060", "dailymed": "amoxicillin",
        "aliases": ["amoxil", "amoxycillin"],
        "drug_class": "Aminopenicillin; beta-lactam antibiotic",
        "summary": "Amoxicillin inhibits bacterial cell-wall synthesis by binding penicillin-binding proteins. It is active against many Gram-positive and some Gram-negative organisms but is destroyed by beta-lactamases.",
        "indications": ["Otitis media, sinusitis, pharyngitis (streptococcal)", "Community-acquired pneumonia", "H. pylori eradication (combination)", "Urinary tract infections (susceptible organisms)"],
        "pros": ["Good oral absorption", "Broad first-line spectrum", "Generally well tolerated", "Safe in pregnancy per guidelines"],
        "cons": ["Penicillin allergy / anaphylaxis", "Beta-lactamase resistance", "Diarrhoea, rash (esp. in EBV infection)", "Contributes to AMR if misused"],
        "adme": {
            "Absorption": "Well absorbed orally (about 70-90%), stable in gastric acid, and unaffected by food.",
            "Distribution": "Protein binding about 17-20%; Vd about 0.3 L/kg; penetrates lung, middle ear, tonsil; limited CSF penetration unless meninges inflamed.",
            "Metabolism": "Minimal; about 10-25% converted to inactive penicilloic acid.",
            "Excretion": "Renal (tubular secretion + filtration); about 60% unchanged in urine within 6 h. t1/2 about 1-1.4 h (longer in renal impairment). Probenecid raises levels.",
        },
        "pk": {"dose_mg": 500, "f": 0.8, "vd_l_kg": 0.3, "t_half_h": 1.2, "tmax_h": 1.5, "note": "Typical adult oral dose."},
        "pk_table": {"Bioavailability": "about 70-90%", "Tmax": "1-2 h", "Half-life": "about 1-1.4 h", "Protein binding": "about 17-20%", "Vd": "about 0.3 L/kg"},
        "tox": {
            "ld50": [{"species": "Rat", "route": "Oral", "value": 5000.0, "unit": "mg/kg (greater than)", "source": "Manufacturer SDS (vendor-dependent)", "url": PUBCHEM_COMPOUND_PAGE + "33613#section=Toxicity"}],
            "clinical": [
                {"marker": "Hypersensitivity / anaphylaxis", "detail": "Cross-reactivity within beta-lactams. Document allergy history.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "amoxicillin"},
                {"marker": "C. difficile-associated diarrhoea", "detail": "May occur during or weeks after therapy.", "source": "FDA label (DailyMed)", "url": DAILYMED_SEARCH + "amoxicillin"},
            ],
        },
        "dosing": {
            "Adult": _d("500 mg every 12 h or 250-500 mg every 8 h (indication-dependent)", "Higher regimens (e.g., 875 mg every 12 h) per label/clinician", "Complete the prescribed course. Adjust for renal impairment."),
            "Pediatric": _d("20-45 mg/kg/day in divided doses every 8-12 h (older than 3 months)", "High-dose regimens up to 80-90 mg/kg/day only by clinician (e.g., otitis media). Under 3 months: max 30 mg/kg/day every 12 h", "Weight-based dosing; reconstitute suspensions exactly per label. Never share antibiotics."),
            "Geriatric": _d("Standard adult dose adjusted for renal function", "Reduce frequency if CrCl is low", "Higher risk of C. difficile and rash. Check drug interactions (warfarin)."),
        },
        "conditions": ["amr"],
    },
}

# --------------------------------------------------------------------------- #
# Global Health knowledge base (educational)
# --------------------------------------------------------------------------- #
GLOBAL_HEALTH: Dict[str, Dict[str, Any]] = {
    "cardiovascular": {
        "title": "Cardiovascular Diseases", "icon": "🫀",
        "burden": "Cardiovascular diseases are the leading cause of death worldwide: WHO estimates about 17.9 million deaths per year (about 32% of global deaths). Most are heart attacks and strokes, and more than three-quarters occur in low- and middle-income countries.",
        "treatments": [
            {"cls": "Statins", "examples": "Atorvastatin, rosuvastatin", "mech": "Inhibit HMG-CoA reductase, lowering LDL-C", "note": "Foundation of primary and secondary prevention"},
            {"cls": "Antiplatelets", "examples": "Aspirin, clopidogrel", "mech": "COX-1 / P2Y12 inhibition reduces thrombosis", "note": "Mainly secondary prevention; bleeding risk"},
            {"cls": "ACE inhibitors / ARBs", "examples": "Lisinopril, losartan", "mech": "RAAS blockade lowers BP and remodelling", "note": "Hypertension, heart failure, CKD"},
            {"cls": "Beta-blockers", "examples": "Metoprolol, bisoprolol", "mech": "Reduce heart rate and myocardial oxygen demand", "note": "Post-MI, heart failure, arrhythmia"},
            {"cls": "Calcium-channel blockers / diuretics", "examples": "Amlodipine, hydrochlorothiazide", "mech": "Vasodilation / natriuresis", "note": "First-line antihypertensives"},
        ],
        "alternatives": ["Ezetimibe (cholesterol absorption inhibitor)", "PCSK9 monoclonal antibodies and inclisiran (siRNA)", "Bempedoic acid for statin-intolerant patients", "SGLT2 inhibitors for heart-failure benefit", "Fixed-dose combination 'polypill' strategies to improve adherence"],
        "prevention": ["Blood-pressure screening and control", "Tobacco cessation and reduced sodium intake", "Regular physical activity and healthy diet", "Weight management and lipid/glucose screening", "Policy measures: trans-fat elimination, salt reduction"],
        "links": [("WHO Fact Sheet: Cardiovascular diseases", "https://www.who.int/news-room/fact-sheets/detail/cardiovascular-diseases-(cvds)")],
        "library": ["atorvastatin", "aspirin"],
    },
    "diabetes": {
        "title": "Type 2 Diabetes", "icon": "🩸",
        "burden": "WHO reports about 830 million people were living with diabetes in 2022, most with type 2 diabetes, and treatment coverage remains low in many low- and middle-income countries. Complications include cardiovascular disease, kidney failure, blindness and amputation.",
        "treatments": [
            {"cls": "Biguanides", "examples": "Metformin", "mech": "Reduces hepatic glucose output; improves insulin sensitivity", "note": "Usual first-line, low cost"},
            {"cls": "SGLT2 inhibitors", "examples": "Empagliflozin, dapagliflozin", "mech": "Increase urinary glucose excretion", "note": "Cardio-renal protection; genital infections, DKA risk"},
            {"cls": "GLP-1 receptor agonists", "examples": "Semaglutide, liraglutide", "mech": "Incretin mimetic; weight loss, glucose-dependent insulin secretion", "note": "CV benefit; GI effects"},
            {"cls": "DPP-4 inhibitors", "examples": "Sitagliptin, linagliptin", "mech": "Prolong endogenous incretin action", "note": "Weight-neutral, well tolerated"},
            {"cls": "Sulfonylureas / insulin", "examples": "Glimepiride; basal/bolus insulin", "mech": "Increase insulin secretion / replace insulin", "note": "Hypoglycaemia and weight gain risk"},
        ],
        "alternatives": ["Dual GIP/GLP-1 agonists (e.g., tirzepatide)", "Thiazolidinediones (pioglitazone) in selected patients", "Emerging once-weekly basal insulins", "Continuous glucose monitoring and digital therapeutics"],
        "prevention": ["Intensive lifestyle programmes (the US Diabetes Prevention Program showed about 58% lower incidence vs placebo)", "Screening of at-risk adults", "Weight management and physical activity", "Metformin for selected high-risk prediabetes"],
        "links": [("WHO Fact Sheet: Diabetes", "https://www.who.int/news-room/fact-sheets/detail/diabetes")],
        "library": ["metformin"],
    },
    "pandemics": {
        "title": "Infectious Pandemics", "icon": "🦠",
        "burden": "Emerging respiratory viruses (influenza, coronaviruses) and other pathogens can spread rapidly across borders. COVID-19 showed the value of rapid diagnostics, antivirals, platform vaccines and equitable access.",
        "treatments": [
            {"cls": "Neuraminidase / endonuclease inhibitors", "examples": "Oseltamivir, baloxavir", "mech": "Block influenza virus release / replication", "note": "Best started early"},
            {"cls": "Protease inhibitors", "examples": "Nirmatrelvir (with ritonavir)", "mech": "Inhibit SARS-CoV-2 main protease", "note": "Major drug-interaction profile"},
            {"cls": "Polymerase inhibitors", "examples": "Remdesivir", "mech": "RNA-dependent RNA polymerase inhibition", "note": "Intravenous; hospital use"},
            {"cls": "Immunomodulators", "examples": "Dexamethasone", "mech": "Reduces hyperinflammation", "note": "Benefit in severe, oxygen-requiring COVID-19"},
            {"cls": "Monoclonal antibodies", "examples": "Pathogen-specific mAbs", "mech": "Neutralise virus", "note": "Variant escape is a challenge"},
        ],
        "alternatives": ["Broad-spectrum antivirals", "Host-directed therapies", "Platform vaccines (mRNA, viral vector, subunit)", "Pan-family vaccine candidates and rapid-response manufacturing"],
        "prevention": ["Genomic surveillance and early warning", "Vaccination and equitable distribution", "Ventilation, hand hygiene, masks in high-risk settings", "One Health: zoonotic spillover monitoring", "Strong primary care and supply chains"],
        "links": [("WHO Fact Sheets (index)", "https://www.who.int/news-room/fact-sheets")],
        "library": [],
    },
    "amr": {
        "title": "Antimicrobial Resistance (AMR)", "icon": "🧫",
        "burden": "The GRAM study (Lancet, 2024) estimated 1.14 million deaths directly attributable to bacterial AMR and 4.71 million associated deaths in 2021, projected to rise sharply by 2050 without action.",
        "treatments": [
            {"cls": "Beta-lactams", "examples": "Amoxicillin, ceftriaxone", "mech": "Inhibit cell-wall synthesis (PBPs)", "note": "Resistance via beta-lactamases"},
            {"cls": "Beta-lactam / inhibitor combos", "examples": "Amoxicillin-clavulanate, ceftazidime-avibactam", "mech": "Protect the beta-lactam from enzymes", "note": "Reserve agents for MDR Gram-negatives"},
            {"cls": "Carbapenems", "examples": "Meropenem, ertapenem", "mech": "Broad beta-lactam", "note": "Carbapenem-resistant organisms are critical priority"},
            {"cls": "Glycopeptides", "examples": "Vancomycin", "mech": "Bind D-Ala-D-Ala, blocking cell-wall cross-linking", "note": "MRSA and Enterococcus"},
            {"cls": "Fluoroquinolones", "examples": "Ciprofloxacin, levofloxacin", "mech": "DNA gyrase / topoisomerase IV inhibition", "note": "Tendinopathy, QT, resistance risk"},
        ],
        "alternatives": ["Siderophore cephalosporins (e.g., cefiderocol)", "Bacteriophage therapy and phage lysins", "Antimicrobial peptides and anti-virulence agents", "Novel Gram-negative chemotypes in development", "Rapid diagnostics to cut empirical broad-spectrum use"],
        "prevention": ["Antimicrobial stewardship in hospitals and community", "Infection prevention and control, WASH", "Vaccines (pneumococcal, typhoid conjugate, etc.)", "Reducing non-human antimicrobial use", "Surveillance (e.g., WHO GLASS)"],
        "links": [("WHO Fact Sheet: Antimicrobial resistance", "https://www.who.int/news-room/fact-sheets/detail/antimicrobial-resistance")],
        "library": ["amoxicillin"],
    },
}

# --------------------------------------------------------------------------- #
# Library helpers
# --------------------------------------------------------------------------- #
def library_keys() -> List[str]:
    return list(LIBRARY.keys())


def library_name(key: str) -> str:
    return LIBRARY[key]["name"]


def find_library_key(query: str) -> Optional[str]:
    q = (query or "").strip().lower()
    if not q:
        return None
    for key, entry in LIBRARY.items():
        names = [key, entry["name"].lower()] + [a.lower() for a in entry.get("aliases", [])]
        if q in names or q == entry["name"].split("(")[0].strip().lower():
            return key
    return None


@lru_cache(maxsize=1)
def _library_canonical_map() -> Dict[str, str]:
    out = {}
    for key, entry in LIBRARY.items():
        mol = Chem.MolFromSmiles(entry["smiles"])
        if mol is not None:
            out[Chem.MolToSmiles(mol, isomericSmiles=False)] = key
    return out


def match_library(mol: Chem.Mol) -> Optional[str]:
    """Match a molecule to the curated library (ignoring stereochemistry)."""
    return _library_canonical_map().get(Chem.MolToSmiles(mol, isomericSmiles=False))


# --------------------------------------------------------------------------- #
# PubChem
# --------------------------------------------------------------------------- #
def _get_json(url: str) -> Optional[dict]:
    try:
        r = requests.get(url, headers=HTTP_HEADERS, timeout=HTTP_TIMEOUT)
        return r.json() if r.status_code == 200 else None
    except (requests.RequestException, ValueError):
        return None


def _post_json(url: str, data: dict) -> Optional[dict]:
    try:
        r = requests.post(url, data=data, headers=HTTP_HEADERS, timeout=HTTP_TIMEOUT)
        return r.json() if r.status_code == 200 else None
    except (requests.RequestException, ValueError):
        return None


def _parse_pubchem_props(payload: Optional[dict]) -> Optional[Dict[str, Any]]:
    try:
        p = payload["PropertyTable"]["Properties"][0]
    except (TypeError, KeyError, IndexError):
        return None
    smiles = next((p[k] for k in ("SMILES", "IsomericSMILES", "ConnectivitySMILES", "CanonicalSMILES") if p.get(k)), None)
    if not smiles:
        return None
    return {
        "cid": p.get("CID"), "title": p.get("Title"), "iupac": p.get("IUPACName"),
        "formula": p.get("MolecularFormula"), "mw": p.get("MolecularWeight"),
        "exact_mass": p.get("ExactMass"), "xlogp": p.get("XLogP"), "tpsa": p.get("TPSA"),
        "smiles": smiles, "inchikey": p.get("InChIKey"),
    }


@lru_cache(maxsize=256)
def fetch_pubchem_by_name(name: str) -> Optional[Dict[str, Any]]:
    quoted = requests.utils.quote(name.strip(), safe="")
    for props in PUBCHEM_PROPERTY_SETS:
        data = _get_json(f"{PUBCHEM_REST}/compound/name/{quoted}/property/{props}/JSON")
        parsed = _parse_pubchem_props(data)
        if parsed:
            return parsed
    return None


@lru_cache(maxsize=256)
def fetch_pubchem_by_smiles(smiles: str) -> Optional[Dict[str, Any]]:
    for props in PUBCHEM_PROPERTY_SETS:
        data = _post_json(f"{PUBCHEM_REST}/compound/smiles/property/{props}/JSON", {"smiles": smiles})
        parsed = _parse_pubchem_props(data)
        if parsed:
            return parsed
    return None


@lru_cache(maxsize=256)
def fetch_pubchem_description(cid: int) -> Optional[Dict[str, str]]:
    data = _get_json(f"{PUBCHEM_REST}/compound/cid/{int(cid)}/description/JSON")
    try:
        for info in data["InformationList"]["Information"]:
            if info.get("Description"):
                return {"text": info["Description"], "source": info.get("DescriptionSourceName", "PubChem")}
    except (TypeError, KeyError):
        pass
    return None


# --------------------------------------------------------------------------- #
# RDKit: parsing + descriptors
# --------------------------------------------------------------------------- #
def parse_smiles(smiles: str) -> Chem.Mol:
    s = (smiles or "").strip()
    if not s:
        raise ValueError("SMILES string is empty.")
    if len(s) > MAX_SMILES_LEN:
        raise ValueError(f"SMILES too long (limit {MAX_SMILES_LEN} characters).")
    mol = Chem.MolFromSmiles(s)
    if mol is None:
        raise ValueError("Invalid SMILES: RDKit could not parse or sanitise this structure.")
    return mol


def compute_descriptors(mol: Chem.Mol) -> Dict[str, Any]:
    heavy = mol.GetNumHeavyAtoms()
    arom = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
    props = {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "mw": Descriptors.MolWt(mol),
        "exact_mass": Descriptors.ExactMolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "mr": Crippen.MolMR(mol),
        "tpsa": rdMolDescriptors.CalcTPSA(mol),
        "hbd": Lipinski.NumHDonors(mol),
        "hba": Lipinski.NumHAcceptors(mol),
        "rotb": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "rings": rdMolDescriptors.CalcNumRings(mol),
        "arom_rings": rdMolDescriptors.CalcNumAromaticRings(mol),
        "heavy_atoms": heavy,
        "total_atoms": Chem.AddHs(mol).GetNumAtoms(),
        "fsp3": rdMolDescriptors.CalcFractionCSP3(mol),
        "charge": Chem.GetFormalCharge(mol),
        "stereocenters": len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)),
        "qed": QED.qed(mol),
    }
    ap = arom / heavy if heavy else 0.0
    props["esol_logs"] = 0.16 - 0.63 * props["logp"] - 0.0062 * props["mw"] + 0.066 * props["rotb"] - 0.74 * ap
    return props


def esol_class(logs: float) -> str:
    if logs > 0: return "Highly soluble"
    if logs > -2: return "Very soluble"
    if logs > -4: return "Soluble"
    if logs > -6: return "Moderately soluble"
    if logs > -10: return "Poorly soluble"
    return "Insoluble"


def lipinski_report(p: Dict[str, Any]) -> Dict[str, Any]:
    rules = [
        ("Molecular weight", p["mw"], "500 Da or less", p["mw"] <= 500),
        ("LogP (Crippen)", p["logp"], "5 or less", p["logp"] <= 5),
        ("H-bond donors", p["hbd"], "5 or less", p["hbd"] <= 5),
        ("H-bond acceptors", p["hba"], "10 or less", p["hba"] <= 10),
    ]
    violations = sum(1 for r in rules if not r[3])
    return {"rules": rules, "violations": violations, "compliant": violations <= 1,
            "label": "Lipinski compliant" if violations <= 1 else "Lipinski non-compliant"}


def druglikeness_filters(p: Dict[str, Any]) -> List[Dict[str, Any]]:
    veber = p["rotb"] <= 10 and p["tpsa"] <= 140
    egan = p["tpsa"] <= 131.6 and -1 <= p["logp"] <= 5.88
    ghose = 160 <= p["mw"] <= 480 and -0.4 <= p["logp"] <= 5.6 and 40 <= p["mr"] <= 130 and 20 <= p["total_atoms"] <= 70
    return [
        {"name": "Veber", "criteria": "RotB 10 or fewer and TPSA 140 or less", "pass": veber},
        {"name": "Egan", "criteria": "TPSA 131.6 or less and LogP between -1 and 5.88", "pass": egan},
        {"name": "Ghose", "criteria": "MW 160-480, LogP -0.4 to 5.6, MR 40-130, atoms 20-70", "pass": ghose},
    ]


def heuristic_adme(p: Dict[str, Any]) -> Dict[str, str]:
    """Rule-of-thumb ADME flags. NOT a substitute for experimental or ML-based predictions."""
    return {
        "GI absorption (Egan)": "Likely high" if (p["tpsa"] <= 131.6 and p["logp"] <= 5.88) else "Likely low",
        "BBB permeation (heuristic)": "Likely" if (p["tpsa"] <= 79 and p["mw"] <= 450 and p["hbd"] <= 3) else "Unlikely",
        "Aqueous solubility (ESOL)": f"{p['esol_logs']:.2f} logS, {esol_class(p['esol_logs'])}",
        "Oral bioavailability (Veber)": "Favourable" if (p["rotb"] <= 10 and p["tpsa"] <= 140) else "Possibly reduced",
    }


@lru_cache(maxsize=4)
def _filter_catalog(name: str) -> FilterCatalog:
    params = FilterCatalogParams()
    params.AddCatalog(getattr(FilterCatalogParams.FilterCatalogs, name))
    return FilterCatalog(params)


def structural_alerts(mol: Chem.Mol) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for label, cat in (("PAINS", "PAINS"), ("Brenk", "BRENK"), ("NIH", "NIH")):
        try:
            matches = _filter_catalog(cat).GetMatches(mol)
            out[label] = sorted({m.GetDescription() for m in matches})
        except Exception:
            out[label] = []
    return out


def elemental_composition(mol: Chem.Mol) -> List[Dict[str, Any]]:
    molh = Chem.AddHs(mol)
    pt = Chem.GetPeriodicTable()
    counts: Dict[str, int] = {}
    for a in molh.GetAtoms():
        counts[a.GetSymbol()] = counts.get(a.GetSymbol(), 0) + 1
    total = sum(pt.GetAtomicWeight(sym) * n for sym, n in counts.items())
    rows = [{"Symbol": sym, "Element": pt.GetElementName(pt.GetAtomicNumber(sym)),
             "Atomic No.": pt.GetAtomicNumber(sym), "Atoms": n,
             "Atomic mass": round(pt.GetAtomicWeight(sym), 3),
             "Mass %": round(100 * pt.GetAtomicWeight(sym) * n / total, 2)} for sym, n in counts.items()]
    return sorted(rows, key=lambda r: -r["Mass %"])


def ghs_acute_category(ld50_mg_kg: float) -> str:
    """GHS acute oral toxicity category from LD50 (mg/kg)."""
    if ld50_mg_kg <= 5: return "Category 1 (fatal)"
    if ld50_mg_kg <= 50: return "Category 2 (fatal)"
    if ld50_mg_kg <= 300: return "Category 3 (toxic)"
    if ld50_mg_kg <= 2000: return "Category 4 (harmful)"
    if ld50_mg_kg <= 5000: return "Category 5 (may be harmful)"
    return "Not classified"


def safe_inchikey(mol: Chem.Mol) -> Optional[str]:
    try:
        from rdkit.Chem import inchi
        return inchi.MolToInchiKey(mol)
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# 3D + 2D generation
# --------------------------------------------------------------------------- #
def generate_3d_sdf(smiles: str, optimize: bool = True, n_confs: int = 10, seed: int = 42, max_iters: int = 2000) -> Dict[str, Any]:
    """ETKDGv3 conformer search + UFF optimisation. Returns the lowest-energy conformer as an SDF block."""
    mol = parse_smiles(smiles)
    if mol.GetNumHeavyAtoms() > MAX_HEAVY_ATOMS_3D:
        raise ValueError(f"Molecule too large for interactive 3D generation (over {MAX_HEAVY_ATOMS_3D} heavy atoms).")
    molh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    n_confs = max(1, int(n_confs))
    cids = list(AllChem.EmbedMultipleConfs(molh, numConfs=n_confs, params=params))
    if not cids:
        params.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(molh, numConfs=n_confs, params=params))
    if not cids:
        raise RuntimeError("3D embedding failed for this structure.")

    energies: List[Optional[float]] = [None] * len(cids)
    converged_flags: List[Optional[bool]] = [None] * len(cids)
    best = 0
    if optimize and AllChem.UFFHasAllMoleculeParams(molh):
        results = AllChem.UFFOptimizeMoleculeConfs(molh, maxIters=max_iters, numThreads=0)
        converged_flags = [r[0] == 0 for r in results]
        energies = [float(r[1]) for r in results]
        best = min(range(len(energies)), key=energies.__getitem__)
    sdf = Chem.MolToMolBlock(molh, confId=cids[best])
    return {"sdf": sdf, "n_conformers": len(cids), "best_energy": energies[best],
            "energies": energies, "converged": converged_flags[best],
            "n_atoms": molh.GetNumAtoms(), "optimized": bool(optimize and energies[best] is not None)}


def mol_to_svg(smiles: str, width: int = 480, height: int = 360, dark: bool = True) -> str:
    mol = parse_smiles(smiles)
    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    if dark:
        try:
            rdMolDraw2D.SetDarkMode(opts)
        except Exception:
            pass
    opts.clearBackground = False
    opts.bondLineWidth = 2
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg[svg.find("<svg"):] if "<svg" in svg else svg


# --------------------------------------------------------------------------- #
# Pharmacokinetic illustration (one-compartment, first-order absorption)
# --------------------------------------------------------------------------- #
def _solve_ka(ke: float, tmax: float) -> float:
    """Solve tmax = ln(ka/ke)/(ka-ke) for ka by bisection; fall back if infeasible."""
    f = lambda ka: math.log(ka / ke) / (ka - ke)
    lo, hi = ke * (1 + 1e-6), ke * 1e4
    if tmax >= f(lo) or tmax <= f(hi):
        return ke * 3.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > tmax:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def simulate_pk(dose_mg: float, f: float, vd_l_kg: float, t_half_h: float, tmax_h: float,
                weight_kg: float = 70.0, points: int = 240) -> Dict[str, List[float]]:
    ke = math.log(2) / t_half_h
    ka = _solve_ka(ke, tmax_h)
    v = vd_l_kg * weight_kg
    horizon = min(48.0, max(8.0, 6 * t_half_h))
    ts = [horizon * i / (points - 1) for i in range(points)]
    cs = []
    for t in ts:
        if abs(ka - ke) < 1e-9:
            c = f * dose_mg / v * ke * t * math.exp(-ke * t)
        else:
            c = f * dose_mg * ka / (v * (ka - ke)) * (math.exp(-ke * t) - math.exp(-ka * t))
        cs.append(max(c, 0.0))
    return {"t": ts, "c": cs, "ka": ka, "ke": ke, "cmax": max(cs), "tmax_sim": ts[cs.index(max(cs))]}


# --------------------------------------------------------------------------- #
# Record builder / resolver
# --------------------------------------------------------------------------- #
def build_record(smiles: str, name: Optional[str] = None, source: str = "User SMILES",
                 pubchem: Optional[Dict[str, Any]] = None, library_key: Optional[str] = None,
                 cid: Optional[int] = None) -> Dict[str, Any]:
    mol = parse_smiles(smiles)
    if library_key is None:
        library_key = match_library(mol)
    curated = LIBRARY.get(library_key) if library_key else None
    if curated:
        name = name or curated["name"]
        cid = cid or curated.get("cid")
    props = compute_descriptors(mol)
    return {
        "name": name or "Custom molecule", "smiles": smiles.strip(), "canonical": Chem.MolToSmiles(mol),
        "source": source, "cid": cid, "pubchem": pubchem, "library_key": library_key, "curated": curated,
        "props": props, "lipinski": lipinski_report(props), "filters": druglikeness_filters(props),
        "heuristic_adme": heuristic_adme(props), "alerts": structural_alerts(mol),
        "elements": elemental_composition(mol), "inchikey": (pubchem or {}).get("inchikey") or safe_inchikey(mol),
        "iupac": (pubchem or {}).get("iupac"),
    }


def resolve_molecule(mode: str, query: str) -> Dict[str, Any]:
    """mode: 'library' | 'name' | 'smiles'."""
    q = (query or "").strip()
    if not q:
        raise ValueError("Please enter a value to analyse.")
    if mode == "library":
        if q not in LIBRARY:
            raise ValueError("Unknown library entry.")
        e = LIBRARY[q]
        return build_record(e["smiles"], e["name"], "Curated library", library_key=q, cid=e["cid"])
    if mode == "name":
        key = find_library_key(q)
        if key:
            e = LIBRARY[key]
            return build_record(e["smiles"], e["name"], "Curated library", library_key=key, cid=e["cid"])
        pc = fetch_pubchem_by_name(q)
        if not pc:
            raise ValueError(f"'{q}' was not found on PubChem (or PubChem is unreachable). Try the SMILES input.")
        return build_record(pc["smiles"], pc.get("title") or q, "PubChem (live)", pubchem=pc, cid=pc.get("cid"))
    if mode == "smiles":
        rec = build_record(q, None, "User SMILES")
        pc = None if rec["library_key"] else fetch_pubchem_by_smiles(rec["canonical"])
        if pc:
            rec = build_record(q, pc.get("title"), "User SMILES + PubChem match", pubchem=pc, cid=pc.get("cid"))
        return rec
    raise ValueError("Unknown input mode.")
