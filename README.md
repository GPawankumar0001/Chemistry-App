# 🧪 SMILES in Chemistry

An interactive Streamlit dashboard for exploring chemical SMILES strings, 3D molecular structures, drug properties and pros/cons.

**Institution:** SRM Institute of Science and Technology (SRM IST), Tiruchirappalli

## Features

- **Home / Dashboard:** overview, quick search and a short guide to SMILES notation
- **Molecule & Drug Explorer:** search by name or SMILES; formula, molecular weight, IUPAC name, LogP, H-bond donors/acceptors, TPSA, Lipinski check and elemental table
- **3D Structure Viewer:** rotatable py3Dmol models (ball-and-stick, stick, space-filling, wireframe) with SDF download
- **Comparative Analysis & Pros/Cons:** side-by-side properties, charts, uses, mechanism of action, benefits and drawbacks
- **About & Credits:** project and institution details

Built-in molecules: Aspirin, Caffeine, Paracetamol, Benzene, Glucose, Ibuprofen. Any other name is looked up on PubChem; any valid SMILES works directly.

## Project structure

```
smiles-in-chemistry/
├── app.py
├── utils.py
├── requirements.txt
├── README.md
└── assets/
    └── style.css
```

## Run locally

```bash
git clone https://github.com/<your-username>/smiles-in-chemistry.git
cd smiles-in-chemistry
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

Python 3.10 or newer is recommended.

## Deploy on Streamlit Community Cloud

1. Push this folder to a public GitHub repository.
2. Go to [share.streamlit.io](https://share.streamlit.io) and sign in with GitHub.
3. Choose **New app**, select the repository and branch, and set the main file to `app.py`.
4. Select **Deploy**. Streamlit installs everything in `requirements.txt` automatically.

## Disclaimer

Drug summaries are simplified for teaching and are not medical advice.

## Credits

Created at SRM Institute of Science and Technology (SRM IST), Tiruchirappalli. Powered by RDKit, PubChem, py3Dmol and Streamlit.
