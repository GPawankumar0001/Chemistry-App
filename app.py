"""SMILES in Chemistry - SRM Institute of Science and Technology (SRM IST), Tiruchirappalli."""
from pathlib import Path

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

import utils

INST = "SRM Institute of Science and Technology (SRM IST), Tiruchirappalli"

st.set_page_config(page_title="SMILES in Chemistry | SRM IST Tiruchirappalli", page_icon="🧪", layout="wide")
css = Path(__file__).parent / "assets" / "style.css"
st.markdown(f"<style>{css.read_text(encoding='utf-8')}</style>", unsafe_allow_html=True)

st.session_state.setdefault("query", "Aspirin")

# ---------- sidebar ----------
PAGES = ["🏠 Home / Dashboard", "🧪 Molecule & Drug Explorer", "🧊 3D Structure Viewer",
         "📊 Comparative Analysis & Pros/Cons", "ℹ️ About & Credits"]
with st.sidebar:
    st.markdown('<div class="side-brand">SMILES in Chemistry</div>'
                f'<div class="side-inst">{INST}</div>', unsafe_allow_html=True)
    st.divider()
    page = st.radio("Navigate", PAGES, label_visibility="collapsed")
    st.divider()
    st.caption("Quick load")
    for name in utils.PRESETS:
        if st.button(name, key=f"side_{name}", use_container_width=True):
            st.session_state["query"] = name


def current():
    """Resolve the active query once per rerun."""
    info = utils.resolve(st.session_state["query"])
    if info is None:
        st.error("No match found. Check the spelling, or paste a valid SMILES string such as CC(=O)OC1=CC=CC=C1C(=O)O.")
    return info


def search_box(key):
    c1, c2 = st.columns([5, 1])
    q = c1.text_input("Chemical name or SMILES", value=st.session_state["query"], key=f"in_{key}",
                      placeholder="Aspirin, Caffeine, or CC(=O)OC1=CC=CC=C1C(=O)O")
    c2.write("")
    c2.write("")
    if c2.button("Search", key=f"go_{key}", use_container_width=True):
        st.session_state["query"] = q
        st.rerun()
    cols = st.columns(len(utils.PRESETS))
    for col, name in zip(cols, utils.PRESETS):
        if col.button(name, key=f"{key}_{name}", use_container_width=True):
            st.session_state["query"] = name
            st.rerun()


# ---------- pages ----------
if page == PAGES[0]:
    st.markdown(f"""
    <div class="hero">
      <span class="inst">{INST}</span>
      <h1>Read molecules the way chemists write them.</h1>
      <p>Type a drug name or a SMILES string and get its formula, drug-likeness properties, an interactive 3D model
      and a plain-language summary of benefits and risks.</p>
      <div class="smiles">CC(=O)OC1=CC=CC=C1C(=O)O<br>CN1C=NC2=C1C(=O)N(C(=O)N2C)C</div>
    </div>""", unsafe_allow_html=True)
    st.write("")
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Built-in molecules", len(utils.PRESETS))
    m2.metric("Properties per molecule", 10)
    m3.metric("3D render styles", 4)
    m4.metric("Live database", "PubChem")
    st.subheader("Search a molecule")
    search_box("home")
    info = utils.resolve(st.session_state["query"])
    if info:
        c1, c2 = st.columns([1, 1])
        c1.image(utils.mol_image(info["smiles"]), caption=info["name"])
        c2.markdown(f"**{info['name']}**")
        c2.code(info["smiles"], language=None)
        c2.caption("Open the Molecule & Drug Explorer for full details.")
    st.subheader("What is SMILES?")
    st.markdown("""
SMILES (Simplified Molecular Input Line Entry System) writes a molecule as one line of text. Atoms are letters,
bonds are symbols, and the string can be turned back into a structure by software.
""")
    st.dataframe(pd.DataFrame([
        ("Atoms", "C, N, O", "Aliphatic carbon, nitrogen, oxygen"),
        ("Bonds", "= and #", "Double and triple bond; single bonds are implied"),
        ("Branches", "CC(=O)O", "Parentheses mark a side group"),
        ("Rings", "C1CCCCC1", "Matching digits close a ring (cyclohexane)"),
        ("Aromatic", "c1ccccc1", "Lowercase letters mark aromatic atoms (benzene)"),
    ], columns=["Feature", "Example", "Meaning"]), hide_index=True, use_container_width=True)

elif page == PAGES[1]:
    st.title("Molecule & Drug Explorer")
    search_box("exp")
    info = current()
    if info:
        props = utils.properties(info["smiles"])
        st.subheader(info["name"])
        c1, c2 = st.columns([1, 1.4])
        c1.image(utils.mol_image(info["smiles"]))
        c2.markdown(f"**IUPAC name:** {info['iupac']}")
        c2.markdown(f"**Source:** {info['source']}")
        c2.markdown("**SMILES**")
        c2.code(info["smiles"], language=None)
        k = st.columns(5)
        for col, key in zip(k, ["Molecular Formula", "Molecular Weight (g/mol)", "LogP", "TPSA (Å²)", "Lipinski Violations"]):
            col.metric(key, props[key])
        k = st.columns(2)
        k[0].metric("H-Bond Donors", props["H-Bond Donors"])
        k[1].metric("H-Bond Acceptors", props["H-Bond Acceptors"])
        t1, t2 = st.columns(2)
        t1.markdown("##### Structural properties")
        t1.dataframe(pd.DataFrame(props.items(), columns=["Property", "Value"]).astype(str), hide_index=True, use_container_width=True)
        t2.markdown("##### Elemental composition (including hydrogens)")
        t2.dataframe(pd.DataFrame(utils.element_table(info["smiles"])), hide_index=True, use_container_width=True)
        ins = utils.INSIGHTS.get(info["name"])
        st.markdown("##### Pharmacology")
        if ins:
            st.markdown(f"**Uses:** {ins['use']}")
            st.markdown(f"**Mechanism of action:** {ins['moa']}")
            p, c = st.columns(2)
            p.markdown('<div class="glass proscons-pro"><b>Pros</b><ul>' + "".join(f"<li>{x}</li>" for x in ins["pros"]) + "</ul></div>", unsafe_allow_html=True)
            c.markdown('<div class="glass proscons-con"><b>Cons</b><ul>' + "".join(f"<li>{x}</li>" for x in ins["cons"]) + "</ul></div>", unsafe_allow_html=True)
        else:
            st.info("Curated usage and pros/cons are available for the six built-in molecules. Properties above are computed for any valid SMILES.")

elif page == PAGES[2]:
    st.title("3D Structure Viewer")
    search_box("3d")
    info = current()
    if info:
        a, b = st.columns([1, 1])
        style = a.selectbox("Style", ["Ball and stick", "Stick", "Space filling", "Wireframe"])
        spin = b.toggle("Auto-rotate", value=False)
        sdf = utils.sdf_block(info["smiles"])
        if sdf is None:
            st.error("A 3D conformer could not be generated for this molecule. Try a simpler structure.")
        else:
            st.caption(f"{info['name']}: drag to rotate, scroll to zoom, right-drag to pan.")
            components.html(utils.viewer_html(sdf, style, spin), height=500)
            st.download_button("Download SDF", sdf, file_name=f"{info['name']}.sdf")

elif page == PAGES[3]:
    st.title("Comparative Analysis & Pros/Cons")
    picks = st.multiselect("Choose molecules to compare", list(utils.PRESETS), default=["Aspirin", "Ibuprofen", "Paracetamol"])
    if len(picks) < 2:
        st.info("Pick at least two molecules to compare.")
    else:
        rows = {n: utils.properties(utils.PRESETS[n]) for n in picks}
        df = pd.DataFrame(rows).T
        st.dataframe(df.astype(str), use_container_width=True)
        num = df[["Molecular Weight (g/mol)", "LogP", "TPSA (Å²)"]].astype(float)
        st.bar_chart(num["Molecular Weight (g/mol)"])
        st.bar_chart(num[["LogP", "TPSA (Å²)"]])
        for n in picks:
            ins = utils.INSIGHTS[n]
            with st.expander(f"{n}: benefits and drawbacks", expanded=False):
                st.markdown(f"**Uses:** {ins['use']}")
                st.markdown(f"**Mechanism:** {ins['moa']}")
                p, c = st.columns(2)
                p.success("**Pros**\n\n" + "\n".join(f"- {x}" for x in ins["pros"]))
                c.error("**Cons**\n\n" + "\n".join(f"- {x}" for x in ins["cons"]))
    st.caption("Educational summaries only. They are not medical advice.")

else:
    st.title("About & Credits")
    st.markdown(f"""
**SMILES in Chemistry** is an educational dashboard created at **{INST}**.

- **Chemistry engine:** RDKit (descriptors, 2D drawing, 3D embedding)
- **Data:** PubChem PUG REST API for name lookup and IUPAC names
- **3D viewer:** py3Dmol
- **Interface:** Streamlit

**Developers:** add your names, department and guide here.

Drug information is a simplified teaching summary and must not be used for medical decisions.
""")

st.markdown(f'<div class="footer">Hosted by <b>{INST}</b> · SMILES in Chemistry</div>', unsafe_allow_html=True)
