import streamlit as st
import streamlit.components.v1 as components
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from utils import MOLECULES_DB, GLOBAL_DISEASES_DB, get_molecule_data, generate_3d_sdf

# Page Configuration
st.set_page_config(
    page_title="SMILES in Chemistry | SRM IST Tiruchirappalli",
    page_icon="🧪",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Load Custom CSS Styling
try:
    with open("assets/style.css") as f:
        st.markdown(f"<style>{f.read()}</style>", unsafe_allow_html=True)
except FileNotFoundError:
    pass

# Sidebar Navigation & Branding
st.sidebar.markdown("""
<div style="text-align: center; padding: 10px;">
    <h2 style="color: #2563eb; margin-bottom: 0px;">🧬 SMILES Engine</h2>
    <p style="font-size: 0.85rem; color: #64748b;"><b>SRM Institute of Science and Technology</b><br>Tiruchirappalli, Tamil Nadu</p>
</div>
<hr style="margin: 10px 0;">
""", unsafe_allow_html=True)

nav_selection = st.sidebar.radio(
    "Navigation Hub",
    [
        "🏠 Home / Dashboard",
        "🧪 Molecule & Drug Explorer",
        "🧊 3D Structure & py3Dmol Viewer",
        "📊 Comparative Analysis & Pros/Cons",
        "🌍 Global Health & Disease Solutions",
        "ℹ️ About & Credits"
    ]
)

st.sidebar.markdown("---")
st.sidebar.markdown("### ⚡ Quick Molecule Load")
quick_load = st.sidebar.selectbox("Select preset compound:", ["-- Choose --"] + list(MOLECULES_DB.keys()))

# Main State Initialization
if "active_mol" not in st.session_state:
    st.session_state.active_mol = "Aspirin"

if quick_load != "-- Choose --":
    st.session_state.active_mol = quick_load

# ----------------------------------------------------
# PAGE 1: HOME / DASHBOARD
# ----------------------------------------------------
if nav_selection == "🏠 Home / Dashboard":
    st.markdown("""
    <div style="background: linear-gradient(135deg, #1e3a8a, #2563eb); padding: 40px; border-radius: 16px; color: white; box-shadow: 0 10px 25px rgba(37,99,235,0.2);">
        <h1 style="color: white; margin-bottom: 10px;">SMILES in Chemistry & Pharmacological Intelligence</h1>
        <h3 style="color: #93c5fd; font-weight: 400; margin-top: 0;">SRM Institute of Science and Technology (SRM IST), Tiruchirappalli</h3>
        <p style="font-size: 1.1rem; max-width: 800px; line-height: 1.6; margin-top: 20px;">
            Welcome to the next-generation computational chemistry platform. Explore chemical SMILES structures, interact with real-time 3D molecular conformations, analyze pharmaceutical properties across age brackets, and review global disease solutions and chemical alternatives.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("<br>", unsafe_allow_html=True)
    
    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown("""
        <div class="metric-card">
            <h3>🔬 Presets</h3>
            <p style="font-size: 1.5rem; font-weight: bold; color: #2563eb;">6+ Core Drugs</p>
            <p style="font-size: 0.85rem; color: #64748b;">Expandable with PubChem REST API & RDKit parsing.</p>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.markdown("""
        <div class="metric-card">
            <h3>🧊 3D Engine</h3>
            <p style="font-size: 1.5rem; font-weight: bold; color: #2563eb;">Interactive py3Dmol</p>
            <p style="font-size: 0.85rem; color: #64748b;">Real-time ball-and-stick, wireframe, and surface rendering.</p>
        </div>
        """, unsafe_allow_html=True)
    with col3:
        st.markdown("""
        <div class="metric-card">
            <h3>💊 Medical Intelligence</h3>
            <p style="font-size: 1.5rem; font-weight: bold; color: #2563eb;">Global Disease Index</p>
            <p style="font-size: 0.85rem; color: #64748b;">Pathology, standard treatments, and chemical alternatives.</p>
        </div>
        """, unsafe_allow_html=True)
    with col4:
        st.markdown("""
        <div class="metric-card">
            <h3>🏛️ Institution</h3>
            <p style="font-size: 1.5rem; font-weight: bold; color: #2563eb;">SRM IST Trichy</p>
            <p style="font-size: 0.85rem; color: #64748b;">Department of Computer Science & Chemistry Research.</p>
        </div>
        """, unsafe_allow_html=True)

    st.markdown("---")
    st.subheader("📌 What is SMILES Notation?")
    st.markdown("""
    **SMILES** (Simplified Molecular-Input Line-Entry System) is a chemical notation using short ASCII character strings to specify the structure of chemical molecules. Developed in the late 1980s, it allows computers to store, parse, and compare chemical topologies effortlessly.
    * Example: **Aspirin** is represented as `CC(=O)Oc1ccccc1C(=O)O`
    * Example: **Caffeine** is represented as `CN1C=NC2=C1C(=O)N(C(=O)N2C)C`
    """)

# ----------------------------------------------------
# PAGE 2: MOLECULE & DRUG EXPLORER
# ----------------------------------------------------
elif nav_selection == "🧪 Molecule & Drug Explorer":
    st.header("🧪 Advanced Molecule & Drug Explorer")
    st.write("Search any chemical by common name or paste a raw SMILES string to compute full biochemical parameters.")
    
    user_query = st.text_input("Enter Molecule Name or SMILES string:", value=st.session_state.active_mol)
    
    mol_info = get_molecule_data(user_query)
    
    if mol_info:
        st.success(f"Successfully loaded data for: **{mol_info['name']}** (Source: {mol_info['source']})")
        
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Molecular Formula", mol_info["formula"])
        c2.metric("Molecular Weight", f"{mol_info['mw']} g/mol")
        c3.metric("Lipophilicity (LogP)", mol_info["logp"])
        c4.metric("Topological Polar Surface Area", f"{mol_info['tpsa']} Å²")
        
        st.markdown("---")
        
        tab_prop, tab_med, tab_dos = st.tabs(["📐 Structural Properties", "🩺 Pharmacological Usage & Pros/Cons", "👶 Age-Level Dosages"])
        
        with tab_prop:
            st.subheader("Chemical & Topological Parameters")
            col_a, col_b = st.columns(2)
            with col_a:
                st.write(f"**IUPAC Name:** `{mol_info['iupac']}`")
                st.write(f"**SMILES String:** `{mol_info['smiles']}`")
                st.write(f"**Category:** {mol_info['category']}")
            with col_b:
                st.write(f"**Hydrogen Bond Donors:** {mol_info['h_donors']}")
                st.write(f"**Hydrogen Bond Acceptors:** {mol_info['h_acceptors']}")
                lipinski_pass = "✅ Pass" if mol_info['mw'] <= 500 and mol_info['logp'] <= 5 and mol_info['h_donors'] <= 5 and mol_info['h_acceptors'] <= 10 else "⚠️ Review"
                st.write(f"**Lipinski's Rule of Five:** {lipinski_pass}")
                
        with tab_med:
            st.subheader("Medical & Therapeutic Profile")
            st.info(f"**Primary Uses:** {mol_info['uses']}")
            st.warning(f"**Mechanism of Action:** {mol_info['mechanism']}")
            
            p_col, c_col = st.columns(2)
            with p_col:
                st.markdown("### ✅ Benefits & Pros")
                for pro in mol_info['pros']:
                    st.markdown(f"- {pro}")
            with c_col:
                st.markdown("### ❌ Drawbacks & Cons")
                for con in mol_info['cons']:
                    st.markdown(f"- {con}")
                    
        with tab_dos:
            st.subheader("Age-Level & Demographic Dosage Guidelines")
            for group, dose in mol_info['dosages'].items():
                st.markdown(f"**{group}** \n> {dose}")
            st.caption(f"🌍 **Global Context:** {mol_info['global_context']}")
            
    else:
        st.error("Compound not recognized in local database or PubChem API. Please check your spelling or verify the SMILES string syntax.")

# ----------------------------------------------------
# PAGE 3: 3D STRUCTURE VIEWER
# ----------------------------------------------------
elif nav_selection == "🧊 3D Structure & py3Dmol Viewer":
    st.header("🧊 Interactive 3D Molecular Conformation Viewer")
    st.write("Visualize real-time 3D ball-and-stick models generated via RDKit UFF optimization and rendered with py3Dmol.")
    
    target_mol = st.selectbox("Select Compound for 3D Rendering:", list(MOLECULES_DB.keys()))
    smiles_val = MOLECULES_DB[target_mol]["smiles"]
    
    col_v1, col_v2 = st.columns([1, 2])
    with col_v1:
        st.markdown(f"### {target_mol}")
        st.write(f"**Formula:** {MOLECULES_DB[target_mol]['formula']}")
        st.write(f"**SMILES:** `{smiles_val}`")
        view_style = st.selectbox("Rendering Style", ["ball and stick", "stick", "sphere", "line"])
        bg_color = st.color_picker("Background Color", "#0f172a")
        
    with col_v2:
        sdf_data = generate_3d_sdf(smiles_val)
        if sdf_data:
            # Render py3Dmol inside Streamlit component
            html_code = f"""
            <div id="container" style="width: 100%; height: 450px; position: relative; border-radius: 12px; overflow: hidden; box-shadow: 0 4px 12px rgba(0,0,0,0.2);"></div>
            <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.3/3Dmol-min.js"></script>
            <script>
              let element = document.getElementById("container");
              let config = {{ backgroundColor: "{bg_color}" }};
              let viewer = $3Dmol.createViewer(element, config);
              let vdata = `{sdf_data}`;
              viewer.addModel(vdata, "sdf");
              viewer.setStyle({{}}, {{"{view_style}": {{"radius": 0.25}}}});
              viewer.zoomTo();
              viewer.render();
            </script>
            """
            components.html(html_code, height=470)
        else:
            st.error("Failed to generate 3D spatial coordinates for this molecule.")

# ----------------------------------------------------
# PAGE 4: COMPARATIVE ANALYSIS & PROS/CONS
# ----------------------------------------------------
elif nav_selection == "📊 Comparative Analysis & Pros/Cons":
    st.header("📊 Multi-Molecule Comparative Analysis")
    st.write("Select multiple compounds to compare physicochemical metrics and pharmacological profiles side by side.")
    
    selected_compounds = st.multiselect("Choose compounds to compare:", list(MOLECULES_DB.keys()), default=["Aspirin", "Paracetamol", "Ibuprofen"])
    
    if len(selected_compounds) > 0:
        comp_data = []
        for c in selected_compounds:
            d = MOLECULES_DB[c]
            comp_data.append({
                "Compound": c,
                "Molecular Weight": d["mw"],
                "LogP (Lipophilicity)": d["logp"],
                "TPSA (Å²)": d["tpsa"],
                "H-Donors": d["h_donors"],
                "H-Acceptors": d["h_acceptors"]
            })
        df_comp = pd.DataFrame(comp_data)
        
        st.dataframe(df_comp, use_container_width=True)
        
        st.markdown("<br>", unsafe_allow_html=True)
        col_ch1, col_ch2 = st.columns(2)
        
        with col_ch1:
            st.subheader("Molecular Weight Comparison")
            fig_mw = px.bar(df_comp, x="Compound", y="Molecular Weight", color="Compound", title="MW g/mol Distribution")
            st.plotly_chart(fig_mw, use_container_width=True)
            
        with col_ch2:
            st.subheader("Lipophilicity (LogP) Profile")
            fig_logp = px.scatter(df_comp, x="Compound", y="LogP (Lipophilicity)", size="Molecular Weight", color="Compound", title="LogP vs Compound")
            st.plotly_chart(fig_logp, use_container_width=True)
            
        st.markdown("---")
        st.subheader("Detailed Benefit & Drawback Breakdowns")
        for c in selected_compounds:
            with st.expander(f"📌 {c}: Benefits, Drawbacks & Clinical Uses"):
                d = MOLECULES_DB[c]
                st.write(f"**Uses:** {d['uses']}")
                st.write(f"**Mechanism:** {d['mechanism']}")
                bc1, bc2 = st.columns(2)
                with bc1:
                    st.markdown("**Pros:**")
                    for p in d['pros']: st.markdown(f"- {p}")
                with bc2:
                    st.markdown("**Cons:**")
                    for cn in d['cons']: st.markdown(f"- {cn}")
    else:
        st.info("Please select at least one compound above to view comparative analytics.")

# ----------------------------------------------------
# PAGE 5: GLOBAL HEALTH & DISEASE SOLUTIONS
# ----------------------------------------------------
elif nav_selection == "🌍 Global Health & Disease Solutions":
    st.header("🌍 Global Health Intelligence & Medical Solutions Engine")
    st.write("Explore known world health conditions, pathological burdens, conventional frontline medicines, and chemical alternatives engineered for global welfare.")
    
    disease_choice = st.selectbox("Select Global Health Condition:", list(GLOBAL_DISEASES_DB.keys()))
    dis_info = GLOBAL_DISEASES_DB[disease_choice]
    
    st.markdown(f"""
    <div style="background: rgba(37,99,235,0.05); border-left: 5px solid #2563eb; padding: 20px; border-radius: 8px; margin-bottom: 20px;">
        <h3 style="margin-top:0; color: #1e3a8a;">{disease_choice}</h3>
        <p><b>Global Burden:</b> {dis_info['global_burden']}</p>
        <p><b>Key Pathology:</b> {dis_info['key_pathology']}</p>
    </div>
    """, unsafe_allow_html=True)
    
    dc1, dc2 = st.columns(2)
    with dc1:
        st.subheader("💊 Standard Frontline Medicines")
        for med in dis_info['standard_medicines']:
            st.markdown(f"- **{med}**")
            
    with dc2:
        st.subheader("🌿 Chemical & Biological Alternatives")
        st.write(dis_info['chemical_alternatives'])
        
    st.markdown("---")
    st.info(f"🛡️ **Strategic Prevention & Management Protocol:** {dis_info['preventative_strategy']}")

# ----------------------------------------------------
# PAGE 6: ABOUT & CREDITS
# ----------------------------------------------------
elif nav_selection == "ℹ️ About & Credits":
    st.header("ℹ️ About the Project & Institution")
    st.markdown("""
    ### 🏛️ SRM Institute of Science and Technology (SRM IST), Tiruchirappalli
    This high-performance computational web application was engineered to bridge computer science, molecular informatics, and clinical chemistry. Developed as an interactive research-grade dashboard to demonstrate real-time SMILES parsing, 3D molecular modeling, and global health solution tracking.
    
    * **Institution:** SRM Institute of Science and Technology (SRM IST), Tiruchirappalli, Tamil Nadu
    * **Core Technologies:** Python, Streamlit, RDKit, py3Dmol, Plotly, PubChem REST API.
    * **Academic Disclaimer:** *Educational, computational, and research summaries only. Not intended to serve as professional medical advice or clinical prescription.*
    """)

# Footer Branding
st.markdown("""
<div class="custom-footer">
    Hosted by <b>SRM Institute of Science and Technology (SRM IST), Tiruchirappalli</b> · SMILES in Chemistry Platform
</div>
""", unsafe_allow_html=True)
