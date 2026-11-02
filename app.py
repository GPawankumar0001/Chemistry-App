"""
SMILES in Chemistry: futuristic cheminformatics lab dashboard
SRM Institute of Science and Technology (SRM IST), Tiruchirappalli
Department of Computer Science & Chemistry Research

Run:  streamlit run app.py
"""
from __future__ import annotations

import json
from html import escape as esc
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

import utils

APP_DIR = Path(__file__).parent
INSTITUTION = "SRM Institute of Science and Technology (SRM IST)"
CITY = "Tiruchirappalli"
DEPARTMENT = "Department of Computer Science & Chemistry Research"

st.set_page_config(
    page_title="SMILES in Chemistry | SRM IST Tiruchirappalli",
    page_icon="🧪", layout="wide", initial_sidebar_state="expanded",
)

# --------------------------------------------------------------------------- #
# Styling + small HTML helpers
# --------------------------------------------------------------------------- #
def load_css() -> None:
    css_path = APP_DIR / "assets" / "style.css"
    if css_path.exists():
        st.markdown(f"<style>{css_path.read_text(encoding='utf-8')}</style>", unsafe_allow_html=True)
    else:
        st.warning("assets/style.css not found. Running with default theme.")


def H(html: str) -> str:
    """Collapse HTML to one line so Markdown never treats indentation as a code block."""
    return "".join(line.strip() for line in html.splitlines())


def show(html: str) -> None:
    st.markdown(H(html), unsafe_allow_html=True)


def badge(text: str, kind: str = "info") -> str:
    return f'<span class="badge {kind}">{esc(text)}</span>'


def link(text: str, url: str) -> str:
    return f'<a href="{esc(url, quote=True)}" target="_blank" rel="noopener noreferrer">{esc(text)}</a>'


def card(title: str, body: str, accent: str = "", icon: str = "") -> str:
    cls = "violet" if accent == "violet" else ""
    return f'<div class="glass-card {cls}"><div class="card-title {cls}">{icon} {esc(title)}</div>{body}</div>'


def bullets(items) -> str:
    return "<ul>" + "".join(f"<li>{esc(i)}</li>" for i in items) + "</ul>"


def tiles(items) -> str:
    out = ""
    for label, value, sub in items:
        out += (f'<div class="tile"><div class="tile-label">{esc(label)}</div>'
                f'<div class="tile-value">{esc(str(value))}</div><div class="tile-sub">{esc(sub)}</div></div>')
    return f'<div class="tile-grid">{out}</div>'


def html_table(headers, rows) -> str:
    """Cells are trusted HTML. Escape dynamic values with esc() before passing."""
    head = "".join(f"<th>{esc(h)}</th>" for h in headers)
    body = "".join("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>" for r in rows)
    return f'<div class="table-wrap"><table class="glass-table"><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>'


def disclaimer(text: str = "") -> str:
    base = ("<strong>Educational use only.</strong> This platform is a teaching and research aid, not medical advice. "
            "Do not use it to diagnose, treat or dose any person. Consult a licensed clinician or pharmacist and the official "
            "product label before taking or giving any medicine.")
    return f'<div class="disclaimer">{base} {esc(text)}</div>'


def plot_layout(fig: go.Figure, height: int = 380, **kw) -> go.Figure:
    fig.update_layout(template="plotly_dark", paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(15,23,42,0.35)",
                      font=dict(family="Inter", color="#e2e8f0"), height=height, margin=dict(l=30, r=30, t=50, b=30), **kw)
    return fig


# --------------------------------------------------------------------------- #
# 3D viewer (custom 3Dmol.js with multi-CDN fallback + py3Dmol option)
# --------------------------------------------------------------------------- #
VIEWER_TEMPLATE = """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>
 html,body{margin:0;padding:0;background:transparent;font-family:Inter,system-ui,sans-serif;}
 .wrap{position:relative;width:100%;height:__HEIGHT__px;border-radius:14px;overflow:hidden;background:__BG__;}
 #viewer{position:relative;width:100%;height:100%;}
 .toolbar{position:absolute;right:10px;top:10px;display:flex;gap:6px;z-index:10;}
 .toolbar button{cursor:pointer;color:#e2e8f0;font-size:12px;padding:6px 10px;border-radius:9px;
   border:1px solid rgba(6,182,212,.6);background:rgba(15,23,42,.7);backdrop-filter:blur(6px);}
 .toolbar button:hover{box-shadow:0 0 12px rgba(6,182,212,.7);}
 #msg{position:absolute;left:0;right:0;top:46%;text-align:center;color:#94a3b8;font-size:14px;z-index:5;}
</style></head>
<body>
<div class="wrap">
 <div id="viewer"></div>
 <div class="toolbar">
  <button id="b-spin" title="Toggle rotation">⟳ Spin</button>
  <button id="b-reset" title="Reset view">⌖ Reset</button>
  <button id="b-png" title="Save PNG">⬇ PNG</button>
 </div>
 <div id="msg">Loading 3D engine…</div>
</div>
<script>
const CFG = __PAYLOAD__;
const SOURCES = [
 "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.1/3Dmol-min.js",
 "https://3Dmol.org/build/3Dmol-min.js",
 "https://cdn.jsdelivr.net/npm/3dmol@2.0.1/build/3Dmol-min.js"
];
let viewer = null, spinning = CFG.spin;
function loadFrom(i){
  if (i >= SOURCES.length){
    document.getElementById('msg').innerText = 'Could not load 3Dmol.js (check internet access). Switch engine in the sidebar.';
    return;
  }
  const s = document.createElement('script');
  s.src = SOURCES[i];
  s.onload = () => { if (typeof $3Dmol !== 'undefined') init(); else loadFrom(i+1); };
  s.onerror = () => loadFrom(i+1);
  document.head.appendChild(s);
}
function applyStyle(){
  const cs = {colorscheme:'Jmol'};
  let st = {};
  if (CFG.style === 'Ball and Stick') st = {stick:{radius:0.14,...cs}, sphere:{scale:0.28,...cs}};
  else if (CFG.style === 'Stick') st = {stick:{radius:0.2,...cs}};
  else if (CFG.style === 'Wireframe') st = {line:{linewidth:2,...cs}};
  else st = {sphere:{...cs}};
  viewer.setStyle({}, st);
}
function init(){
  document.getElementById('msg').style.display = 'none';
  const el = document.getElementById('viewer');
  viewer = $3Dmol.createViewer(el, {backgroundColor: CFG.bg, antialias: true});
  viewer.addModel(CFG.sdf, 'sdf');
  applyStyle();
  if (CFG.surface) viewer.addSurface($3Dmol.SurfaceType.VDW, {opacity:0.35, color:'#06b6d4'}, {});
  if (CFG.labels) viewer.addPropertyLabels('elem', {}, {fontSize:12, fontColor:'white', backgroundColor:'#0f172a', backgroundOpacity:0.55, showBackground:true, alignment:'center'});
  viewer.zoomTo(); viewer.render();
  if (spinning) viewer.spin('y', 0.6);
  document.getElementById('b-spin').onclick = () => { spinning = !spinning; viewer.spin(spinning ? 'y' : false); };
  document.getElementById('b-reset').onclick = () => { viewer.zoomTo(); viewer.render(); };
  document.getElementById('b-png').onclick = () => {
    const a = document.createElement('a'); a.href = viewer.pngURI(); a.download = 'molecule.png'; a.click();
  };
  window.addEventListener('resize', () => { viewer.resize(); viewer.render(); });
}
loadFrom(0);
</script>
</body></html>"""


def custom_viewer_html(sdf: str, style: str, bg: str, spin: bool, labels: bool, surface: bool, height: int) -> str:
    payload = json.dumps({"sdf": sdf, "style": style, "bg": bg, "spin": spin, "labels": labels, "surface": surface})
    payload = payload.replace("</", "<\\/")
    return (VIEWER_TEMPLATE.replace("__PAYLOAD__", payload).replace("__HEIGHT__", str(height)).replace("__BG__", bg))


def py3dmol_html(sdf: str, style: str, bg: str, spin: bool, labels: bool, surface: bool, height: int) -> str:
    import py3Dmol
    view = py3Dmol.view(width="100%", height=height)
    view.addModel(sdf, "sdf")
    styles = {
        "Ball and Stick": {"stick": {"radius": 0.14, "colorscheme": "Jmol"}, "sphere": {"scale": 0.28, "colorscheme": "Jmol"}},
        "Stick": {"stick": {"radius": 0.2, "colorscheme": "Jmol"}},
        "Wireframe": {"line": {"linewidth": 2, "colorscheme": "Jmol"}},
        "Space-filling": {"sphere": {"colorscheme": "Jmol"}},
    }
    view.setStyle(styles.get(style, styles["Ball and Stick"]))
    if surface:
        view.addSurface(py3Dmol.VDW, {"opacity": 0.35, "color": "#06b6d4"})
    if labels:
        view.addPropertyLabels("elem", {}, {"fontSize": 12, "fontColor": "white", "backgroundOpacity": 0.5})
    view.setBackgroundColor(bg)
    view.zoomTo()
    if spin:
        view.spin(True)
    return view._make_html()


@st.cache_data(show_spinner=False)
def cached_3d(smiles: str, optimize: bool, n_confs: int) -> dict:
    return utils.generate_3d_sdf(smiles, optimize=optimize, n_confs=n_confs)


@st.cache_data(show_spinner=False)
def cached_svg(smiles: str) -> str:
    return utils.mol_to_svg(smiles)


# --------------------------------------------------------------------------- #
# State
# --------------------------------------------------------------------------- #
MODE_LIB, MODE_NAME, MODE_SMILES = "Curated library", "PubChem name search", "SMILES string"


def init_state() -> None:
    st.session_state.setdefault("input_mode", MODE_LIB)
    st.session_state.setdefault("lib_choice", "aspirin")
    st.session_state.setdefault("name_query", "")
    st.session_state.setdefault("smiles_query", "")
    if "record" not in st.session_state:
        st.session_state["record"] = utils.resolve_molecule("library", "aspirin")


def load_library_molecule(key: str) -> None:
    """Callback used by Global Health buttons."""
    st.session_state["input_mode"] = MODE_LIB
    st.session_state["lib_choice"] = key
    st.session_state["record"] = utils.resolve_molecule("library", key)
    st.session_state["jump_notice"] = f"Loaded {utils.library_name(key)} into the lab."


# --------------------------------------------------------------------------- #
# Sidebar
# --------------------------------------------------------------------------- #
def render_sidebar() -> dict:
    with st.sidebar:
        show(f"""
        <div class="sidebar-brand">
          <div class="logo">🧪</div>
          <div class="name">SMILES in Chemistry</div>
          <div class="inst">{esc(INSTITUTION)}</div>
          <div class="city">{esc(CITY)}</div>
          <div class="inst" style="margin-top:6px">{esc(DEPARTMENT)}</div>
        </div>""")

        st.markdown("### 🔎 Molecule input")
        mode = st.radio("Input method", [MODE_LIB, MODE_NAME, MODE_SMILES], key="input_mode", label_visibility="collapsed")
        analyze = False
        if mode == MODE_LIB:
            st.selectbox("Curated molecule", utils.library_keys(), key="lib_choice", format_func=utils.library_name)
            current = st.session_state.get("record", {}).get("library_key")
            if st.session_state["lib_choice"] != current or st.session_state["record"]["source"] != "Curated library":
                analyze = True
        elif mode == MODE_NAME:
            st.text_input("Compound name", key="name_query", placeholder="e.g., ciprofloxacin, dopamine, caffeine")
            analyze = st.button("⚡ Search PubChem and analyse", use_container_width=True)
        else:
            st.text_area("SMILES", key="smiles_query", height=90, placeholder="e.g., CC(=O)OC1=CC=CC=C1C(=O)O")
            analyze = st.button("⚡ Parse and analyse", use_container_width=True)

        if analyze:
            try:
                with st.spinner("Analysing molecule…"):
                    if mode == MODE_LIB:
                        st.session_state["record"] = utils.resolve_molecule("library", st.session_state["lib_choice"])
                    elif mode == MODE_NAME:
                        st.session_state["record"] = utils.resolve_molecule("name", st.session_state["name_query"])
                    else:
                        st.session_state["record"] = utils.resolve_molecule("smiles", st.session_state["smiles_query"])
            except ValueError as err:
                st.error(str(err))
            except Exception as err:  # defensive: never crash the dashboard
                st.error(f"Unexpected error: {err}")

        st.markdown("---")
        st.markdown("### 🎛️ 3D viewer")
        viewer = {
            "engine": st.radio("Engine", ["Custom 3Dmol.js", "py3Dmol"], horizontal=True, key="v_engine"),
            "style": st.selectbox("Render style", ["Ball and Stick", "Stick", "Wireframe", "Space-filling"], key="v_style"),
            "bg": st.color_picker("Background colour", "#070b19", key="v_bg"),
            "spin": st.toggle("Auto-rotate", value=False, key="v_spin"),
            "labels": st.toggle("Atom labels", value=False, key="v_labels"),
            "surface": st.toggle("VdW surface", value=False, key="v_surface"),
            "height": st.slider("Viewer height (px)", 380, 760, 540, 20, key="v_height"),
            "optimize": st.toggle("UFF geometry optimisation", value=True, key="v_opt"),
            "n_confs": st.slider("Conformers searched", 1, 30, 10, key="v_confs"),
        }
        st.markdown("---")
        st.caption("Data: RDKit · PubChem PUG-REST · curated educational library. Not for clinical use.")
    return viewer


# --------------------------------------------------------------------------- #
# Header / footer
# --------------------------------------------------------------------------- #
def render_header() -> None:
    show(f"""
    <div class="hero">
      <div class="hero-eyebrow">Computational Chemistry · Cheminformatics · Pharmaceutical Intelligence</div>
      <h1 class="hero-title">SMILES in Chemistry</h1>
      <div class="hero-inst">{esc(INSTITUTION)}, {esc(CITY)}</div>
      <div class="hero-sub">{esc(DEPARTMENT)}</div>
      <div class="hero-badges">
        {badge("RDKit engine", "info")}{badge("PubChem live search", "violet")}{badge("3D UFF-optimised", "info")}{badge("ADME · Toxicity · Dosing", "violet")}
      </div>
    </div>""")


def render_footer() -> None:
    show(f"""
    <div class="footer">
      <div class="f-inst">{esc(INSTITUTION)}</div>
      <div class="f-inst" style="font-size:.9rem">{esc(CITY)}, Tamil Nadu, India</div>
      <div class="f-dept">{esc(DEPARTMENT)}</div>
      <div class="f-meta">SMILES in Chemistry · Built with Streamlit, RDKit, 3Dmol.js, Plotly and PubChem PUG-REST.<br>
      For education and research only. Not medical advice. Verify all data against primary sources and official labelling.</div>
    </div>""")


def render_banner(rec: dict) -> None:
    p, lip = rec["props"], rec["lipinski"]
    chips = [badge(rec["source"], "violet"), badge(p["formula"], "info"),
             badge(lip["label"], "pass" if lip["compliant"] else "fail")]
    if rec.get("cid"):
        chips.append(f'<span class="badge info">{link("PubChem CID " + str(rec["cid"]), utils.PUBCHEM_COMPOUND_PAGE + str(rec["cid"]))}</span>')
    if rec.get("curated"):
        chips.append(badge(rec["curated"]["drug_class"].split(";")[0], "warn"))
    show(f'<div class="glass-card"><div class="card-title">🧬 Active molecule</div>'
         f'<h2 style="margin:0 0 10px">{esc(rec["name"])}</h2><div class="hero-badges" style="margin:0">{"".join(chips)}</div></div>')
    if st.session_state.get("jump_notice"):
        st.success(st.session_state.pop("jump_notice"))


# --------------------------------------------------------------------------- #
# Charts
# --------------------------------------------------------------------------- #
def radar_fig(p: dict) -> go.Figure:
    axes = [("MW / 500", p["mw"] / 500), ("LogP / 5", p["logp"] / 5), ("HBD / 5", p["hbd"] / 5),
            ("HBA / 10", p["hba"] / 10), ("TPSA / 140", p["tpsa"] / 140), ("RotB / 10", p["rotb"] / 10)]
    names = [a[0] for a in axes]
    vals = [min(max(a[1], 0), 1.6) for a in axes]
    fig = go.Figure()
    fig.add_trace(go.Scatterpolar(r=[1] * 7, theta=names + names[:1], name="Rule limit",
                                  line=dict(color="#f59e0b", dash="dash"), mode="lines"))
    fig.add_trace(go.Scatterpolar(r=vals + vals[:1], theta=names + names[:1], fill="toself", name="Molecule",
                                  line=dict(color="#06b6d4", width=3), fillcolor="rgba(139,92,246,0.30)"))
    fig.update_layout(polar=dict(bgcolor="rgba(15,23,42,0.35)", radialaxis=dict(range=[0, 1.6], gridcolor="rgba(148,163,184,.25)"),
                                 angularaxis=dict(gridcolor="rgba(148,163,184,.25)")),
                      title="Drug-likeness radar (1.0 = Lipinski limit)", legend=dict(orientation="h", y=-0.1))
    return plot_layout(fig, 400)


def element_fig(rows: list) -> go.Figure:
    palette = ["#06b6d4", "#8b5cf6", "#ec4899", "#10b981", "#f59e0b", "#3b82f6", "#ef4444", "#a3e635"]
    fig = go.Figure(go.Pie(labels=[r["Symbol"] for r in rows], values=[r["Mass %"] for r in rows], hole=0.58,
                           marker=dict(colors=palette, line=dict(color="#0f172a", width=2)), textinfo="label+percent"))
    fig.update_layout(title="Elemental mass composition (%)", showlegend=False)
    return plot_layout(fig, 360)


def pk_fig(sim: dict, name: str) -> go.Figure:
    fig = go.Figure(go.Scatter(x=sim["t"], y=sim["c"], mode="lines", fill="tozeroy",
                               line=dict(color="#06b6d4", width=3), fillcolor="rgba(139,92,246,0.25)"))
    fig.update_layout(title=f"Illustrative plasma concentration: {name}", xaxis_title="Time (h)", yaxis_title="Concentration (mg/L)")
    return plot_layout(fig, 380)


# --------------------------------------------------------------------------- #
# Tabs
# --------------------------------------------------------------------------- #
def tab_overview(rec: dict) -> None:
    p, cur = rec["props"], rec.get("curated")
    left, right = st.columns([1, 1.15], gap="large")
    with left:
        try:
            show(card("2D structure", f'<div class="mol-svg">{cached_svg(rec["canonical"])}</div>', icon="⬡"))
        except Exception as err:
            st.error(f"Could not draw structure: {err}")
        st.code(rec["canonical"], language="text")
    with right:
        show(tiles([("Mol. weight", f'{p["mw"]:.2f}', "g/mol"), ("Exact mass", f'{p["exact_mass"]:.4f}', "monoisotopic"),
                    ("LogP", f'{p["logp"]:.2f}', "Crippen"), ("TPSA", f'{p["tpsa"]:.1f}', "square angstrom"),
                    ("H-bond donors", p["hbd"], "Lipinski"), ("H-bond acceptors", p["hba"], "Lipinski"),
                    ("Rotatable bonds", p["rotb"], "flexibility"), ("QED", f'{p["qed"]:.2f}', "0 to 1")]))
        ident = [("Molecular formula", p["formula"]), ("IUPAC name", rec.get("iupac") or "Not retrieved (use PubChem search)"),
                 ("InChIKey", rec.get("inchikey") or "Unavailable"), ("Heavy atoms", p["heavy_atoms"]),
                 ("Stereocentres", p["stereocenters"]), ("Formal charge", p["charge"])]
        kv = "".join(f'<div class="k">{esc(k)}</div><div class="v">{esc(str(v))}</div>' for k, v in ident)
        show(card("Identifiers", f'<div class="kv">{kv}</div>', "violet", "🆔"))

    if cur:
        show(card("Summary", f'<p>{esc(cur["summary"])}</p><div class="note">Drug class: {esc(cur["drug_class"])}</div>', icon="📖"))
        c1, c2, c3 = st.columns(3)
        c1.markdown(H(card("Indications", bullets(cur["indications"]), icon="🎯")), unsafe_allow_html=True)
        c2.markdown(H(card("Advantages", bullets(cur["pros"]), icon="✅")), unsafe_allow_html=True)
        c3.markdown(H(card("Limitations", bullets(cur["cons"]), "violet", "⚠️")), unsafe_allow_html=True)
    elif rec.get("cid"):
        desc = utils.fetch_pubchem_description(int(rec["cid"]))
        if desc:
            show(card("PubChem description", f'<p>{esc(desc["text"])}</p><div class="note">Source: {esc(desc["source"])}</div>', icon="📖"))
    else:
        st.info("This molecule is outside the curated library. Computed descriptors and heuristics are shown in the other tabs.")

    report = pd.DataFrame([
        ("Name", rec["name"]), ("Canonical SMILES", rec["canonical"]), ("Formula", p["formula"]),
        ("MW", round(p["mw"], 3)), ("Exact mass", round(p["exact_mass"], 5)), ("LogP", round(p["logp"], 3)),
        ("TPSA", round(p["tpsa"], 2)), ("HBD", p["hbd"]), ("HBA", p["hba"]), ("RotB", p["rotb"]),
        ("QED", round(p["qed"], 3)), ("Lipinski violations", rec["lipinski"]["violations"]),
    ], columns=["Property", "Value"])
    st.download_button("⬇ Download property report (CSV)", report.to_csv(index=False).encode("utf-8"),
                       file_name=f"{rec['name'].split()[0].lower()}_report.csv", mime="text/csv")


def tab_3d(rec: dict, v: dict) -> None:
    try:
        with st.spinner("Embedding 3D coordinates (ETKDG) and optimising (UFF)…"):
            data = cached_3d(rec["canonical"], v["optimize"], v["n_confs"])
    except (ValueError, RuntimeError) as err:
        st.error(f"3D generation failed: {err}")
        return
    c = st.columns(4)
    c[0].metric("Atoms (with H)", data["n_atoms"])
    c[1].metric("Conformers", data["n_conformers"])
    c[2].metric("UFF energy", f'{data["best_energy"]:.2f} kcal/mol' if data["best_energy"] is not None else "n/a")
    c[3].metric("Converged", {True: "Yes", False: "Max iters", None: "n/a"}[data["converged"]])

    html = None
    if v["engine"] == "py3Dmol":
        try:
            html = py3dmol_html(data["sdf"], v["style"], v["bg"], v["spin"], v["labels"], v["surface"], v["height"])
        except Exception as err:
            st.warning(f"py3Dmol engine failed ({err}); using custom 3Dmol.js renderer.")
    if html is None:
        html = custom_viewer_html(data["sdf"], v["style"], v["bg"], v["spin"], v["labels"], v["surface"], v["height"])
    st.markdown('<div class="viewer-frame">', unsafe_allow_html=True)
    components.html(html, height=v["height"] + 8, scrolling=False)
    st.markdown("</div>", unsafe_allow_html=True)
    st.caption("Drag to rotate · scroll to zoom · right-drag or Shift+drag to pan. The geometry is the lowest-energy UFF conformer found.")
    st.download_button("⬇ Download 3D structure (SDF)", data["sdf"].encode("utf-8"),
                       file_name=f"{rec['name'].split()[0].lower()}_3d.sdf", mime="chemical/x-mdl-sdfile")


def tab_properties(rec: dict) -> None:
    p, lip = rec["props"], rec["lipinski"]
    show(f'<div class="hero-badges" style="margin-bottom:12px">{badge(lip["label"], "pass" if lip["compliant"] else "fail")}'
         f'{badge(str(lip["violations"]) + " rule violation(s)", "warn" if lip["violations"] else "pass")}</div>')
    c1, c2 = st.columns([1.1, 1], gap="large")
    with c1:
        rows = [(esc(n), esc(f"{v:.2f}" if isinstance(v, float) else str(v)), esc(t),
                 badge("Pass" if ok else "Fail", "pass" if ok else "fail")) for n, v, t, ok in lip["rules"]]
        show(card("Lipinski's Rule of Five", html_table(["Parameter", "Value", "Threshold", "Status"], rows), icon="📏"))
        frows = [(esc(f["name"]), esc(f["criteria"]), badge("Pass" if f["pass"] else "Fail", "pass" if f["pass"] else "fail")) for f in rec["filters"]]
        show(card("Additional drug-likeness filters", html_table(["Filter", "Criteria", "Status"], frows), "violet", "🧮"))
    with c2:
        st.plotly_chart(radar_fig(p), use_container_width=True)

    more = [("Molar refractivity", f'{p["mr"]:.2f}'), ("Fraction Csp3", f'{p["fsp3"]:.2f}'), ("Rings / aromatic rings", f'{p["rings"]} / {p["arom_rings"]}'),
            ("Total atoms (with H)", p["total_atoms"]), ("ESOL logS", f'{p["esol_logs"]:.2f} ({utils.esol_class(p["esol_logs"])})')]
    show(card("Additional descriptors", '<div class="kv">' + "".join(f'<div class="k">{esc(k)}</div><div class="v">{esc(str(v))}</div>' for k, v in more) + "</div>", icon="🔬"))

    e1, e2 = st.columns([1, 1.1], gap="large")
    with e1:
        st.plotly_chart(element_fig(rec["elements"]), use_container_width=True)
    with e2:
        st.dataframe(pd.DataFrame(rec["elements"]), use_container_width=True, hide_index=True)


def tab_adme(rec: dict) -> None:
    cur = rec.get("curated")
    if cur:
        letters = {"Absorption": "A", "Distribution": "D", "Metabolism": "M", "Excretion": "E"}
        cells = "".join(f'<div class="adme-card"><div class="adme-letter">{letters[k]}</div><h4>{k}</h4><p>{esc(v)}</p></div>' for k, v in cur["adme"].items())
        show(f'<div class="adme-grid">{cells}</div>')
        show(card("Pharmacokinetic parameters", html_table(["Parameter", "Typical value"], [(esc(k), esc(v)) for k, v in cur["pk_table"].items()]), icon="⏱️"))
        pk = cur["pk"]
        st.markdown("#### 📈 Illustrative concentration-time simulation")
        s1, s2 = st.columns(2)
        dose = s1.number_input("Single oral dose (mg)", 1.0, 5000.0, float(pk["dose_mg"]), step=1.0, key=f"pk_dose_{rec['library_key']}")
        weight = s2.number_input("Body weight (kg)", 3.0, 150.0, 70.0, step=1.0, key=f"pk_wt_{rec['library_key']}")
        sim = utils.simulate_pk(dose, pk["f"], pk["vd_l_kg"], pk["t_half_h"], pk["tmax_h"], weight)
        st.plotly_chart(pk_fig(sim, rec["name"]), use_container_width=True)
        st.caption(f'Model Cmax about {sim["cmax"]:.3g} mg/L at about {sim["tmax_sim"]:.1f} h. {pk["note"]} '
                   "One-compartment, first-order model with literature-approximate parameters. For teaching only; not for dosing decisions.")
    else:
        st.info("Experimental ADME text is available for the curated library. Below are computed heuristics for this molecule.")
    rows = [(esc(k), esc(v)) for k, v in rec["heuristic_adme"].items()]
    show(card("Computed ADME heuristics", html_table(["Endpoint", "Estimate"], rows) +
              '<div class="note">Rule-of-thumb estimates from descriptors (Egan, Veber, ESOL). They are not experimental measurements.</div>', "violet", "🧪"))


def tab_toxicity(rec: dict) -> None:
    cur = rec.get("curated")
    show(disclaimer("Animal LD50 values do not translate directly to human risk."))
    if cur:
        ld = cur["tox"]["ld50"]
        if ld:
            rows = []
            for r in ld:
                cat = utils.ghs_acute_category(r["value"]) if r["route"] == "Oral" else "n/a"
                rows.append((esc(r["species"]), esc(r["route"]), esc(f'{r["value"]:g} {r["unit"]}'), esc(cat), link(r["source"], r["url"])))
            show(card("Acute toxicity (LD50)", html_table(["Species", "Route", "LD50", "GHS acute oral (derived)", "Source"], rows) +
                      '<div class="note">GHS category derived from the numeric value for orientation only. Vendor SDS values vary by salt form.</div>', icon="☠️"))
        else:
            st.info("No numeric LD50 is curated for this molecule. See the linked sources below.")
        rows = [(esc(m["marker"]), esc(m["detail"]), link(m["source"], m["url"])) for m in cur["tox"]["clinical"]]
        show(card("Clinical safety markers", html_table(["Marker", "Detail", "Source"], rows), "violet", "🩺"))
    alerts = rec["alerts"]
    arows = [(esc(k), esc(", ".join(v)) if v else badge("No alerts", "pass")) for k, v in alerts.items()]
    show(card("Structural alerts (RDKit FilterCatalog)", html_table(["Catalogue", "Matches"], arows) +
              '<div class="note">PAINS = assay interference, Brenk = unwanted/reactive/toxic fragments, NIH = reactive/toxic substructures. Alerts flag review, not proof of toxicity.</div>', icon="🚨"))
    if rec.get("cid"):
        cid = rec["cid"]
        show(card("Primary toxicity sources", f'<ul><li>{link("PubChem toxicity section", utils.PUBCHEM_COMPOUND_PAGE + str(cid) + "#section=Toxicity")}</li>'
                  f'<li>{link("PubChem safety and hazards", utils.PUBCHEM_COMPOUND_PAGE + str(cid) + "#section=Safety-and-Hazards")}</li></ul>', "violet", "🔗"))


def tab_dosing(rec: dict) -> None:
    cur = rec.get("curated")
    show(disclaimer("Doses shown summarise typical label-style ranges for teaching. Individual dosing depends on indication, weight, kidney/liver function, interactions and local guidelines."))
    if not cur:
        st.info("Demographic dosing tables are only provided for curated library drugs. Select one in the sidebar or consult official labelling.")
        return
    icons = {"Adult": "🧑", "Pediatric": "🧒", "Geriatric": "🧓"}
    rows = [(f'<strong>{icons[g]} {g}</strong>', esc(d["dose"]), esc(d["limit"]), esc(d["safety"])) for g, d in cur["dosing"].items()]
    show(card(f'Demographic dosing matrix: {cur["name"]}', html_table(["Population", "Typical dose", "Limit / ceiling", "Safety considerations"], rows), icon="📊"))
    src = (f'<ul><li>{link("Official labels on DailyMed", utils.DAILYMED_SEARCH + cur["dailymed"])}</li>'
           f'<li>{link("DrugBank", utils.DRUGBANK_PAGE + cur["drugbank"])}</li></ul>')
    show(card("Verify against primary sources", src, "violet", "🔗"))


def tab_global() -> None:
    show(disclaimer("Population-level statistics are rounded from WHO and peer-reviewed sources linked below."))
    key = st.radio("Global health challenge", list(utils.GLOBAL_HEALTH.keys()), horizontal=True,
                   format_func=lambda k: f'{utils.GLOBAL_HEALTH[k]["icon"]} {utils.GLOBAL_HEALTH[k]["title"]}')
    g = utils.GLOBAL_HEALTH[key]
    show(card(f'{g["title"]}: global burden', f'<p>{esc(g["burden"])}</p>', icon=g["icon"]))
    rows = [(f'<strong>{esc(t["cls"])}</strong>', esc(t["examples"]), esc(t["mech"]), esc(t["note"])) for t in g["treatments"]]
    show(card("Standard pharmaceutical treatments", html_table(["Class", "Examples", "Mechanism", "Notes"], rows), icon="💊"))
    c1, c2 = st.columns(2)
    c1.markdown(H(card("Chemical and therapeutic alternatives", bullets(g["alternatives"]), "violet", "🧬")), unsafe_allow_html=True)
    c2.markdown(H(card("Preventative strategies", bullets(g["prevention"]), icon="🛡️")), unsafe_allow_html=True)
    show(card("Primary sources", "<ul>" + "".join(f"<li>{link(t, u)}</li>" for t, u in g["links"]) + "</ul>", "violet", "🔗"))
    if g["library"]:
        st.markdown("**Explore related molecules in the lab:**")
        cols = st.columns(len(g["library"]))
        for col, lk in zip(cols, g["library"]):
            col.button(f"🔬 Open {utils.library_name(lk)}", key=f"gh_{key}_{lk}", on_click=load_library_molecule, args=(lk,), use_container_width=True)


def tab_guide() -> None:
    rows = [("C, N, O, S, P, F, Cl, Br, I", "Organic-subset atoms; aromatic atoms are lowercase (c, n, o)", "<code>CCO</code> ethanol"),
            ("= and #", "Double and triple bonds", "<code>C=O</code>, <code>C#N</code>"),
            ("( )", "Branches", "<code>CC(=O)O</code> acetic acid"),
            ("Ring digits", "Ring closure labels", "<code>C1CCCCC1</code> cyclohexane"),
            ("[ ]", "Explicit atoms, charge, isotopes, H count", "<code>[NH4+]</code>"),
            ("@ and @@", "Tetrahedral chirality", "<code>N[C@@H](C)C(=O)O</code> L-alanine"),
            ("/ and \\", "Double-bond geometry (E/Z)", "<code>F/C=C/F</code>")]
    show(card("SMILES quick reference", html_table(["Symbol", "Meaning", "Example"], rows), icon="📚"))
    show(card("Try these in the SMILES input", '<ul><li><code>CC(=O)OC1=CC=CC=C1C(=O)O</code> aspirin</li><li><code>CN1C=NC2=C1C(=O)N(C(=O)N2C)C</code> caffeine</li>'
              '<li><code>C1=CC=CC=C1</code> benzene</li><li><code>CCN(CC)CC</code> triethylamine</li></ul>', "violet", "🧪"))


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main() -> None:
    load_css()
    init_state()
    viewer = render_sidebar()
    render_header()
    rec = st.session_state["record"]
    render_banner(rec)

    tabs = st.tabs(["🧬 Overview", "🔮 3D Lab", "⚗️ Properties", "💊 ADME / PK", "☠️ Toxicity",
                    "📊 Dosing Matrix", "🌍 Global Health", "📚 SMILES Guide"])
    with tabs[0]: tab_overview(rec)
    with tabs[1]: tab_3d(rec, viewer)
    with tabs[2]: tab_properties(rec)
    with tabs[3]: tab_adme(rec)
    with tabs[4]: tab_toxicity(rec)
    with tabs[5]: tab_dosing(rec)
    with tabs[6]: tab_global()
    with tabs[7]: tab_guide()
    render_footer()


if __name__ == "__main__":
    main()
