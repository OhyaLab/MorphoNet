"""MorphoNet Streamlit app with Sigma.js (WebGL) network viewer.

How to run (locally):
  streamlit run PATH/App.py
"""

from __future__ import annotations
import json
from io import BytesIO
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
import streamlit as st
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from html import escape as html_escape

                               
                                   
                               
import os
import sys

def _candidate_base_dirs() -> list[Path]:
    """Locations to check for app resources.

    PyInstaller 6 onedir builds usually place bundled files inside
    dist/MorphoNet/_internal, while users may prefer to keep data next to
    MorphoNet.exe in dist/MorphoNet/data. We support both.
    """
    candidates: list[Path] = []

    for env_name in ("MORPHONET_EXE_DIR", "MORPHONET_BUNDLED_DIR"):
        v = os.environ.get(env_name)
        if v:
            candidates.append(Path(v).resolve())

    if getattr(sys, "frozen", False):
        candidates.append(Path(sys.executable).resolve().parent)

    meipass = getattr(sys, "_MEIPASS", None)
    if meipass:
        candidates.append(Path(meipass).resolve())

    here = Path(__file__).resolve().parent
    candidates.extend([here, here / "_internal"])

                                               
    unique: list[Path] = []
    seen: set[str] = set()
    for c in candidates:
        key = str(c)
        if key not in seen:
            unique.append(c)
            seen.add(key)
    return unique

def resource_path(*parts: str) -> Path:
    """Return an existing resource path, checking external and bundled locations."""
    for base in _candidate_base_dirs():
        candidate = base.joinpath(*parts)
        if candidate.exists():
            return candidate
                                                                                    
    return _candidate_base_dirs()[0].joinpath(*parts)

APP_DIR = resource_path()
DATA_DIR = resource_path("data")
ASSETS_DIR = resource_path("assets")

def data_file(*parts: str) -> Path:
    path = resource_path("data", *parts)
    if not path.exists():
        searched = [str(base / "data" / Path(*parts)) for base in _candidate_base_dirs()]
        raise FileNotFoundError("Missing data file. Checked:\n" + "\n".join(searched))
    return path

def load_text_asset(filename: str) -> str:
    path = resource_path("assets", filename)
    if not path.exists():
        searched = [str(base / "assets" / filename) for base in _candidate_base_dirs()]
        raise FileNotFoundError(
            "Missing offline JavaScript asset. Checked:\n" + "\n".join(searched)
        )
    return path.read_text(encoding="utf-8")

                               
              
                               

                               
                                            
                               
NONESS_FUNCTIONAL_GROUPS: list[str] = [
    "Actin cytoskeleton organization",
    "Autophagy",
    "Carbohydrate metabolism",
    "Carboxylic acid biosynthesis",
    "Cell cycle regulation",
    "Cell differentiation",
    "Cell wall",
    "Ion homeostasis",
    "Mitosis",
    "mRNA splicing",
    "Nuclear transport",
    "Nucleotide metabolism",
    "Protein targeting",
    "Ribosome structure",
    "RNP complex biogenesis",
    "Transcriptional regulation",
    "tRNA processing",
    "Ubiquitin-mediated proteolysis",
]

ESS_FUNCTIONAL_GROUPS: list[str] = [
    "90S preribosome complex",
    "Chromosome condensation",
    "DNA replication initiation",
    "Homologous recombination repair",
    "Nucleobase metabolism",
    "Nucleocytoplasmic transport",
    "Phospholipid metabolism",
    "Ribosomal large subunit biogenesis",
    "Ribosomal small subunit biogenesis",
    "RNA polymerase activity",
    "RNA polymerase II initiation",
    "RNA transport",
    "snRNP complex",
    "Translational nucleic acid binding",
    "Ubiquitin-mediated proteolysis",
]

@st.cache_data(show_spinner=False)
def load_data(base_dir: Path):
                      
    info_df = pd.read_csv(data_file('GeneInfo_Both.csv'), index_col=0)
    info_df = info_df.rename(columns={"SGD.ID": "SGD ID", "Gene.name": "Standard name", "Name_description": "Name description"})

                                                 
                                                                             
                                                                             
                                   
    nodes_both = pd.read_csv(data_file('Network_Nodes_Both.csv'), index_col=0)
    if "Study" not in nodes_both.columns:
        raise ValueError("Network_Nodes_Both.csv must contain a 'Study' column with Essential/NonEssential values.")

    nodes_both["Study"] = nodes_both["Study"].astype(str).str.strip()
    nodes_df_Essen = nodes_both[nodes_both["Study"].str.casefold() == "essential"].copy()
    nodes_df_NonEssen = nodes_both[nodes_both["Study"].str.casefold() == "nonessential"].copy()

                                                     
                                                                             
                                                                            
                                 
    edges_both = pd.read_csv(data_file('Network_Edges_Both.csv'))
    if "Study" not in edges_both.columns:
        raise ValueError("Network_Edges_Both.csv must contain a 'Study' column with Essential/NonEssential values.")

    required_edge_cols = {"x0", "y0", "x1", "y1"}
    missing_edge_cols = required_edge_cols.difference(edges_both.columns)
    if missing_edge_cols:
        raise ValueError(
            "Network_Edges_Both.csv is missing required coordinate columns: "
            + ", ".join(sorted(missing_edge_cols))
        )

    edges_both["Study"] = edges_both["Study"].astype(str).str.strip()
    edges_df_Essen = edges_both[edges_both["Study"].str.casefold() == "essential"].copy()
    edges_df_NonEssen = edges_both[edges_both["Study"].str.casefold() == "nonessential"].copy()

                  
    z_values_df= pd.read_csv(data_file('ZValues_Both.csv'), index_col=0)
    q_values_df= pd.read_csv(data_file('qValues_Both.csv'), index_col=0)

                                                                           
    parameter_desc_df = pd.read_csv(data_file('ParametersDescription.csv'))

                          
    Corr_df_Essen= pd.read_csv(data_file('Essential_CorMat.csv'), index_col=0)
    Corr_df_NonEssen= pd.read_csv(data_file('NonEssential_CorMat.csv'), index_col=0)

    return (
        info_df,
        nodes_df_Essen,
        nodes_df_NonEssen,
        edges_df_Essen,
        edges_df_NonEssen,
        z_values_df,
        q_values_df,
        parameter_desc_df,
        Corr_df_Essen,
        Corr_df_NonEssen,
    )

                               
         
                               
def map_to_systematic_name(input_genes: List[str], info_df: pd.DataFrame) -> Tuple[List[str], List[str]]:
    """Map input identifiers (Systematic / SGD ID / Standard name) to Systematic name (index)."""

    id_mapping = pd.Series(info_df.index.values, index=info_df.index).to_dict()
    sgd_mapping = pd.Series(info_df.index.values, index=info_df["SGD ID"]).to_dict()
    std_name_mapping = pd.Series(info_df.index.values, index=info_df["Standard name"]).to_dict()
    combined = {**id_mapping, **sgd_mapping, **std_name_mapping}

    mapped = [combined[g] for g in input_genes if g in combined]
    unique = list(dict.fromkeys(mapped))                  
    duplicates = [g for g in mapped if mapped.count(g) > 1]
    return unique, sorted(set(duplicates))


def extract_correlations(input_genes: List[str], corr_df: pd.DataFrame) -> Dict[Tuple[str, str], float]:
    """Extract all pairwise correlations among input genes.

    Correlation matrices store values only in the *upper triangle* (and keep
    diagonal/lower triangle as NaN). Depending on the (g1,g2) ordering, a lookup
    may hit NaN even though the correlation exists in the opposite direction.
    This helper therefore tries both (g1,g2) and (g2,g1) and keeps the first
    non-NaN value.
    """
    out: Dict[Tuple[str, str], float] = {}
    for i, g1 in enumerate(input_genes):
        for g2 in input_genes[i + 1 :]:
            r = np.nan
            if g1 in corr_df.index and g2 in corr_df.columns:
                r = corr_df.loc[g1, g2]
            if np.isnan(r) and (g2 in corr_df.index and g1 in corr_df.columns):
                r = corr_df.loc[g2, g1]
            if not np.isnan(r):
                out[(g1, g2)] = float(r)
    return out


def save_svg(fig) -> BytesIO:
    """Save a matplotlib figure as SVG into an in-memory buffer.

    Illustrator can fail to display SVG <text> elements if the font is missing or
    if the SVG uses font features it can’t interpret. To make the SVG robust, we
    embed all text as vector paths (svg.fonttype="path").
    """
    buf = BytesIO()
                                                                                          
    old_fonttype = mpl.rcParams.get("svg.fonttype", "path")
    mpl.rcParams["svg.fonttype"] = "path"
    try:
        fig.savefig(buf, format="svg", bbox_inches="tight")
    finally:
        mpl.rcParams["svg.fonttype"] = old_fonttype
    buf.seek(0)
    return buf

def save_svg_editable(fig) -> BytesIO:
    """Save a matplotlib figure as SVG with *editable* text.
    This keeps text as SVG <text> elements (svg.fonttype="none"), which is what
    Adobe Illustrator needs for text to remain editable.
    """
    buf = BytesIO()
    old_fonttype = mpl.rcParams.get("svg.fonttype", "none")
    old_family = mpl.rcParams.get("font.family", None)
    mpl.rcParams["svg.fonttype"] = "none"
                                                                                
    mpl.rcParams["font.family"] = ["Arial", "DejaVu Sans", "sans-serif"]
    try:
        fig.savefig(buf, format="svg", bbox_inches="tight")
    finally:
        mpl.rcParams["svg.fonttype"] = old_fonttype
        if old_family is not None:
            mpl.rcParams["font.family"] = old_family
    buf.seek(0)
    return buf


def plot_sigma_graph_for_export(graph: Dict) -> plt.Figure:
    """Render the Sigma graph to a matplotlib figure for Illustrator-friendly export.

    - Uses the node/edge coordinates already in the Sigma graph dict.
    - Draws labels only for highlighted (red) nodes and in-graph group labels.
    """
    nodes = [n for n in graph.get("nodes", []) if not n.get("hidden", False)]
    edges = [e for e in graph.get("edges", []) if not e.get("hidden", False)]

    fig, ax = plt.subplots(figsize=(10, 10), dpi=200)
    ax.set_aspect("equal", adjustable="box")
    ax.axis("off")

    if not nodes:
        return fig

    xs = np.array([float(n.get("x", 0.0)) for n in nodes], dtype=float)
    ys = np.array([float(n.get("y", 0.0)) for n in nodes], dtype=float)
    xmin, xmax = float(xs.min()), float(xs.max())
    ymin, ymax = float(ys.min()), float(ys.max())
    span = max(xmax - xmin, ymax - ymin, 1e-9)
    pad = 0.10 * span
    ax.set_xlim(xmin - pad, xmax + pad)
    ax.set_ylim(ymin - pad, ymax + pad)

                             
    if edges:
                                           
        pos = {n["id"]: (float(n.get("x", 0.0)), float(n.get("y", 0.0))) for n in nodes}
        segs = []
        for e in edges:
            s = e.get("source")
            t = e.get("target")
            if s in pos and t in pos:
                segs.append([pos[s], pos[t]])
        if segs:
            lc = mpl.collections.LineCollection(
                segs,
                linewidths=max(0.2, 0.0015 * span * 72 / 10),                              
                colors=["#D3D3D3"],
                alpha=0.85,
            )
            ax.add_collection(lc)

           
    colors = [n.get("color", "#000000") for n in nodes]
                                                                      
                                                                      
    base_sizes = np.array([float(n.get("size", 3.0)) for n in nodes], dtype=float)
    ms = (base_sizes * (220.0 / max(span, 1e-9))) ** 2
    ax.scatter(xs, ys, s=ms, c=colors, linewidths=0)

                                    
    label_fs = max(6.0, 0.020 * span * 100)             
    dy = 0.030 * span
    for n in nodes:
        c = str(n.get("color", "")).upper()
        if c != "#FF3B30":
            continue
        txt = n.get("label", n.get("id", ""))
        ax.text(
            float(n.get("x", 0.0)),
            float(n.get("y", 0.0)) + dy,
            str(txt),
            ha="center",
            va="bottom",
            fontsize=label_fs,
            color="#FF3B30",
            fontfamily="Arial",
        )

                                             
    grp_labels = graph.get("group_labels", []) or []
    grp_fs = max(7.0, 0.022 * span * 100)
    grp_dy = 0.035 * span
    for gl in grp_labels:
        try:
            gtxt = str(gl.get("text", "")).strip()
            if not gtxt:
                continue
            ax.text(
                float(gl.get("x", 0.0)),
                float(gl.get("y", 0.0)) - grp_dy,
                gtxt,
                ha="center",
                va="top",
                fontsize=grp_fs,
                fontweight="bold",
                color=str(gl.get("color", "#111111")),
                fontfamily="Arial",
            )
        except Exception:
            continue

    return fig


def save_tiff(fig, dpi: int = 300) -> BytesIO:
    """Save a matplotlib figure as TIFF into an in-memory buffer."""
    buf = BytesIO()
                                                             
    fig.savefig(buf, format="tiff", dpi=dpi, bbox_inches="tight", facecolor="white")
    buf.seek(0)
                                                               
    return buf


def similarity_network_json_bytes(
    genes: List[str],
    correlations: Dict[Tuple[str, str], float],
    layout: Dict[str, np.ndarray],
    info_df: pd.DataFrame,
) -> bytes:
    """Create a JSON export (nodes + edges) for similarity networks."""
    nodes = []
    for g in genes:
        pos = layout.get(g)
        x = float(pos[0]) if pos is not None else 0.0
        y = float(pos[1]) if pos is not None else 0.0
        std = None
        if g in info_df.index and "Standard name" in info_df.columns:
            v = info_df.loc[g, "Standard name"]
            std = None if pd.isna(v) else str(v)
        nodes.append({"id": g, "label": (std.lower() if std else g), "x": x, "y": y, "sys": g, "std": std})

    edges = []
    for (g1, g2), r in correlations.items():
        edges.append({"id": f"{g1}__{g2}", "source": g1, "target": g2, "r": float(r)})

    payload = {"nodes": nodes, "edges": edges, "metadata": {"type": "morphological_similarity", "n_genes": len(genes)}}
    return json.dumps(payload, indent=2).encode("utf-8")


def extract_correlations_threshold(
    input_genes: List[str],
    corr_df: pd.DataFrame,
    threshold: float,
) -> Dict[Tuple[str, str], float]:
    """Extract correlations among input genes and keep those with abs(r) >= threshold.
    Handles upper-triangle storage by trying both (g1,g2) and (g2,g1).
    """
    out: Dict[Tuple[str, str], float] = {}
    for i, g1 in enumerate(input_genes):
        for g2 in input_genes[i + 1 :]:
            r = np.nan
            if g1 in corr_df.index and g2 in corr_df.columns:
                r = corr_df.loc[g1, g2]
            if np.isnan(r) and (g2 in corr_df.index and g1 in corr_df.columns):
                r = corr_df.loc[g2, g1]
            if not np.isnan(r) and abs(float(r)) >= float(threshold):
                out[(g1, g2)] = float(r)
    return out


def generate_layout(
    genes: List[str],
    correlations: Dict[Tuple[str, str], float],
    seed: int,
):
    """Generate a stable spring layout for the similarity network."""
    G = nx.Graph()
    for g in genes:
        G.add_node(g)
    for (g1, g2), r in correlations.items():
        if abs(r) > 0:
            G.add_edge(g1, g2, weight=float(r))
    return nx.spring_layout(G, seed=seed)


def plot_similarity_network(
    input_genes: List[str],
    correlations: Dict[Tuple[str, str], float],
    info_df: pd.DataFrame,
    layout=None,
    show_r_values: bool = False,
    edge_label_fmt: str = "{:.2f}",
):
    """Matplotlib similarity network: edges colored by r, with colorbar."""
    G = nx.Graph()

    for g in input_genes:
        G.add_node(g)

    for (g1, g2), r in correlations.items():
        G.add_edge(g1, g2, weight=float(r))

    if layout is None:
        layout = nx.spring_layout(G)

    if G.edges:
        edges, weights = zip(*nx.get_edge_attributes(G, "weight").items())
    else:
        edges, weights = ([], [])

    cmap = plt.cm.coolwarm
    norm = mcolors.Normalize(vmin=-1, vmax=1)
    edge_colors = [cmap(norm(w)) for w in weights] if weights else []

    fig, ax = plt.subplots()

    node_size = max(2000 // max(len(G.nodes()), 1), 200)

    if edges:
        nx.draw_networkx_edges(G, layout, edge_color=edge_colors, width=2, ax=ax)

    nx.draw_networkx_nodes(
        G,
        layout,
        node_color="white",
        node_size=node_size,
        edgecolors="black",
        ax=ax,
    )

    if show_r_values and edges:
        edge_labels = {
            e: edge_label_fmt.format(w)
            for e, w in nx.get_edge_attributes(G, "weight").items()
        }
        nx.draw_networkx_edge_labels(
            G,
            layout,
            edge_labels=edge_labels,
            font_size=7,
            rotate=False,
            bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7),
            ax=ax,
        )

    labels = {
        node: (str(info_df.loc[node, "Standard name"]).lower() if node in info_df.index else node)
        for node in G.nodes()
    }

    nx.draw_networkx_labels(G, layout, labels, font_size=8, ax=ax)

    if edges:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Correlation (r)")

    ax.set_axis_off()
    plt.tight_layout()
    return fig, layout

def convert_df_to_file(df: pd.DataFrame, file_type: str = "csv") -> bytes:
    buffer = BytesIO()
    if file_type == "csv":
        df.to_csv(buffer, sep=",", index=True)
    elif file_type == "tsv":
        df.to_csv(buffer, sep="\t", index=True)
    buffer.seek(0)
    return buffer.getvalue()


def _scale_positions(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize x,y roughly into [-1,1] for consistent sigma zoom behavior."""
    out = df.copy()
    for col in ("x", "y"):
        v = out[col].astype(float).values
        if np.nanmax(v) == np.nanmin(v):
            out[col] = 0.0
        else:
            out[col] = 2 * (v - np.nanmin(v)) / (np.nanmax(v) - np.nanmin(v)) - 1
    return out



def _read_uploaded_gene_list(uploaded) -> List[str]:
    """Read a 1-column gene list from Streamlit UploadedFile (CSV/TSV/TXT)."""
    if uploaded is None:
        return []
    try:
        name = getattr(uploaded, "name", "") or ""
        if name.lower().endswith(".csv"):
            df = pd.read_csv(uploaded, header=None)
        elif name.lower().endswith(".tsv"):
            df = pd.read_csv(uploaded, sep="	", header=None)
        else:
                                    
            df = pd.read_csv(uploaded, sep="", header=None, engine="python")
        return df[0].astype(str).str.upper().tolist()
    except Exception as e:
        st.error(f"Error processing uploaded file: {e}")
        return []


def _panel_genes(panel_key: str, info_df: pd.DataFrame) -> Tuple[List[str], List[str], List[str]]:
    """Green-box inputs: return (raw_unique_inputs, mapped_unique_systematic, duplicates_after_mapping)."""
    txt = st.text_input(
        "Search a gene (comma separated)",
       placeholder="e.g., CDC25, FKS1",
    key=f"{panel_key}_txt",
    ).upper().strip()

    up = st.file_uploader(
        "Upload a gene list (CSV/TSV/TXT)",
        key=f"{panel_key}_upload",
    )

    genes = [g.strip() for g in txt.split(",") if g.strip()] if txt else []
    genes += _read_uploaded_gene_list(up)

    genes = list(dict.fromkeys([g for g in genes if str(g).strip()]))
    mapped, dups = map_to_systematic_name(genes, info_df)
    return genes, mapped, dups



def build_sigma_graph(
    nodes_df: pd.DataFrame,
    info_df: pd.DataFrame,
    highlighted: List[str] | None = None,
    correlations: Dict[Tuple[str, str], float] | None = None,
    edges_df: pd.DataFrame | None = None,
    show_edges_precomputed: bool = True,
    visible_groups: set[str] | None = None,
    group_colors: Dict[str, str] | None = None,
    group_col: str = "GO",
    color_col: str = "Col",
    display_name_col: str | None = None,
    color_by_group: bool = False,
    colored_groups: set[str] | None = None,
    default_node_color: str = "#B0B0B0",
    size_scale: float = 1.0,
) -> Dict:
    """Create a JSON-serializable graph for Sigma.js (graphology).

    Styling policy (simple for v1):
      - nodes: white
      - edges: gray
      - highlighted nodes: slightly larger
    """
    highlighted = highlighted or []

                                                                
    nodes_df_raw = nodes_df.copy()
    nodes_df = _scale_positions(nodes_df)

    def _sf(v: float) -> float:
        try:
            v = float(v)
            if np.isfinite(v):
                return v
        except Exception:
            pass
        return 0.0

    nodes = []
    for sys_name, row in nodes_df.iterrows():
        label = sys_name
        std = None
        if sys_name in info_df.index and "Standard name" in info_df.columns:
            std = str(info_df.loc[sys_name, "Standard name"])
            label = std.lower() if std and std != "nan" else sys_name

        is_hi = sys_name in highlighted

                                                           
        grp = None
        if group_col in nodes_df_raw.columns:
            try:
                grp = nodes_df_raw.loc[sys_name, group_col]
            except Exception:
                grp = None
        if grp is None or str(grp) == "nan":
            grp = None
        grp = str(grp) if grp is not None else None

        hidden = False
        if visible_groups is not None and grp is not None and grp not in visible_groups:
            hidden = True
        elif visible_groups is not None and grp is None:
                                                                            
            hidden = True

        base_color = default_node_color
        if (
            color_by_group
            and grp is not None
            and group_colors
            and grp in group_colors
            and (colored_groups is None or grp in colored_groups)
        ):
            base_color = str(group_colors[grp])
        disp_name = None
        if display_name_col and display_name_col in nodes_df_raw.columns:
            try:
                disp_name = nodes_df_raw.loc[sys_name, display_name_col]
            except Exception:
                disp_name = None
        if disp_name is not None and str(disp_name) != "nan" and str(disp_name).strip():
            label = str(disp_name)

        if is_hi:
            base_color = "#FF3B30"

        nodes.append(
            {
                "id": sys_name,
                "label": label,
                "x": _sf(row.get("x", 0.0)),
                "y": _sf(row.get("y", 0.0)),
                "size": (6.0 if is_hi else 3.0) * float(size_scale),
                "color": base_color,
                "baseColor": base_color,
                "labelColor": base_color,
                "baseLabelColor": base_color,
                "baseSize": (6.0 if is_hi else 3.0) * float(size_scale),
                "hidden": bool(hidden),
                "group": grp,
                "sys": sys_name,
                "std": std,
            }
        )

    edges = []
    hidden_nodes = {n['id'] for n in nodes if n.get('hidden')}

                                                                  
    if show_edges_precomputed and edges_df is not None and len(edges_df) > 0:
        cols = {c.lower(): c for c in edges_df.columns}

        def _add_edge(s: str, t: str, key: str | None = None):
            if s is None or t is None:
                return
            s = str(s)
            t = str(t)
            if s == t:
                return
                                                                             
            node_ids = {n["id"] for n in nodes}
            if s not in node_ids or t not in node_ids:
                return
            if s in hidden_nodes or t in hidden_nodes:
                return
            edges.append(
                {
                    "id": key or f"{s}__{t}__pre",
                    "source": s,
                    "target": t,
                    "weight": 1.0,
                    "size": 0.05,
                    "color": "#D3D3D3",
                }
            )

                                        
        if ("source" in cols and "target" in cols) or ("s" in cols and "t" in cols):
            s_col = cols.get("source", cols.get("s"))
            t_col = cols.get("target", cols.get("t"))
            for i, row in edges_df.iterrows():
                _add_edge(row.get(s_col), row.get(t_col), key=str(i))

                                                                                 
        elif all(k in cols for k in ("x0", "y0", "x1", "y1")):
                                                          
            coord_map = {}
            for sys_name, r in nodes_df_raw.iterrows():
                try:
                    coord_map[(round(float(r["x"]), 6), round(float(r["y"]), 6))] = sys_name
                except Exception:
                    continue

            x0c, y0c, x1c, y1c = cols["x0"], cols["y0"], cols["x1"], cols["y1"]
            for i, row in edges_df.iterrows():
                try:
                    s = coord_map.get((round(float(row[x0c]), 6), round(float(row[y0c]), 6)))
                    t = coord_map.get((round(float(row[x1c]), 6), round(float(row[y1c]), 6)))
                    _add_edge(s, t, key=str(i))
                except Exception:
                    continue
                                                                 
    if correlations:
        for (g1, g2), r in correlations.items():
            if g1 in hidden_nodes or g2 in hidden_nodes:
                continue
            edges.append(
                {
                    "id": f"{g1}__{g2}",
                    "source": g1,
                    "target": g2,
                    "weight": float(r),
                    "size": 1.2,
                    "color": "#D3D3D3",
                }
            )


                                                                            
    group_labels = []
    if colored_groups:
        for grp_name in sorted(set(colored_groups)):
            pts = [(n["x"], n["y"]) for n in nodes if (n.get("group") == grp_name and not n.get("hidden", False))]
            if not pts:
                continue

                                                                                      
            xs = [p[0] for p in pts]
            ys = [p[1] for p in pts]
            cx = float(sum(xs) / len(xs))
            y_bottom = float(min(ys))

            col = None
            if group_colors and grp_name in group_colors:
                col = str(group_colors[grp_name])
            if not col or col.lower() == "nan":
                col = "#111111"

            group_labels.append({"text": str(grp_name), "x": cx, "y": y_bottom, "color": col})

    return {"nodes": nodes, "edges": edges, "group_labels": group_labels}



def sigma_viewer_component(
    graph: Dict,
    height: int = 650,
    default_threshold: float = 0.5,
    show_threshold_slider: bool = False,
    title: str = "",
) -> None:
    """Embed a Sigma.js viewer in Streamlit.

    Notes:
      - Avoids Python f-strings inside the HTML/JS template (token replacement only).
      - Adds a simple on-page error console to catch silent JS errors.
      - Exports SVG (full-graph extent) instead of PNG.
    """
    import uuid

    uid = uuid.uuid4().hex[:10]
    container_id = f"sigma_container_{uid}"
    thr_id = f"sigma_thr_{uid}"
    thr_val_id = f"sigma_thrval_{uid}"
    err_id = f"sigma_err_{uid}"

    graph_json = json.dumps(graph)
    sigma_js = load_text_asset("sigma.min.js")
    graphology_js = load_text_asset("graphology.umd.min.js")

    html = r"""
<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8"/>
<script>__SIGMA_JS__</script>
<script>__GRAPHOLOGY_JS__</script>
<style>
  body { margin: 0; background: transparent; }
  .mn-wrap { font-family: Arial, sans-serif; }
  .mn-title { color:#111; font-weight:600; margin: 0 0 8px 0; }
  .controls {
    display: flex;
    gap: 10px;
    align-items: end;
    flex-wrap: wrap;
    margin: 0 0 10px 0;
  }
  .controls input[type="text"]{
    padding: 8px 10px;
    border-radius: 10px;
    border: 1px solid rgba(0,0,0,0.2);
    min-width: 260px;
  }
  .controls button{
    padding: 8px 10px;
    border-radius: 10px;
    border: 1px solid rgba(0,0,0,0.2);
    background: #fff;
    cursor: pointer;
  }
  .controls button:hover{ background:#f4f4f4; }
  .sliderwrap{ display: flex; gap: 8px; align-items:center; }
  .sliderwrap input[type="range"]{ width: 220px; }
  #__CONTAINER_ID__ { width: 100%; height: __HEIGHT__px; background:#FFFFFF; border-radius: 18px; border: 1px solid rgba(0,0,0,0.2); }
  .foot { margin-top: 6px; color: #444; font-size: 12px; display:flex; justify-content: space-between; gap:10px; flex-wrap: wrap; }
  .mono { font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, "Liberation Mono", "Courier New", monospace; }
  #__ERR_ID__ { color:#b00020; font-size:12px; margin-top:6px; white-space: pre-wrap; }

.grpLabel{
  position: absolute;
  transform: translate(-50%, 0%);
  font-size: 12px;
  font-weight: 700;
  color: #111;
  text-shadow: 0 1px 0 rgba(255,255,255,0.85);
  white-space: nowrap;
}

  /* Export dropdown (pure HTML/CSS, no external deps) */
  .exportMenu { position: relative; display: inline-block; }
  .exportMenu > summary {
    list-style: none;
    padding: 8px 10px;
    border-radius: 10px;
    border: 1px solid rgba(0,0,0,0.2);
    background: #fff;
    cursor: pointer;
    user-select: none;
  }
  .exportMenu > summary::-webkit-details-marker { display: none; }
  .exportMenu[open] > summary { background: #f4f4f4; }
  .exportMenuContent{
    position: absolute;
    top: calc(100% + 6px);
    right: 0;
    min-width: 120px;
    background: #fff;
    border: 1px solid rgba(0,0,0,0.2);
    border-radius: 12px;
    box-shadow: 0 8px 22px rgba(0,0,0,0.12);
    padding: 6px;
    display: flex;
    flex-direction: column;
    gap: 6px;
    z-index: 9999;
  }
  .exportMenuContent button{
    text-align: left;
    width: 100%;
    padding: 8px 10px;
    border-radius: 10px;
    border: 1px solid rgba(0,0,0,0.15);
    background: #fff;
    cursor: pointer;
  }
  .exportMenuContent button:hover{ background:#f4f4f4; }

</style>
</head>
<body>

<div class="mn-wrap">
  __TITLE_BLOCK__

  <div class="controls">
    __SLIDER_BLOCK__

    <div style="display:flex; gap:8px; align-items:end;">
      <button id="btnReset">Reset</button>

      <details class="exportMenu" id="exportMenu">
        <summary>Export ▾</summary>
        <div class="exportMenuContent">
          <button id="btnExportSVG_AI" type="button">SVG</button>
          <button id="btnExportTIFF" type="button">TIFF</button>
          <button id="btnExportJSON" type="button">JSON</button>
        </div>
      </details>
    </div>
  </div>

  <div class="sigmaWrap" style="position:relative;">
  <div id="__CONTAINER_ID__"></div>
  <div id="groupLabels" style="position:absolute; left:0; top:0; width:100%; height:100%; pointer-events:none;"></div>
</div>

  <div class="foot">
    <div id="hoverInfo" class="mono">Hover a node to see details. Click a node to highlight neighbors.</div>
    <div>Scroll = zoom · Drag = pan</div>
  </div>

  <div id="__ERR_ID__"></div>
</div>

<script>
(function(){
  const errBox = document.getElementById("__ERR_ID__");
  function logErr(msg){
    try{
      errBox.textContent = (errBox.textContent ? errBox.textContent + "\n" : "") + msg;
    }catch(e){}
  }
  window.addEventListener("error", (e) => logErr("JS error: " + (e.message || e.error)));
  window.addEventListener("unhandledrejection", (e) => logErr("Promise rejection: " + (e.reason || e)));

  const rawData = __GRAPH_DATA__;

  if (!rawData || !rawData.nodes) {
    logErr("Graph data missing or malformed.");
    return;
  }
  logErr(""); // clear

  const Graph = graphology.Graph;
  const graph = new Graph();

  // Nodes
  rawData.nodes.forEach(n => {
    try{
      const id = String(n.id);
      if (!graph.hasNode(id)) {
        graph.addNode(id, {
          label: (n.label != null ? String(n.label) : id),
          x: Number(n.x),
          y: Number(n.y),
          size: Number(n.size || 3),
          color: (n.color || "#000000"),
          baseColor: (n.baseColor || n.color || "#000000"),
          labelColor: (n.labelColor || "#000000"),
          baseLabelColor: (n.baseLabelColor || n.labelColor || "#000000"),
          baseSize: Number(n.baseSize || n.size || 3),
          sys: n.sys || id,
          std: n.std || null
        });
      }
    } catch (e){
      logErr("Node error: " + e);
    }
  });

  // Edges (optional)
  if (rawData.edges && rawData.edges.length){
    rawData.edges.forEach(e => {
      try{
        const s = String(e.source), t = String(e.target);
        if (!graph.hasNode(s) || !graph.hasNode(t)) return;
        const key = e.id ? String(e.id) : (s + "__" + t);
        if (!graph.hasEdge(key)) {
          graph.addEdgeWithKey(key, s, t, {
            weight: Number(e.weight || 0),
            size: Number(e.size || 1),
            color: e.color || "#D3D3D3",
          });
        }
      } catch (err){
        logErr("Edge error: " + err);
      }
    });
  }

  const container = document.getElementById("__CONTAINER_ID__");
  if (!container) {
    logErr("Container not found.");
    return;
  }

  let renderer;
  try{
    renderer = new Sigma(graph, container, {
      renderEdgeLabels: false,
      zIndex: true,
      labelColor: { attribute: "labelColor", color: "#000000" },
      nodeReducer: (node, data) => {
        // Use per-node labelColor when provided (e.g., highlighted genes in red).
        const res = Object.assign({}, data);
        if (data && data.labelColor) res.labelColor = data.labelColor;
        return res;
      },
      // Make labels readable on a black background.
      labelRenderedSizeThreshold: 999,
      defaultLabelColor: "#000000",
      // Keep edges visible even for large graphs.
      minEdgeThickness: 0.2,
      maxEdgeThickness: 1.2,
    });
  }catch(e){
    logErr("Sigma init failed: " + e);
    return;
  }

  // Show labels only when zooming in.
  // Sigma camera ratio: larger = zoomed out, smaller = zoomed in.
  const camera = renderer.getCamera();
  let _lastThr = null;
  function updateLabelVisibility(){
    const ratio = camera.getState().ratio;
    // Tune thresholds if needed:
    // - zoomed out: hide all labels
    // - medium: show only larger nodes
    // - zoomed in: show all labels
    let thr;
    if (ratio > 1.2) thr = 999;       // hide all
    else if (ratio > 0.7) thr = 6;    // show only larger labels
    else thr = 0;                     // show all labels
    if (_lastThr !== thr){
      _lastThr = thr;
      try { renderer.setSetting("labelRenderedSizeThreshold", thr); } catch(e){}
      try { renderer.refresh(); } catch(e){}
    }
  }
  updateLabelVisibility();

// Group labels (functional group names) rendered as HTML overlay, anchored to graph coordinates
const groupLabelHost = document.getElementById("groupLabels");
const groupLabels = (rawData.group_labels || []);
const labelEls = [];

function clearGroupLabels(){
  if (!groupLabelHost) return;
  while (groupLabelHost.firstChild) groupLabelHost.removeChild(groupLabelHost.firstChild);
  labelEls.length = 0;
}


function initGroupLabels(){
  clearGroupLabels();
  if (!groupLabelHost || !groupLabels || !groupLabels.length) return;
  groupLabels.forEach((gl) => {
    const d = document.createElement("div");
    d.className = "grpLabel";
    d.textContent = String(gl.text || "");
    if (gl.color) {
      d.style.color = String(gl.color);
    }
    groupLabelHost.appendChild(d);
    labelEls.push({ el: d, x: Number(gl.x || 0), y: Number(gl.y || 0) });
  });
}

function updateGroupLabelPositions(){
  if (!renderer || !labelEls.length) return;
  const cam = renderer.getCamera();
  // graphToViewport accounts for camera state
  labelEls.forEach(({el, x, y}) => {
    const p = renderer.graphToViewport({x, y});
    el.style.left = p.x + "px";
    el.style.top  = p.y + "px";
  });
}

initGroupLabels();
updateGroupLabelPositions();
camera.on("updated", () => {
  updateLabelVisibility();
  updateGroupLabelPositions();
});

  const dimNode = "rgba(0,0,0,0.15)";
  const baseNode = "#000000";
  const neighNode = "rgba(0,0,0,0.70)";
  const activeNode = "#FF3B30";

  function resetColors() {
    graph.forEachNode((k) => {
      const baseColor = graph.getNodeAttribute(k, "baseColor") || baseNode;
      const baseLabelColor = graph.getNodeAttribute(k, "baseLabelColor") || "#000000";
      const baseSize = graph.getNodeAttribute(k, "baseSize") || 3;
      graph.setNodeAttribute(k, "color", baseColor);
      graph.setNodeAttribute(k, "labelColor", baseLabelColor);
      graph.setNodeAttribute(k, "size", baseSize);
      graph.setNodeAttribute(k, "hidden", false);
    });
    graph.forEachEdge((k) => graph.setEdgeAttribute(k, "hidden", false));
  }

  function resetView() {
    resetColors();
    renderer.getCamera().animatedReset({ duration: 500 });
    renderer.refresh();
  }

  function highlightNeighbors(nodeKey) {
    const neigh = new Set(graph.neighbors(nodeKey));
    neigh.add(nodeKey);

    graph.forEachNode((k) => {
      if (neigh.has(k)) {
        if (k === nodeKey) {
          graph.setNodeAttribute(k, "color", activeNode);
          graph.setNodeAttribute(k, "labelColor", activeNode);
          graph.setNodeAttribute(k, "size", 7);
        } else {
          graph.setNodeAttribute(k, "color", neighNode);
          graph.setNodeAttribute(k, "labelColor", neighNode);
          graph.setNodeAttribute(k, "size", 4.5);
        }
      } else {
        graph.setNodeAttribute(k, "color", dimNode);
        graph.setNodeAttribute(k, "labelColor", dimNode);
        graph.setNodeAttribute(k, "size", 2);
      }
    });

    graph.forEachEdge((e, attrs, s, t) => {
      graph.setEdgeAttribute(e, "hidden", !(neigh.has(s) && neigh.has(t)));
    });

    renderer.refresh();
  }

  // Hover info
  const hoverEl = document.getElementById("hoverInfo");
  renderer.on("enterNode", (e) => {
    const a = graph.getNodeAttributes(e.node);
    const sys = a.sys || e.node;
    const std = a.std ? (" (" + a.std + ")") : "";
    hoverEl.textContent = sys + std;
  });
  renderer.on("leaveNode", () => {
    hoverEl.textContent = "Hover a node to see details. Click a node to highlight neighbors.";
  });

  renderer.on("clickNode", (e) => highlightNeighbors(e.node));

  // Threshold slider (optional)
  const thrEl = document.getElementById("__THR_ID__");
  const thrValEl = document.getElementById("__THR_VAL_ID__");
  function applyThreshold(t) {
    if (!thrEl) return;
    graph.forEachEdge((k, attrs) => {
      const w = Math.abs(attrs.weight || 0);
      graph.setEdgeAttribute(k, "hidden", w < t);
    });
    renderer.refresh();
  }
  if (thrEl) {
    applyThreshold(parseFloat(thrEl.value));
    thrEl.addEventListener("input", () => {
      const t = parseFloat(thrEl.value);
      if (thrValEl) thrValEl.textContent = t.toFixed(2);
      applyThreshold(t);
    });
  }

  // Export JSON
  document.getElementById("btnExportJSON").addEventListener("click", () => {
    const exportObj = { nodes: [], edges: [] };
    graph.forEachNode((k, attrs) => exportObj.nodes.push({ id: k, ...attrs }));
    graph.forEachEdge((k, attrs, s, t) => exportObj.edges.push({ id: k, source: s, target: t, ...attrs }));
    const blob = new Blob([JSON.stringify(exportObj, null, 2)], { type: "application/json" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "morphonet_graph.json";
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
    closeExportMenu();
  });

  // Export SVG (full graph extent, Illustrator-friendly sizing)
  function escapeXml(s){
    return String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/\'/g,"&apos;");
  }

    function exportSVG(mode){
    // mode:
    //   "editable" -> emits SVG <text> labels (editable in Illustrator when imported as text)
    //   "ai_visible" -> embeds labels as a raster overlay so Illustrator always shows them (NOT editable)
    let xmin=Infinity, xmax=-Infinity, ymin=Infinity, ymax=-Infinity;
    graph.forEachNode((k, a) => {
      const x = a.x, y = a.y;
      if (x<xmin) xmin=x;
      if (x>xmax) xmax=x;
      if (y<ymin) ymin=y;
      if (y>ymax) ymax=y;
    });
    if (!isFinite(xmin)) return;

    const pad = 0.10 * Math.max(xmax-xmin, ymax-ymin, 1e-6);
    xmin -= pad; xmax += pad; ymin -= pad; ymax += pad;

    const vbW = xmax - xmin;
    const vbH = ymax - ymin;

    // Flip y so the exported SVG matches the on-screen orientation
    const yFlip = (y) => (ymin + ymax) - y;

    // --- Size tuning (graph units -> readable in Illustrator) ---
    const scaleUnit = Math.max(vbW, vbH);
    const edgeW   = 0.0015 * scaleUnit;
    const nodeR   = 0.0070 * scaleUnit;
    const hiNodeR = 0.0100 * scaleUnit;
    const labelSz = 0.0150 * scaleUnit;
    const labelDy = 0.0150 * scaleUnit;
    const grpSz   = 0.0280 * scaleUnit;
    const grpDy   = 0.0150 * scaleUnit;

    // Choose an explicit pixel canvas size so Illustrator doesn't import as tiny.
    const outWpx = 1800;
    const outHpx = Math.round(outWpx * (vbH / vbW));

    const lines = [];
    const circles = [];
    const labels = [];
    const grpText = [];

    // edges
    graph.forEachEdge((ek, ea, s, t) => {
      if (ea.hidden) return;
      const sA = graph.getNodeAttributes(s);
      const tA = graph.getNodeAttributes(t);
      if (sA.hidden || tA.hidden) return;
      const x1=sA.x, y1=yFlip(sA.y), x2=tA.x, y2=yFlip(tA.y);
      const col = escapeXml(ea.color || "#111111");
      const w = ea.size ? (edgeW * ea.size) : edgeW;
      lines.push(`<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${col}" stroke-width="${w}" stroke-linecap="round" />`);
    });

    // nodes + (optional) highlighted labels
    graph.forEachNode((k, a) => {
      if (a.hidden) return;
      const x=a.x, y=yFlip(a.y);
      const baseC = (a.baseColor || a.color || "#999999");
      const col = escapeXml(baseC);
      const r = (String(baseC).toUpperCase()==="#FF3B30") ? hiNodeR : nodeR;
      circles.push(`<circle cx="${x}" cy="${y}" r="${r}" fill="${col}" stroke="none"/>`);

      if (mode === "editable"){
        // Only label highlighted (red) nodes
        const isRed = String(baseC).toUpperCase()==="#FF3B30";
        if (isRed){
          const txt = escapeXml(a.label || k);
          labels.push(`<text x="${x}" y="${y + labelDy}" fill="#FF3B30" font-family="Arial, sans-serif" font-size="${labelSz}" font-weight="400" text-anchor="middle">${txt}</text>`);
        }
      }
    });

    // group labels (optional)
    if (mode === "editable"){
      (groupLabels || []).forEach((gl) => {
        try{
          const gx = Number(gl.x || 0);
          const gy = Number(gl.y || 0);
          const raw = String(gl.text || "");
          if (!raw) return;
          const gtxt = escapeXml(raw);
          const gcol = gl.color ? escapeXml(gl.color) : "#111111";
          grpText.push(
            `<text x="${gx}" y="${yFlip(gy - grpDy)}" fill="${gcol}" font-family="Arial, sans-serif" font-size="${grpSz}" font-weight="700" text-anchor="middle">${gtxt}</text>`
          );
        } catch(e){}
      });
    }

    // Optional: label raster overlay for Illustrator visibility (NOT editable)
    let labelsPngHref = "";
    if (mode === "ai_visible"){
      try{
        const hasAnyLabel = true; // overlay includes both highlighted labels + group labels
        if (hasAnyLabel){
          const c = document.createElement("canvas");
          c.width = outWpx; c.height = outHpx;
          const ctx = c.getContext("2d");
          ctx.clearRect(0,0,outWpx,outHpx);

          const sx = outWpx / vbW;
          const sy = outHpx / vbH;
          function drawTextAt(xu, yu, text, color, sizeUser, weight){
            const px = (xu - xmin) * sx;
            const py = (yu - ymin) * sy;
            const fs = Math.max(10, sizeUser * sx);
            ctx.font = `${weight || "400"} ${fs}px Arial, sans-serif`;
            ctx.fillStyle = color || "#111111";
            ctx.textAlign = "center";
            ctx.textBaseline = "middle";
            ctx.fillText(text, px, py);
          }

          // highlighted node labels
          graph.forEachNode((k, a) => {
            if (a.hidden) return;
            const baseC = (a.baseColor || a.color || "").toUpperCase();
            const isRed = (baseC === "#FF3B30");
            if (!isRed) return;
            const txt = String(a.label || k);
            drawTextAt(a.x, yFlip(a.y + labelDy), txt, "#FF3B30", labelSz, "400");
          });

          // group labels
          (groupLabels || []).forEach((gl) => {
            try{
              const gx = Number(gl.x || 0);
              const gy = Number(gl.y || 0);
              const raw = String(gl.text || "");
              if (!raw) return;
              const gcol = gl.color ? String(gl.color) : "#111111";
              drawTextAt(gx, yFlip(gy - grpDy), raw, gcol, grpSz, "700");
            } catch(e){}
          });

          labelsPngHref = c.toDataURL("image/png");
        }
      } catch(e){
        labelsPngHref = "";
      }
    }

    const svg =
`<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="${outWpx}px" height="${outHpx}px" viewBox="${xmin} ${ymin} ${vbW} ${vbH}" preserveAspectRatio="xMidYMid meet">
  <rect x="${xmin}" y="${ymin}" width="${vbW}" height="${vbH}" fill="white"/>
  ${lines.join("\n")}
  ${circles.join("\n")}
  ${labelsPngHref ? `<image href="${labelsPngHref}" xlink:href="${labelsPngHref}" x="${xmin}" y="${ymin}" width="${vbW}" height="${vbH}" preserveAspectRatio="none"/>` : ""}
  ${labels.join("\n")}
  ${grpText.join("\n")}
</svg>`;

    const blob = new Blob([svg], { type: "image/svg+xml;charset=utf-8" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = (mode === "ai_visible") ? "Morphonet_graph.svg" : "morphonet_graph_editable.svg";
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(url);
    closeExportMenu();
  }
  document.getElementById("btnExportSVG_AI").addEventListener("click", () => exportSVG("ai_visible"));



// TIFF encoder (uncompressed RGBA, single-strip, little-endian)
  function encodeTIFF_RGBA(rgba, width, height){
    // rgba: Uint8ClampedArray length width*height*4
    const entries = [
      // tag, type, count, value/offset placeholder
      {tag:256, type:4, count:1, value:width},           // ImageWidth (LONG)
      {tag:257, type:4, count:1, value:height},          // ImageLength (LONG)
      {tag:258, type:3, count:4, value:null},            // BitsPerSample (SHORT[4]) -> offset
      {tag:259, type:3, count:1, value:1},               // Compression = 1 (none)
      {tag:262, type:3, count:1, value:2},               // PhotometricInterpretation = 2 (RGB)
      {tag:273, type:4, count:1, value:null},            // StripOffsets (LONG) -> offset
      {tag:277, type:3, count:1, value:4},               // SamplesPerPixel = 4
      {tag:278, type:4, count:1, value:height},          // RowsPerStrip = height
      {tag:279, type:4, count:1, value:width*height*4},  // StripByteCounts
      {tag:284, type:3, count:1, value:1},               // PlanarConfiguration = 1 (chunky)
      {tag:338, type:3, count:1, value:1},               // ExtraSamples = 1 (alpha)
    ];

    const numEntries = entries.length;
    const ifdOffset = 8;
    const ifdSize = 2 + numEntries*12 + 4;

    const bitsOffset = ifdOffset + ifdSize;  // where BitsPerSample array will live
    const bitsSize = 8; // 4 SHORTs

    const imageOffset = bitsOffset + bitsSize;
    const imageSize = width*height*4;

    const totalSize = imageOffset + imageSize;
    const buf = new ArrayBuffer(totalSize);
    const view = new DataView(buf);
    const u8 = new Uint8Array(buf);

    // Header: "II", 42, IFD offset
    view.setUint8(0, 0x49); view.setUint8(1, 0x49);
    view.setUint16(2, 42, true);
    view.setUint32(4, ifdOffset, true);

    // IFD
    view.setUint16(ifdOffset, numEntries, true);
    let p = ifdOffset + 2;

    function writeEntry(tag, type, count, value){
      view.setUint16(p, tag, true);
      view.setUint16(p+2, type, true);
      view.setUint32(p+4, count, true);
      // type sizes: 3=SHORT(2),4=LONG(4)
      if (type === 3 && count === 1){
        view.setUint16(p+8, value, true);
        view.setUint16(p+10, 0, true);
      } else if (type === 4 && count === 1){
        view.setUint32(p+8, value, true);
      } else {
        // offset
        view.setUint32(p+8, value, true);
      }
      p += 12;
    }

    for (const e of entries){
      if (e.tag === 258) writeEntry(e.tag, e.type, e.count, bitsOffset);
      else if (e.tag === 273) writeEntry(e.tag, e.type, e.count, imageOffset);
      else writeEntry(e.tag, e.type, e.count, e.value);
    }
    view.setUint32(p, 0, true); // next IFD offset

    // BitsPerSample: 8,8,8,8
    view.setUint16(bitsOffset+0, 8, true);
    view.setUint16(bitsOffset+2, 8, true);
    view.setUint16(bitsOffset+4, 8, true);
    view.setUint16(bitsOffset+6, 8, true);

    // Image data (top-left origin). Canvas ImageData is already RGBA row-major.
    u8.set(rgba, imageOffset);

    return new Uint8Array(buf);
  }

  function closeExportMenu(){
    const m = document.getElementById("exportMenu");
    if (m) m.open = false;
  }

  // Close export menu when clicking outside
  document.addEventListener("click", (ev) => {
    const m = document.getElementById("exportMenu");
    if (!m) return;
    const target = ev.target;
    if (m.open && target && !m.contains(target)) m.open = false;
  });

  // Export TIFF (rasterize SVG to avoid WebGL preserveDrawingBuffer issues)
  const btnTiff = document.getElementById("btnExportTIFF");
  if (btnTiff){
    btnTiff.addEventListener("click", async () => {
      try{
        // Reuse the same SVG we export for the SVG option
        let xmin=Infinity, xmax=-Infinity, ymin=Infinity, ymax=-Infinity;
        graph.forEachNode((k, a) => {
          const x = a.x, y = a.y;
          if (x<xmin) xmin=x;
          if (x>xmax) xmax=x;
          if (y<ymin) ymin=y;
          if (y>ymax) ymax=y;
        });
        if (!isFinite(xmin)) { closeExportMenu(); return; }

        const pad = 0.10 * Math.max(xmax-xmin, ymax-ymin, 1e-6);
        xmin -= pad; xmax += pad; ymin -= pad; ymax += pad;

        const vbW = xmax - xmin;
        const vbH = ymax - ymin;

    const yFlip = (y) => (ymin + ymax) - y;

        // Size tuning (match SVG export so TIFF is readable when opened in Illustrator)
        const scaleUnit = Math.max(vbW, vbH);
        const edgeW   = 0.0015 * scaleUnit;
        const nodeR   = 0.0070 * scaleUnit;
        const labelSz = 0.0150 * scaleUnit;
        const labelDy = 0.0300 * scaleUnit;
        const grpSz   = 0.0160 * scaleUnit;
        const grpDy   = 0.0350 * scaleUnit;

        const lines = [];
        graph.forEachEdge((k, attrs, s, t) => {
          if (attrs.hidden) return;
          const aS = graph.getNodeAttributes(s);
          const aT = graph.getNodeAttributes(t);
          lines.push(
            `<line x1="${aS.x}" y1="${yFlip(aS.y)}" x2="${aT.x}" y2="${yFlip(aT.y)}" stroke="#D3D3D3" stroke-width="${edgeW}" opacity="0.8"/>`
          );
        });

        const circles = [];
        graph.forEachNode((k, a) => {
          if (a.hidden) return;
          circles.push(
            `<circle cx="${a.x}" cy="${yFlip(a.y)}" r="${nodeR}" fill="${a.color || "#000000"}" />`
          );
        });

        const labels = [];
        graph.forEachNode((k, a) => {
          if (a.hidden) return;
          const baseC = (a.baseColor || a.color || "").toUpperCase();
          const lblC  = (a.labelColor || a.baseLabelColor || "").toUpperCase();
          const isRed = (baseC === "#FF3B30") || (lblC === "#FF3B30");
          if (!isRed) return;
          const txt = (function escapeXml(s){
            return String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/'/g,"&apos;");
          })(a.label || k);
          labels.push(
            `<text x="${a.x}" y="${yFlip(a.y + labelDy)}" fill="#FF3B30" font-family="Arial" font-size="${labelSz}" text-anchor="middle">${txt}</text>`
          );
        });

const grpText = [];
        (groupLabels || []).forEach((gl) => {
          try{
            const gx = Number(gl.x || 0);
            const gy = Number(gl.y || 0);
            const raw = gl.text || "";
            const gtxt = (function escapeXml2(s){
              return String(s).replace(/&/g,"&amp;").replace(/</g,"&lt;").replace(/>/g,"&gt;").replace(/"/g,"&quot;").replace(/'/g,"&apos;");
            })(raw);
            if (!gtxt) return;
            const gcol = gl.color ? String(gl.color) : "#111111";
            grpText.push(
              `<text x="${gx}" y="${yFlip(gy - grpDy)}" fill="${gcol}" font-family="Arial" font-size="${grpSz}" font-weight="700" text-anchor="middle">${gtxt}</text>`
            );
          } catch(e){}
        });

        // Choose a reasonable raster size
        const MAX_DIM = 2500; // pixels
        const scale = MAX_DIM / Math.max(vbW, vbH);
        const outW = Math.max(1, Math.round(vbW * scale));
        const outH = Math.max(1, Math.round(vbH * scale));

        const svg =
`<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="${outW}" height="${outH}" viewBox="${xmin} ${ymin} ${vbW} ${vbH}">
  <rect x="${xmin}" y="${ymin}" width="${vbW}" height="${vbH}" fill="white"/>
  ${lines.join("\n")}
  ${circles.join("\n")}
  ${labels.join("\n")}
  ${grpText.join("\n")}
</svg>`;

        const svgBlob = new Blob([svg], { type: "image/svg+xml;charset=utf-8" });
        const svgUrl = URL.createObjectURL(svgBlob);

        // Rasterize SVG into a 2D canvas
        const img = new Image();
        const loaded = new Promise((resolve, reject) => {
          img.onload = () => resolve(true);
          img.onerror = (e) => reject(e);
        });
        img.src = svgUrl;
        await loaded;

        const tmp = document.createElement("canvas");
        tmp.width = outW; tmp.height = outH;
        const ctx = tmp.getContext("2d");
        ctx.fillStyle = "#FFFFFF";
        ctx.fillRect(0,0,outW,outH);
        ctx.drawImage(img, 0, 0, outW, outH);

        URL.revokeObjectURL(svgUrl);

        const imageData = ctx.getImageData(0,0,outW,outH);
        const tiffBytes = encodeTIFF_RGBA(imageData.data, outW, outH);
        const blob = new Blob([tiffBytes], {type: "image/tiff"});
        const url = URL.createObjectURL(blob);
        const a = document.createElement("a");
        a.href = url;
        a.download = "morphonet_view.tiff";
        document.body.appendChild(a);
        a.click();
        a.remove();
        URL.revokeObjectURL(url);
      } catch(e){
        logErr("TIFF export failed: " + (e && e.message ? e.message : e));
      }
      closeExportMenu();
    });
  }

  // Reset button
  document.getElementById("btnReset").addEventListener("click", () => resetView());

  // Initial camera reset & refresh (some environments need it)
  setTimeout(() => {
    try { renderer.getCamera().animatedReset({ duration: 0 }); } catch(e){}
    try { renderer.refresh(); } catch(e){}
  }, 50);

})();
</script>
</body>
</html>
"""

    title_block = f'<div class="mn-title">{title}</div>' if title else ""
    if show_threshold_slider:
        slider_block = f"""
<div class="sliderwrap">
  <div style="font-size:12px; color:#333; margin-bottom:4px;">|r| threshold</div>
  <input type="range" id="{thr_id}" min="0" max="1" step="0.01" value="{default_threshold}">
  <span id="{thr_val_id}" class="mono">{default_threshold:.2f}</span>
</div>
"""
    else:
        slider_block = ""

    html = (
        html.replace("__GRAPH_DATA__", graph_json)
        .replace("__CONTAINER_ID__", container_id)
        .replace("__HEIGHT__", str(height))
        .replace("__TITLE_BLOCK__", title_block)
        .replace("__SLIDER_BLOCK__", slider_block)
        .replace("__THR_ID__", thr_id)
        .replace("__THR_VAL_ID__", thr_val_id)
        .replace("__ERR_ID__", err_id)
        .replace("__SIGMA_JS__", sigma_js)
        .replace("__GRAPHOLOGY_JS__", graphology_js)
    )

    st.components.v1.html(html, height=height + 170, scrolling=False)

                               
     
                               
st.set_page_config(page_title="MorphoNet", layout="wide")


def install_tab_close_shutdown_handler() -> None:
    """Ask for confirmation when the browser tab is closed, then notify the launcher to stop Streamlit.

    Chrome/Edge do not allow custom text in the close warning; they show their own generic message.
    """
    shutdown_url = os.environ.get("MORPHONET_SHUTDOWN_URL", "").strip()
    if not shutdown_url:
        return

    html = f"""
<script>
(function() {{
  const shutdownUrl = {json.dumps(shutdown_url)};
  let sent = false;

  function sendShutdown() {{
    if (sent) return;
    sent = true;
    try {{
      navigator.sendBeacon(shutdownUrl, new Blob(["close"], {{type: "text/plain"}}));
    }} catch (e) {{
      try {{ fetch(shutdownUrl, {{method: "POST", keepalive: true, mode: "no-cors"}}); }} catch (_) {{}}
    }}
  }}

  window.addEventListener("beforeunload", function(e) {{
    e.preventDefault();
    e.returnValue = "";  // Required for Chrome/Edge to show the built-in confirmation warning.
  }});

  window.addEventListener("pagehide", sendShutdown);
  window.addEventListener("unload", sendShutdown);
}})();
</script>
"""
    st.components.v1.html(html, height=0, width=0)

install_tab_close_shutdown_handler()

                               
          
                               
def get_app_logo_path() -> Path | None:
    """Return the MorphoNet logo path if it is available.

    Put the logo file at assets/MorphoNet_logo.tiff for the offline app/exe.
    For development runs, a Logo(2).tiff file next to the script is also accepted.
    """
    candidates = [
        resource_path("assets", "MorphoNet_logo.tiff"),
        resource_path("assets", "Logo(2).tiff"),
        resource_path("MorphoNet_logo.tiff"),
        resource_path("Logo(2).tiff"),
    ]
    for p in candidates:
        if p.exists():
            return p
    return None

_logo_path = get_app_logo_path()

_logo_col, _title_col = st.columns([1, 7], vertical_alignment="center")
with _logo_col:
    if _logo_path is not None:
        st.image(str(_logo_path), use_container_width=True)
with _title_col:
    st.title("MorphoNet")

                                                              
try:
    (
        info_df,
        nodes_df_Essen,
        nodes_df_NonEssen,
        edges_df_Essen,
        edges_df_NonEssen,
        z_values_df,
        q_values_df,
        parameter_desc_df,
        Corr_df_Essen,
        Corr_df_NonEssen,
    ) = load_data(APP_DIR)
except Exception as e:
    st.error(f"Could not load data with the current working directory and hard-coded paths.\n\nError: {e}")
    st.stop()

                                                                                                 
                                                            
if "global_genes_raw" not in st.session_state:
    st.session_state["global_genes_raw"] = []
if "global_genes_mapped" not in st.session_state:
    st.session_state["global_genes_mapped"] = []
if "global_genes_dups" not in st.session_state:
    st.session_state["global_genes_dups"] = []


tab1, tab2, tab3, tab4, tab5 = st.tabs([
    "Morphology–Gene Function Networks",
    "Morphological Similarity Networks",
    "Statistical Profiles",
    "Morphological Profile",
    "About",
])



def display_filtered_values(values_df: pd.DataFrame, file_name_prefix: str):
    st.info("Use the search box under Morphology–Gene Function Networks tab.")
    raw = st.session_state.get("global_genes_raw", [])
    mapped = st.session_state.get("global_genes_mapped", [])
    dups = st.session_state.get("global_genes_dups", [])

    st.info(f"Number of input genes: {len(raw)} · Found: {len(mapped)}")

    if dups:
        st.warning("There are duplicate gene names in the input")

    if mapped:
        filtered = values_df.loc[mapped].copy()
        filtered = filtered.merge(info_df[["Standard name"]], left_index=True, right_index=True)
        filtered = filtered[["Standard name"] + [c for c in filtered.columns if c != "Standard name"]]

                                               
        try:
            export_box = st.popover("Export", use_container_width=False)
        except Exception:
                                                   
            export_box = st.expander("Export", expanded=False)

        with export_box:
            st.download_button(
                label="CSV",
                data=convert_df_to_file(filtered, "csv"),
                file_name=f"{file_name_prefix}.csv",
                mime="text/csv",
                key=f"dl_{file_name_prefix}_csv",
                use_container_width=True,
            )
            st.download_button(
                label="TSV",
                data=convert_df_to_file(filtered, "tsv"),
                file_name=f"{file_name_prefix}.tsv",
                mime="text/tab-separated-values",
                key=f"dl_{file_name_prefix}_tsv",
                use_container_width=True,
            )

        st.dataframe(filtered)


def display_values_table_only(values_df: pd.DataFrame, file_name_prefix: str):
    """Render only the filtered table + export buttons (no info/warning blocks).
    """
    mapped = st.session_state.get("global_genes_mapped", [])
    if not mapped:
        return

    filtered = values_df.loc[mapped].copy()
    filtered = filtered.merge(info_df[["Standard name"]], left_index=True, right_index=True)
    filtered = filtered[["Standard name"] + [c for c in filtered.columns if c != "Standard name"]]

                                           
    try:
        export_box = st.popover("Export", use_container_width=False)
    except Exception:
                                               
        export_box = st.expander("Export", expanded=False)

    with export_box:
        st.download_button(
            label="CSV",
            data=convert_df_to_file(filtered, "csv"),
            file_name=f"{file_name_prefix}.csv",
            mime="text/csv",
            key=f"dl_{file_name_prefix}_csv",
            use_container_width=True,
        )
        st.download_button(
            label="TSV",
            data=convert_df_to_file(filtered, "tsv"),
            file_name=f"{file_name_prefix}.tsv",
            mime="text/tab-separated-values",
            key=f"dl_{file_name_prefix}_tsv",
            use_container_width=True,
        )

    st.dataframe(filtered)


                               
                               
                               
def normalize_parameter_description_table(parameter_desc_df: pd.DataFrame) -> pd.DataFrame:
    """Return a two-column table indexed by CalMorph parameter name.

    The app accepts common column names, but also falls back to the first two
    columns so it can work with slightly different ParametersDescription.csv
    formats.
    """
    desc = parameter_desc_df.copy()
    if desc.empty:
        return pd.DataFrame(columns=["Description"])

    col_lookup = {str(c).strip().casefold(): c for c in desc.columns}
    param_candidates = [
        "parameter", "parameters", "parameter name", "parameter_name",
        "calmorph parameter", "calmorph_parameter", "name", "id"
    ]
    desc_candidates = [
        "description", "parameter description", "parameter_description",
        "definition", "detail", "details"
    ]

    param_col = next((col_lookup[c] for c in param_candidates if c in col_lookup), None)
    desc_col = next((col_lookup[c] for c in desc_candidates if c in col_lookup), None)

    if param_col is None:
        param_col = desc.columns[0]
    if desc_col is None:
        desc_col = desc.columns[1] if len(desc.columns) > 1 else desc.columns[0]

    out = desc[[param_col, desc_col]].copy()
    out.columns = ["Parameter", "Description"]
                                                                                
                                                                      
    out["Parameter"] = out["Parameter"].astype(str).str.strip().str.replace(".", "-", regex=False)
    out = out.dropna(subset=["Parameter"])
    out = out[out["Parameter"].astype(str).str.len() > 0]
    out = out.drop_duplicates(subset=["Parameter"], keep="first")
    return out.set_index("Parameter")


def display_calmorph_parameter_name(parameter: str) -> str:
    """Return user-facing CalMorph parameter name with dots replaced by dashes."""
    return str(parameter).replace(".", "-")


def get_gene_display_names(sys_gene: str, info_df: pd.DataFrame) -> tuple[str, str, str]:
    """Return systematic, standard, and SGD ID names for a mapped gene."""
    std = ""
    sgd = ""
    if sys_gene in info_df.index:
        if "Standard name" in info_df.columns and not pd.isna(info_df.loc[sys_gene, "Standard name"]):
            std = str(info_df.loc[sys_gene, "Standard name"])
        if "SGD ID" in info_df.columns and not pd.isna(info_df.loc[sys_gene, "SGD ID"]):
            sgd = str(info_df.loc[sys_gene, "SGD ID"])
    return sys_gene, std, sgd


def get_gene_name_description(sys_gene: str, info_df: pd.DataFrame) -> str:
    """Return the gene name-description field from GeneInfo_Both.csv when available."""
    if sys_gene not in info_df.index:
        return ""

                                                                               
                                                                             
    normalized_cols = {str(c).strip().casefold(): c for c in info_df.columns}
    candidate_names = [
        "name description",
        "name.description",
        "name_description",
        "gene description",
        "gene.description",
        "gene_description",
        "description",
        "desc",
    ]
    desc_col = next((normalized_cols[c] for c in candidate_names if c in normalized_cols), None)
    if desc_col is None:
        return ""

    val = info_df.loc[sys_gene, desc_col]
    if pd.isna(val):
        return ""
    return str(val)


def format_q_value_sci2(value) -> str:
    """Format q values in scientific notation with two decimals."""
    try:
        if pd.isna(value):
            return ""
        return f"{float(value):.2e}"
    except Exception:
        return ""


def format_q_values_for_display(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for col in out.columns:
        if "q value" in str(col).strip().casefold():
            out[col] = out[col].map(format_q_value_sci2)
    return out




def parse_single_gene_query(raw_query: str) -> tuple[str, bool]:
    """Return (single_gene, has_multiple_genes).

    The Morphological Profile tab accepts only one gene at a time.
    Treat comma-, semicolon-, or whitespace-separated entries as multiple genes.
    """
    txt = str(raw_query or "").strip().upper()
    if not txt:
        return "", False
                                                             
    if "," in txt or ";" in txt:
        parts = [g.strip() for g in txt.replace(";", ",").split(",") if g.strip()]
        return (parts[0] if parts else ""), len(parts) > 1
                                                                                   
    parts = txt.split()
    if len(parts) > 1:
        return parts[0], True
    return txt, False


def extract_significant_morphological_profile(
    sys_gene: str,
    q_threshold: float,
    z_values_df: pd.DataFrame,
    q_values_df: pd.DataFrame,
    parameter_desc_df: pd.DataFrame,
) -> pd.DataFrame:
    """Collect significant CalMorph parameters and their Z/q values for one gene."""
    if sys_gene not in q_values_df.index:
        return pd.DataFrame(columns=["Parameter", "Z value", "q value", "Description"])
    if sys_gene not in z_values_df.index:
        return pd.DataFrame(columns=["Parameter", "Z value", "q value", "Description"])

    q_row = pd.to_numeric(q_values_df.loc[sys_gene], errors="coerce")
    z_row = pd.to_numeric(z_values_df.loc[sys_gene], errors="coerce")

                                                                  
    common_params = [c for c in q_values_df.columns if c in z_values_df.columns]
    q_row = q_row.reindex(common_params)
    z_row = z_row.reindex(common_params)
    sig_params = q_row.index[(q_row <= float(q_threshold)) & q_row.notna()].tolist()

    if not sig_params:
        return pd.DataFrame(columns=["Parameter", "Z value", "q value", "Description"])

    desc_norm = normalize_parameter_description_table(parameter_desc_df)
    out = pd.DataFrame({
        "Parameter": sig_params,
        "Z value": [float(z_row[p]) if pd.notna(z_row[p]) else np.nan for p in sig_params],
        "q value": [float(q_row[p]) if pd.notna(q_row[p]) else np.nan for p in sig_params],
    })
                                                                                  
                                                       
    out["Parameter"] = out["Parameter"].map(display_calmorph_parameter_name)
    out["Description"] = out["Parameter"].map(desc_norm["Description"] if "Description" in desc_norm.columns else {})
    out["Description"] = out["Description"].fillna("")
    out = out.sort_values(["q value", "Parameter"], ascending=[True, True]).reset_index(drop=True)
    return out

def _canonical_calmorph_parameter_name(parameter: str) -> str:
    """Normalize CalMorph parameter names for user input matching.

    Users may type either the raw CalMorph name with dots or the app display
    version with dashes, for example DCV17.1_C or DCV17-1_C.
    """
    return str(parameter or "").strip().replace("-", ".").casefold()


def find_qvalue_parameter_column(parameter_raw: str, q_values_df: pd.DataFrame) -> str | None:
    """Return the qValues_Both.csv column matching a user-entered parameter name."""
    target = _canonical_calmorph_parameter_name(parameter_raw)
    if not target:
        return None

                                                                                           
    for col in q_values_df.columns:
        col_str = str(col).strip()
        if col_str == str(parameter_raw).strip() or display_calmorph_parameter_name(col_str) == str(parameter_raw).strip():
            return col

    for col in q_values_df.columns:
        if _canonical_calmorph_parameter_name(col) == target:
            return col
    return None


def find_mutants_for_calmorph_parameter(
    parameter_raw: str,
    q_threshold: float,
    q_values_df: pd.DataFrame,
    info_df: pd.DataFrame,
) -> tuple[pd.DataFrame, str | None]:
    """Return mutants whose q value passes the threshold for one CalMorph parameter."""
    param_col = find_qvalue_parameter_column(parameter_raw, q_values_df)
    out_cols = ["Systematic name", "Standard name", "SGD ID", "Name description", "q value"]
    if param_col is None:
        return pd.DataFrame(columns=out_cols), None

    q_series = pd.to_numeric(q_values_df[param_col], errors="coerce")
    sig_series = q_series[(q_series <= float(q_threshold)) & q_series.notna()].sort_values()

    rows = []
    for sys_gene, q_val in sig_series.items():
        systematic, standard, sgd_id = get_gene_display_names(str(sys_gene), info_df)
        rows.append({
            "Systematic name": systematic,
            "Standard name": standard,
            "SGD ID": sgd_id,
            "Name description": get_gene_name_description(str(sys_gene), info_df),
            "q value": float(q_val),
        })

    return pd.DataFrame(rows, columns=out_cols), str(param_col)


def parse_calmorph_parameter_queries(raw_parameters: str) -> list[str]:
    parts = [p.strip() for p in str(raw_parameters or "").split(",")]
    return [p for p in parts if p]


def find_mutants_for_calmorph_parameters(
    parameters_raw: list[str],
    q_threshold: float,
    q_values_df: pd.DataFrame,
    z_values_df: pd.DataFrame,
    info_df: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str], list[str]]:
    matched_parameters: list[str] = []
    missing_parameters: list[str] = []

    for parameter_raw in parameters_raw:
        matched_parameter = find_qvalue_parameter_column(parameter_raw, q_values_df)
        if matched_parameter is None:
            missing_parameters.append(parameter_raw)
        else:
            matched_parameters.append(str(matched_parameter))

    value_cols: list[str] = []
    for matched_parameter in matched_parameters:
        display_parameter = display_calmorph_parameter_name(matched_parameter)
        value_cols.append(f"Z value ({display_parameter})")
        value_cols.append(f"q value ({display_parameter})")

    out_cols = ["Systematic name", "Standard name", "SGD ID", "Name description"] + value_cols + ["Maximum q value"]

    if not matched_parameters:
        return pd.DataFrame(columns=out_cols), matched_parameters, missing_parameters

    significant_sets: list[set[str]] = []
    q_lookup: dict[str, pd.Series] = {}
    z_lookup: dict[str, pd.Series] = {}

    for matched_parameter in matched_parameters:
        q_series = pd.to_numeric(q_values_df[matched_parameter], errors="coerce")
        q_lookup[matched_parameter] = q_series

        if matched_parameter in z_values_df.columns:
            z_lookup[matched_parameter] = pd.to_numeric(z_values_df[matched_parameter], errors="coerce")
        else:
            z_lookup[matched_parameter] = pd.Series(np.nan, index=q_values_df.index)

        sig_genes = set(q_series[(q_series <= float(q_threshold)) & q_series.notna()].index.astype(str))
        significant_sets.append(sig_genes)

    shared_genes = set.intersection(*significant_sets) if significant_sets else set()
    rows = []
    for sys_gene in sorted(shared_genes):
        systematic, standard, sgd_id = get_gene_display_names(str(sys_gene), info_df)
        row = {
            "Systematic name": systematic,
            "Standard name": standard,
            "SGD ID": sgd_id,
            "Name description": get_gene_name_description(str(sys_gene), info_df),
        }
        q_values_for_gene = []
        for matched_parameter in matched_parameters:
            display_parameter = display_calmorph_parameter_name(matched_parameter)
            z_val = z_lookup[matched_parameter].loc[sys_gene] if sys_gene in z_lookup[matched_parameter].index else np.nan
            q_val = float(q_lookup[matched_parameter].loc[sys_gene])
            row[f"Z value ({display_parameter})"] = float(z_val) if pd.notna(z_val) else np.nan
            row[f"q value ({display_parameter})"] = q_val
            q_values_for_gene.append(q_val)
        row["Maximum q value"] = max(q_values_for_gene) if q_values_for_gene else np.nan
        rows.append(row)

    out = pd.DataFrame(rows, columns=out_cols)
    if not out.empty:
        out = out.sort_values(["Maximum q value", "Systematic name"], ascending=[True, True]).reset_index(drop=True)
    return out, matched_parameters, missing_parameters


def classify_calmorph_parameter(parameter: str) -> tuple[str | None, str | None]:
    """Map a CalMorph parameter to a cell component and cell-cycle stage."""
    p = str(parameter).strip()
    component = None
    if p.startswith("A"):
        component = "actin"
    elif p.startswith("C"):
        component = "cell_wall"
    elif p.startswith("D"):
        component = "nucleus"

    stage = None
    if p.endswith("_A1B"):
        stage = "S/G2"
    elif p.endswith("_A"):
        stage = "G1"
    elif p.endswith("_C"):
        stage = "M"
    return component, stage


def summarize_profile_colors(profile_df: pd.DataFrame) -> dict[str, dict[str, str]]:
    """Return component colors for G1, S/G2, and M according to Z-value direction.

    Default = gray when no significant parameter maps to that stage/component.
    Cell wall = green, actin = red, nucleus = blue.
    Positive mean Z = darker shade; negative mean Z = lighter shade.
    """
    stage_component_z: dict[str, dict[str, list[float]]] = {
        "G1": {"actin": [], "cell_wall": [], "nucleus": []},
        "S/G2": {"actin": [], "cell_wall": [], "nucleus": []},
        "M": {"actin": [], "cell_wall": [], "nucleus": []},
    }
    if profile_df is not None and len(profile_df) > 0:
        for _, row in profile_df.iterrows():
            comp, stage = classify_calmorph_parameter(row.get("Parameter", ""))
            if comp is None or stage is None:
                continue
            z = pd.to_numeric(row.get("Z value", np.nan), errors="coerce")
            if pd.notna(z):
                stage_component_z[stage][comp].append(float(z))

    palette = {
        "cell_wall": {"positive": "#1B5E20", "negative": "#C8E6C9"},         
        "actin": {"positive": "#B71C1C", "negative": "#FFCDD2"},           
        "nucleus": {"positive": "#0D47A1", "negative": "#BBDEFB"},          
    }

    colors = {stage: {comp: "#BDBDBD" for comp in comps} for stage, comps in stage_component_z.items()}
    for stage, comps in stage_component_z.items():
        for comp, zs in comps.items():
            if not zs:
                continue
            mean_z = float(np.nanmean(zs))
            if mean_z > 0:
                colors[stage][comp] = palette[comp]["positive"]
            elif mean_z < 0:
                colors[stage][comp] = palette[comp]["negative"]
    return colors

def profile_stage_component_parameters(profile_df: pd.DataFrame) -> dict[str, dict[str, list[dict[str, object]]]]:
    """Collect significant parameters by stage/component for SVG labels."""
    out: dict[str, dict[str, list[dict[str, object]]]] = {
        "G1": {"actin": [], "cell_wall": [], "nucleus": []},
        "S/G2": {"actin": [], "cell_wall": [], "nucleus": []},
        "M": {"actin": [], "cell_wall": [], "nucleus": []},
    }
    if profile_df is None or profile_df.empty:
        return out
    for _, row in profile_df.iterrows():
        comp, stage = classify_calmorph_parameter(row.get("Parameter", ""))
        if comp is None or stage is None:
            continue
        z = pd.to_numeric(row.get("Z value", np.nan), errors="coerce")
        q = pd.to_numeric(row.get("q value", np.nan), errors="coerce")
        if pd.isna(z):
            continue
        out[stage][comp].append({
            "parameter": str(row.get("Parameter", "")),
            "z": float(z),
            "q": float(q) if pd.notna(q) else np.nan,
        })
    for stage in out:
        for comp in out[stage]:
            out[stage][comp].sort(key=lambda d: (np.inf if pd.isna(d.get("q", np.nan)) else float(d.get("q", np.inf))))
    return out


SVG_TEMPLATES = {
    "G1": '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<!-- Created with Inkscape (http://www.inkscape.org/) -->\n\n<svg\n   width="38.770222mm"\n   height="30.69809mm"\n   viewBox="0 0 38.770221 30.698089"\n   version="1.1"\n   id="svg1"\n   xml:space="preserve"\n   inkscape:version="1.4.3 (0d15f75, 2025-12-25)"\n   sodipodi:docname="G1.svg"\n   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"\n   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"\n   xmlns="http://www.w3.org/2000/svg"\n   xmlns:svg="http://www.w3.org/2000/svg"><sodipodi:namedview\n     id="namedview1"\n     pagecolor="#ffffff"\n     bordercolor="#000000"\n     borderopacity="0.25"\n     inkscape:showpageshadow="2"\n     inkscape:pageopacity="0.0"\n     inkscape:pagecheckerboard="0"\n     inkscape:deskcolor="#d1d1d1"\n     inkscape:document-units="mm"\n     inkscape:zoom="3.8954759"\n     inkscape:cx="61.738285"\n     inkscape:cy="52.62515"\n     inkscape:window-width="1920"\n     inkscape:window-height="991"\n     inkscape:window-x="-9"\n     inkscape:window-y="-9"\n     inkscape:window-maximized="1"\n     inkscape:current-layer="layer1" /><defs\n     id="defs1"><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(29.944918,0,0,31.072665,124.11605,234.76371)"\n       spreadMethod="pad"\n       id="radialGradient3"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.202018,0,0,18.88752,168.72075,213.58708)"\n       spreadMethod="pad"\n       id="radialGradient6"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop4" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(14.748164,-6.5909696,2.7657445,8.1230594,145.03468,226.82087)"\n       spreadMethod="pad"\n       id="radialGradient9"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop7" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.4817462,-2.9884689,1.97163,-4.688695,495.58499,36.967648)"\n       spreadMethod="pad"\n       id="radialGradient12"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop10" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop11" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop12" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath12"><path\n         d="m 492.82501,41.714304 h 5.516 v -9.495701 h -5.516 z"\n         id="path12" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath14"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path14" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath57"><path\n         d="M -21.860801,-2.0250001e-4 H 20.816201 V -62.223204 h -42.677002 z"\n         id="path57" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath58"><path\n         d="m 620.04667,389.236 h 60.288 V 473.2 h -60.288 z"\n         id="path58" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath61"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path61" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath62"><path\n         d="M -42.421399,0 H 41.153604 V -18.934 H -42.421399 Z"\n         id="path62" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath63"><path\n         d="M 596.15867,447.332 H 708.592 v 26.24533 H 596.15867 Z"\n         id="path63" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath66"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path66" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath67"><path\n         d="M -17.009698,-5.0000001e-4 H 20.600303 V -48.434501 h -37.610001 z"\n         id="path67" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath68"><path\n         d="M 626.80133,408.33467 H 677.448 V 473.772 h -50.64667 z"\n         id="path68" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath71"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path71" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,136.87935,213.61111)"\n       spreadMethod="pad"\n       id="radialGradient75"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop74" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop75" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,169.12798,224.81072)"\n       spreadMethod="pad"\n       id="radialGradient77"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,143.58212,216.58798)"\n       spreadMethod="pad"\n       id="radialGradient79"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop78" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop79" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3591334,0,0,2.3591334,148.14712,241.9406)"\n       spreadMethod="pad"\n       id="radialGradient81"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop80" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop81" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,163.32091,225.16177)"\n       spreadMethod="pad"\n       id="radialGradient83"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop82" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop83" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,153.9732,234.14904)"\n       spreadMethod="pad"\n       id="radialGradient85"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop84" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop85" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,156.257,213.85752)"\n       spreadMethod="pad"\n       id="radialGradient87"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop86" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop87" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.088001,0,0,-17.088001,457.0029,611.94415)"\n       spreadMethod="pad"\n       id="radialGradient345"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop343" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop344" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop345" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(10.386942,0,0,-10.386942,482.45648,623.58997)"\n       spreadMethod="pad"\n       id="radialGradient348"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop346" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop347" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop348" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.416008,3.6246166,1.578266,-4.4671698,468.94009,616.31226)"\n       spreadMethod="pad"\n       id="radialGradient351"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop349" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop350" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop351" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.888271,-2.191292,1.500142,-3.4379811,474.04248,618.56189)"\n       spreadMethod="pad"\n       id="radialGradient354"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop352" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop353" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop354" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath354"><path\n         d="m 471.94301,622.04202 h 4.197 v -6.962 h -4.197 z"\n         id="path354" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath356"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path356" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.129417,0,0,-17.681181,455.73108,563.27264)"\n       spreadMethod="pad"\n       id="radialGradient360"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop358" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop359" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop360" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.3039217,0,0,-5.4747686,453.74182,559.15039)"\n       spreadMethod="pad"\n       id="radialGradient363"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop361" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop362" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop363" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(9.7889442,0,0,-9.7889442,474.04147,535.25586)"\n       spreadMethod="pad"\n       id="radialGradient365"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop364" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop365" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(13.859965,2.4953344,2.3838298,-14.508271,529.60748,590.08429)"\n       spreadMethod="pad"\n       id="radialGradient398"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop396" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop397" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop398" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(16.999365,0,0,-17.794517,498.96738,577.06146)"\n       spreadMethod="pad"\n       id="radialGradient401"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop399" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop400" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop401" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.2615695,0,0,-5.5076814,497.80457,574.14221)"\n       spreadMethod="pad"\n       id="radialGradient404"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop402" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop403" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop404" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(3.8265915,0.6943604,0.6633328,-4.0055819,530.68677,590.08441)"\n       spreadMethod="pad"\n       id="radialGradient407"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop405" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop406" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop407" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.9543619,-2.4196851,1.5915818,-3.7604923,516.02515,583.75952)"\n       spreadMethod="pad"\n       id="radialGradient410"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop408" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop409" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop410" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath410"><path\n         d="m 513.83501,587.57602 h 4.379 v -7.636 h -4.379 z"\n         id="path410" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath412"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path412" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,469.0929,479.88107)"\n       spreadMethod="pad"\n       id="radialGradient419"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop417" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop418" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop419" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,467.53519,476.02194)"\n       spreadMethod="pad"\n       id="radialGradient422"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop420" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop421" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop422" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,507.55173,439.42313)"\n       spreadMethod="pad"\n       id="radialGradient427"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop425" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop426" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop427" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,505.99402,435.564)"\n       spreadMethod="pad"\n       id="radialGradient430"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop428" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop429" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop430" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(19.683309,-3.5213403,3.3854092,20.473634,151.17283,330.66865)"\n       spreadMethod="pad"\n       id="radialGradient436"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop434" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop435" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop436" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.141743,0,0,25.111086,107.65912,349.04613)"\n       spreadMethod="pad"\n       id="radialGradient439"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop437" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop438" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop439" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4722473,0,0,7.7722738,114.72071,492.57303)"\n       spreadMethod="pad"\n       id="radialGradient442"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop440" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop441" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop442" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.4343553,-0.97986038,0.94203584,5.6525553,161.41859,470.0759)"\n       spreadMethod="pad"\n       id="radialGradient445"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop443" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop444" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop445" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0816238,-2.5609398,1.6952205,-3.98002,520.40234,194.37172)"\n       spreadMethod="pad"\n       id="radialGradient448"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop446" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop447" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop448" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath448"><path\n         d="m 518.07001,198.41101 h 4.663 v -8.082 h -4.663 z"\n         id="path448" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath450"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path450" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.094729,0,0,24.870851,124.129,154.78121)"\n       spreadMethod="pad"\n       id="radialGradient456"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop454" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop455" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop456" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4606481,0,0,7.7009657,113.04409,52.125052)"\n       spreadMethod="pad"\n       id="radialGradient459"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop457" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop458" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop459" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.150868,0,0,-18.834442,451.0462,147.19133)"\n       spreadMethod="pad"\n       id="radialGradient462"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop460" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop461" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop462" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(11.033005,0,0,-11.448516,478.08298,160.02737)"\n       spreadMethod="pad"\n       id="radialGradient465"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop463" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop464" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop465" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.9394798,3.9950626,1.6764334,-4.9237266,463.72586,152.00581)"\n       spreadMethod="pad"\n       id="radialGradient468"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop466" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop467" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop468" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0057206,-2.4152482,1.5934501,-3.7893524,469.14563,154.48541)"\n       spreadMethod="pad"\n       id="radialGradient471"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop469" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop470" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop471" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath471"><path\n         d="m 466.91501,158.32201 h 4.458 v -7.675 h -4.458 z"\n         id="path471" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath473"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path473" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,160.06101,208.15524)"\n       spreadMethod="pad"\n       id="radialGradient77-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,142.82225,469.66275)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,144.99955,477.22245)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,152.21191,481.3324)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,131.30639,474.06771)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.75035,480.51167)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-38"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-7" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,138.83947,489.85995)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-19"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-86" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.16259,44.279802)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,126.77659,60.979642)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,121.51251,31.391882)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-0" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,92.47587,51.717192)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-7" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,108.08659,64.332832)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-84" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,128.59835,171.61115)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-5" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,122.17299,137.41595)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,102.38731,35.583232)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-4"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(29.944918,0,0,31.072665,124.11605,234.76371)"\n       spreadMethod="pad"\n       id="radialGradient3-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop1-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop2-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop3-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.202018,0,0,18.88752,168.72075,213.58708)"\n       spreadMethod="pad"\n       id="radialGradient6-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop4-7" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop5-3" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop6-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(14.748164,-6.5909696,2.7657445,8.1230594,145.03468,226.82087)"\n       spreadMethod="pad"\n       id="radialGradient9-5"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop7-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop8-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop9-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.4817462,-2.9884689,1.97163,-4.688695,495.58499,36.967648)"\n       spreadMethod="pad"\n       id="radialGradient12-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop10-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop11-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop12-6" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath12-6"><path\n         d="m 492.82501,41.714304 h 5.516 v -9.495701 h -5.516 z"\n         id="path12-0" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath14-3"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path14-8" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath57-0"><path\n         d="M -21.860801,-2.0250001e-4 H 20.816201 V -62.223204 h -42.677002 z"\n         id="path57-1" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath58-2"><path\n         d="m 620.04667,389.236 h 60.288 V 473.2 h -60.288 z"\n         id="path58-5" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath61-0"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path61-9" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath62-4"><path\n         d="M -42.421399,0 H 41.153604 V -18.934 H -42.421399 Z"\n         id="path62-7" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath63-8"><path\n         d="M 596.15867,447.332 H 708.592 v 26.24533 H 596.15867 Z"\n         id="path63-3" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath66-5"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path66-1" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath67-2"><path\n         d="M -17.009698,-5.0000001e-4 H 20.600303 V -48.434501 h -37.610001 z"\n         id="path67-0" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath68-1"><path\n         d="M 626.80133,408.33467 H 677.448 V 473.772 h -50.64667 z"\n         id="path68-6" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath71-4"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path71-0" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,136.87935,213.61111)"\n       spreadMethod="pad"\n       id="radialGradient75-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop74-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop75-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,169.12798,224.81072)"\n       spreadMethod="pad"\n       id="radialGradient77-9"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,143.58212,216.58798)"\n       spreadMethod="pad"\n       id="radialGradient79-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop78-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop79-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3591334,0,0,2.3591334,148.14712,241.9406)"\n       spreadMethod="pad"\n       id="radialGradient81-9"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop80-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop81-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,163.32091,225.16177)"\n       spreadMethod="pad"\n       id="radialGradient83-0"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop82-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop83-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,153.9732,234.14904)"\n       spreadMethod="pad"\n       id="radialGradient85-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop84-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop85-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,156.257,213.85752)"\n       spreadMethod="pad"\n       id="radialGradient87-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop86-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop87-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.088001,0,0,-17.088001,457.0029,611.94415)"\n       spreadMethod="pad"\n       id="radialGradient345-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop343-0" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop344-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop345-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(10.386942,0,0,-10.386942,482.45648,623.58997)"\n       spreadMethod="pad"\n       id="radialGradient348-4"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop346-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop347-6" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop348-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.416008,3.6246166,1.578266,-4.4671698,468.94009,616.31226)"\n       spreadMethod="pad"\n       id="radialGradient351-1"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop349-0" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop350-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop351-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.888271,-2.191292,1.500142,-3.4379811,474.04248,618.56189)"\n       spreadMethod="pad"\n       id="radialGradient354-6"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop352-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop353-3" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop354-4" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath354-8"><path\n         d="m 471.94301,622.04202 h 4.197 v -6.962 h -4.197 z"\n         id="path354-4" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath356-9"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path356-9" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.129417,0,0,-17.681181,455.73108,563.27264)"\n       spreadMethod="pad"\n       id="radialGradient360-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop358-5" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop359-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop360-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.3039217,0,0,-5.4747686,453.74182,559.15039)"\n       spreadMethod="pad"\n       id="radialGradient363-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop361-3" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop362-7" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop363-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(9.7889442,0,0,-9.7889442,474.04147,535.25586)"\n       spreadMethod="pad"\n       id="radialGradient365-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop364-8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop365-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(13.859965,2.4953344,2.3838298,-14.508271,529.60748,590.08429)"\n       spreadMethod="pad"\n       id="radialGradient398-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop396-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop397-0" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop398-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(16.999365,0,0,-17.794517,498.96738,577.06146)"\n       spreadMethod="pad"\n       id="radialGradient401-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop399-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop400-9" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop401-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.2615695,0,0,-5.5076814,497.80457,574.14221)"\n       spreadMethod="pad"\n       id="radialGradient404-9"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop402-7" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop403-2" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop404-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(3.8265915,0.6943604,0.6633328,-4.0055819,530.68677,590.08441)"\n       spreadMethod="pad"\n       id="radialGradient407-8"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop405-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop406-8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop407-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.9543619,-2.4196851,1.5915818,-3.7604923,516.02515,583.75952)"\n       spreadMethod="pad"\n       id="radialGradient410-0"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop408-7" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop409-8" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop410-1" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath410-5"><path\n         d="m 513.83501,587.57602 h 4.379 v -7.636 h -4.379 z"\n         id="path410-8" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath412-6"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path412-1" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,469.0929,479.88107)"\n       spreadMethod="pad"\n       id="radialGradient419-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop417-4" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop418-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop419-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,467.53519,476.02194)"\n       spreadMethod="pad"\n       id="radialGradient422-8"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop420-6" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop421-2" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop422-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,507.55173,439.42313)"\n       spreadMethod="pad"\n       id="radialGradient427-5"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop425-3" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop426-9" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop427-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,505.99402,435.564)"\n       spreadMethod="pad"\n       id="radialGradient430-4"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop428-6" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop429-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop430-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(19.683309,-3.5213403,3.3854092,20.473634,151.17283,330.66865)"\n       spreadMethod="pad"\n       id="radialGradient436-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop434-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop435-1" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop436-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.141743,0,0,25.111086,107.65912,349.04613)"\n       spreadMethod="pad"\n       id="radialGradient439-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop437-6" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop438-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop439-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4722473,0,0,7.7722738,114.72071,492.57303)"\n       spreadMethod="pad"\n       id="radialGradient442-5"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop440-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop441-0" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop442-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.4343553,-0.97986038,0.94203584,5.6525553,161.41859,470.0759)"\n       spreadMethod="pad"\n       id="radialGradient445-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop443-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop444-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop445-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0816238,-2.5609398,1.6952205,-3.98002,520.40234,194.37172)"\n       spreadMethod="pad"\n       id="radialGradient448-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop446-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop447-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop448-3" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath448-2"><path\n         d="m 518.07001,198.41101 h 4.663 v -8.082 h -4.663 z"\n         id="path448-5" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath450-2"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path450-5" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.094729,0,0,24.870851,124.129,154.78121)"\n       spreadMethod="pad"\n       id="radialGradient456-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop454-6" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop455-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop456-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4606481,0,0,7.7009657,113.04409,52.125052)"\n       spreadMethod="pad"\n       id="radialGradient459-2"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop457-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop458-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop459-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.150868,0,0,-18.834442,451.0462,147.19133)"\n       spreadMethod="pad"\n       id="radialGradient462-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop460-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop461-6" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop462-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(11.033005,0,0,-11.448516,478.08298,160.02737)"\n       spreadMethod="pad"\n       id="radialGradient465-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop463-2" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop464-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop465-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.9394798,3.9950626,1.6764334,-4.9237266,463.72586,152.00581)"\n       spreadMethod="pad"\n       id="radialGradient468-4"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop466-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop467-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop468-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0057206,-2.4152482,1.5934501,-3.7893524,469.14563,154.48541)"\n       spreadMethod="pad"\n       id="radialGradient471-5"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop469-0" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop470-8" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop471-3" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath471-9"><path\n         d="m 466.91501,158.32201 h 4.458 v -7.675 h -4.458 z"\n         id="path471-3" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath473-9"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path473-6" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,160.06101,208.15524)"\n       spreadMethod="pad"\n       id="radialGradient77-2-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,142.82225,469.66275)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,144.99955,477.22245)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-35"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-76" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-65" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,152.21191,481.3324)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-1-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-8-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-9-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,131.30639,474.06771)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-3-4"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-1-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-2-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.75035,480.51167)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-38-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-7-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-4-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,138.83947,489.85995)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-19-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-86-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-5-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.16259,44.279802)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-4-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-8-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,126.77659,60.979642)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,121.51251,31.391882)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-6-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-0-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-6-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,92.47587,51.717192)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-3-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-7-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-8-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,108.08659,64.332832)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-5-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-9-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-84-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,128.59835,171.61115)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-1-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-5-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-4-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,122.17299,137.41595)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-6-0" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-0-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,102.38731,35.583232)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-4-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-1-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-2-9" /></radialGradient></defs><g\n     inkscape:label="Layer 1"\n     inkscape:groupmode="layer"\n     id="layer1"\n     transform="translate(-65.082567,-72.465926)"><g\n       inkscape:label="Layer 1"\n       id="layer1-2"\n       transform="translate(-163.1037,62.614794)" /><g\n       id="g4"\n       transform="matrix(0.63213006,0,0,0.62189589,12.43297,57.672216)"><path\n         id="path456"\n         d="m 111.46638,27.139752 c -14.84134,2.44 -25.424003,13.77867 -23.64267,25.32134 v 0 c 1.786667,11.55066 15.26,18.932 30.10133,16.496 v 0 c 14.83467,-2.43467 25.412,-13.77467 23.632,-25.32934 v 0 c -1.55333,-10.05866 -11.97466,-16.95466 -24.44,-16.95333 v 0 c -1.84666,0 -3.73733,0.15067 -5.65066,0.46533"\n         style="fill:none;stroke:#006200;stroke-width:1.33333;stroke-opacity:1" /><path\n         id="path77-7-4-5-0"\n         d="m 137.05176,41.900952 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-2-5);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-93"\n         d="m 126.66576,58.600792 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-0-8);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-93-9"\n         d="m 92.36504,49.338342 c -1.28801,0.0587 -2.286671,1.17334 -2.225341,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376001,0.884 1.254671,1.48133 2.257341,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-0-3-8);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-93-1"\n         d="m 102.27648,33.204382 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-0-4-6);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-93-6"\n         d="m 107.97576,61.953982 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-0-5-1);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-84"\n         d="m 121.40168,29.013032 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-6-3);stroke:none;stroke-width:1.33333" /><path\n         id="path459"\n         d="m 112.77038,46.078422 c -4.84,0.23733 -8.64134,3.132 -8.488,6.47067 v 0 c 0.032,0.75333 0.26666,1.46666 0.664,2.116 v 0 c 1.35333,2.224 4.628,3.68933 8.37466,3.50533 v 0 c 3.74934,-0.17733 6.87467,-1.95733 8.016,-4.30533 v 0 c 0.33734,-0.67867 0.50267,-1.412 0.468,-2.164 v 0 c -0.14666,-3.20134 -3.856,-5.63867 -8.412,-5.63867 v 0 c -0.20533,0 -0.41333,0.005 -0.62266,0.016"\n         style="fill:url(#radialGradient459-2);stroke:none;stroke-width:1.33333" /></g></g></svg>\n',
    "S/G2": '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<!-- Created with Inkscape (http://www.inkscape.org/) -->\n\n<svg\n   width="54.103817mm"\n   height="42.456009mm"\n   viewBox="0 0 54.103818 42.45601"\n   version="1.1"\n   id="svg1"\n   xml:space="preserve"\n   inkscape:version="1.4.3 (0d15f75, 2025-12-25)"\n   sodipodi:docname="SG2.svg"\n   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"\n   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"\n   xmlns="http://www.w3.org/2000/svg"\n   xmlns:svg="http://www.w3.org/2000/svg"><sodipodi:namedview\n     id="namedview1"\n     pagecolor="#ffffff"\n     bordercolor="#000000"\n     borderopacity="0.25"\n     inkscape:showpageshadow="2"\n     inkscape:pageopacity="0.0"\n     inkscape:pagecheckerboard="0"\n     inkscape:deskcolor="#d1d1d1"\n     inkscape:document-units="mm"\n     inkscape:zoom="2.7545174"\n     inkscape:cx="5.26408"\n     inkscape:cy="98.383841"\n     inkscape:window-width="1920"\n     inkscape:window-height="991"\n     inkscape:window-x="-9"\n     inkscape:window-y="-9"\n     inkscape:window-maximized="1"\n     inkscape:current-layer="layer1" /><defs\n     id="defs1"><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(29.944918,0,0,31.072665,124.11605,234.76371)"\n       spreadMethod="pad"\n       id="radialGradient3"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.202018,0,0,18.88752,168.72075,213.58708)"\n       spreadMethod="pad"\n       id="radialGradient6"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop4" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(14.748164,-6.5909696,2.7657445,8.1230594,145.03468,226.82087)"\n       spreadMethod="pad"\n       id="radialGradient9"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop7" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.4817462,-2.9884689,1.97163,-4.688695,495.58499,36.967648)"\n       spreadMethod="pad"\n       id="radialGradient12"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop10" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop11" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop12" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath12"><path\n         d="m 492.82501,41.714304 h 5.516 v -9.495701 h -5.516 z"\n         id="path12" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath14"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path14" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath57"><path\n         d="M -21.860801,-2.0250001e-4 H 20.816201 V -62.223204 h -42.677002 z"\n         id="path57" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath58"><path\n         d="m 620.04667,389.236 h 60.288 V 473.2 h -60.288 z"\n         id="path58" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath61"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path61" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath62"><path\n         d="M -42.421399,0 H 41.153604 V -18.934 H -42.421399 Z"\n         id="path62" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath63"><path\n         d="M 596.15867,447.332 H 708.592 v 26.24533 H 596.15867 Z"\n         id="path63" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath66"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path66" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath67"><path\n         d="M -17.009698,-5.0000001e-4 H 20.600303 V -48.434501 h -37.610001 z"\n         id="path67" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath68"><path\n         d="M 626.80133,408.33467 H 677.448 V 473.772 h -50.64667 z"\n         id="path68" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath71"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path71" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,136.87935,213.61111)"\n       spreadMethod="pad"\n       id="radialGradient75"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop74" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop75" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,169.12798,224.81072)"\n       spreadMethod="pad"\n       id="radialGradient77"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,143.58212,216.58798)"\n       spreadMethod="pad"\n       id="radialGradient79"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop78" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop79" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3591334,0,0,2.3591334,148.14712,241.9406)"\n       spreadMethod="pad"\n       id="radialGradient81"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop80" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop81" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,163.32091,225.16177)"\n       spreadMethod="pad"\n       id="radialGradient83"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop82" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop83" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,153.9732,234.14904)"\n       spreadMethod="pad"\n       id="radialGradient85"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop84" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop85" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,156.257,213.85752)"\n       spreadMethod="pad"\n       id="radialGradient87"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop86" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop87" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.088001,0,0,-17.088001,457.0029,611.94415)"\n       spreadMethod="pad"\n       id="radialGradient345"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop343" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop344" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop345" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(10.386942,0,0,-10.386942,482.45648,623.58997)"\n       spreadMethod="pad"\n       id="radialGradient348"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop346" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop347" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop348" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.416008,3.6246166,1.578266,-4.4671698,468.94009,616.31226)"\n       spreadMethod="pad"\n       id="radialGradient351"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop349" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop350" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop351" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.888271,-2.191292,1.500142,-3.4379811,474.04248,618.56189)"\n       spreadMethod="pad"\n       id="radialGradient354"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop352" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop353" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop354" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath354"><path\n         d="m 471.94301,622.04202 h 4.197 v -6.962 h -4.197 z"\n         id="path354" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath356"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path356" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.129417,0,0,-17.681181,455.73108,563.27264)"\n       spreadMethod="pad"\n       id="radialGradient360"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop358" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop359" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop360" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.3039217,0,0,-5.4747686,453.74182,559.15039)"\n       spreadMethod="pad"\n       id="radialGradient363"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop361" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop362" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop363" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(9.7889442,0,0,-9.7889442,474.04147,535.25586)"\n       spreadMethod="pad"\n       id="radialGradient365"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop364" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop365" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(13.859965,2.4953344,2.3838298,-14.508271,529.60748,590.08429)"\n       spreadMethod="pad"\n       id="radialGradient398"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop396" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop397" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop398" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(16.999365,0,0,-17.794517,498.96738,577.06146)"\n       spreadMethod="pad"\n       id="radialGradient401"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop399" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop400" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop401" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.2615695,0,0,-5.5076814,497.80457,574.14221)"\n       spreadMethod="pad"\n       id="radialGradient404"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop402" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop403" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop404" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(3.8265915,0.6943604,0.6633328,-4.0055819,530.68677,590.08441)"\n       spreadMethod="pad"\n       id="radialGradient407"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop405" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop406" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop407" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.9543619,-2.4196851,1.5915818,-3.7604923,516.02515,583.75952)"\n       spreadMethod="pad"\n       id="radialGradient410"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop408" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop409" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop410" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath410"><path\n         d="m 513.83501,587.57602 h 4.379 v -7.636 h -4.379 z"\n         id="path410" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath412"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path412" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,469.0929,479.88107)"\n       spreadMethod="pad"\n       id="radialGradient419"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop417" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop418" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop419" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,467.53519,476.02194)"\n       spreadMethod="pad"\n       id="radialGradient422"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop420" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop421" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop422" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,507.55173,439.42313)"\n       spreadMethod="pad"\n       id="radialGradient427"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop425" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop426" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop427" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,505.99402,435.564)"\n       spreadMethod="pad"\n       id="radialGradient430"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop428" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop429" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop430" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(19.683309,-3.5213403,3.3854092,20.473634,151.17283,330.66865)"\n       spreadMethod="pad"\n       id="radialGradient436"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop434" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop435" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop436" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.141743,0,0,25.111086,107.65912,349.04613)"\n       spreadMethod="pad"\n       id="radialGradient439"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop437" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop438" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop439" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4722473,0,0,7.7722738,114.72071,492.57303)"\n       spreadMethod="pad"\n       id="radialGradient442"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop440" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop441" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop442" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.4343553,-0.97986038,0.94203584,5.6525553,161.41859,470.0759)"\n       spreadMethod="pad"\n       id="radialGradient445"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop443" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop444" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop445" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0816238,-2.5609398,1.6952205,-3.98002,520.40234,194.37172)"\n       spreadMethod="pad"\n       id="radialGradient448"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop446" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop447" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop448" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath448"><path\n         d="m 518.07001,198.41101 h 4.663 v -8.082 h -4.663 z"\n         id="path448" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath450"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path450" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.094729,0,0,24.870851,124.129,154.78121)"\n       spreadMethod="pad"\n       id="radialGradient456"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop454" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop455" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop456" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4606481,0,0,7.7009657,113.04409,52.125052)"\n       spreadMethod="pad"\n       id="radialGradient459"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop457" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop458" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop459" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.150868,0,0,-18.834442,451.0462,147.19133)"\n       spreadMethod="pad"\n       id="radialGradient462"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop460" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop461" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop462" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(11.033005,0,0,-11.448516,478.08298,160.02737)"\n       spreadMethod="pad"\n       id="radialGradient465"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop463" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop464" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop465" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.9394798,3.9950626,1.6764334,-4.9237266,463.72586,152.00581)"\n       spreadMethod="pad"\n       id="radialGradient468"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop466" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop467" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop468" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0057206,-2.4152482,1.5934501,-3.7893524,469.14563,154.48541)"\n       spreadMethod="pad"\n       id="radialGradient471"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop469" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop470" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop471" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath471"><path\n         d="m 466.91501,158.32201 h 4.458 v -7.675 h -4.458 z"\n         id="path471" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath473"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path473" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,160.06101,208.15524)"\n       spreadMethod="pad"\n       id="radialGradient77-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,142.82225,469.66275)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,144.99955,477.22245)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,152.21191,481.3324)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,131.30639,474.06771)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.75035,480.51167)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-38"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-7" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,138.83947,489.85995)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-19"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-86" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.16259,44.279802)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,126.77659,60.979642)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,121.51251,31.391882)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-0" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,92.47587,51.717192)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-7" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,108.08659,64.332832)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-84" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,128.59835,171.61115)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-5" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,122.17299,137.41595)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,102.38731,35.583232)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-4"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(29.944918,0,0,31.072665,124.11605,234.76371)"\n       spreadMethod="pad"\n       id="radialGradient3-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop1-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop2-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop3-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.202018,0,0,18.88752,168.72075,213.58708)"\n       spreadMethod="pad"\n       id="radialGradient6-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop4-7" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop5-3" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop6-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(14.748164,-6.5909696,2.7657445,8.1230594,145.03468,226.82087)"\n       spreadMethod="pad"\n       id="radialGradient9-5"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop7-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop8-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop9-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.4817462,-2.9884689,1.97163,-4.688695,495.58499,36.967648)"\n       spreadMethod="pad"\n       id="radialGradient12-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop10-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop11-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop12-6" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath12-6"><path\n         d="m 492.82501,41.714304 h 5.516 v -9.495701 h -5.516 z"\n         id="path12-0" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath14-3"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path14-8" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath57-0"><path\n         d="M -21.860801,-2.0250001e-4 H 20.816201 V -62.223204 h -42.677002 z"\n         id="path57-1" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath58-2"><path\n         d="m 620.04667,389.236 h 60.288 V 473.2 h -60.288 z"\n         id="path58-5" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath61-0"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path61-9" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath62-4"><path\n         d="M -42.421399,0 H 41.153604 V -18.934 H -42.421399 Z"\n         id="path62-7" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath63-8"><path\n         d="M 596.15867,447.332 H 708.592 v 26.24533 H 596.15867 Z"\n         id="path63-3" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath66-5"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path66-1" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath67-2"><path\n         d="M -17.009698,-5.0000001e-4 H 20.600303 V -48.434501 h -37.610001 z"\n         id="path67-0" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath68-1"><path\n         d="M 626.80133,408.33467 H 677.448 V 473.772 h -50.64667 z"\n         id="path68-6" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath71-4"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path71-0" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,136.87935,213.61111)"\n       spreadMethod="pad"\n       id="radialGradient75-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop74-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop75-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,169.12798,224.81072)"\n       spreadMethod="pad"\n       id="radialGradient77-9"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,143.58212,216.58798)"\n       spreadMethod="pad"\n       id="radialGradient79-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop78-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop79-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3591334,0,0,2.3591334,148.14712,241.9406)"\n       spreadMethod="pad"\n       id="radialGradient81-9"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop80-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop81-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,163.32091,225.16177)"\n       spreadMethod="pad"\n       id="radialGradient83-0"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop82-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop83-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,153.9732,234.14904)"\n       spreadMethod="pad"\n       id="radialGradient85-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop84-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop85-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,156.257,213.85752)"\n       spreadMethod="pad"\n       id="radialGradient87-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop86-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop87-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.088001,0,0,-17.088001,457.0029,611.94415)"\n       spreadMethod="pad"\n       id="radialGradient345-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop343-0" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop344-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop345-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(10.386942,0,0,-10.386942,482.45648,623.58997)"\n       spreadMethod="pad"\n       id="radialGradient348-4"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop346-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop347-6" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop348-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.416008,3.6246166,1.578266,-4.4671698,468.94009,616.31226)"\n       spreadMethod="pad"\n       id="radialGradient351-1"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop349-0" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop350-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop351-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.888271,-2.191292,1.500142,-3.4379811,474.04248,618.56189)"\n       spreadMethod="pad"\n       id="radialGradient354-6"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop352-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop353-3" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop354-4" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath354-8"><path\n         d="m 471.94301,622.04202 h 4.197 v -6.962 h -4.197 z"\n         id="path354-4" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath356-9"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path356-9" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.129417,0,0,-17.681181,455.73108,563.27264)"\n       spreadMethod="pad"\n       id="radialGradient360-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop358-5" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop359-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop360-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.3039217,0,0,-5.4747686,453.74182,559.15039)"\n       spreadMethod="pad"\n       id="radialGradient363-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop361-3" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop362-7" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop363-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(9.7889442,0,0,-9.7889442,474.04147,535.25586)"\n       spreadMethod="pad"\n       id="radialGradient365-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop364-8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop365-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(13.859965,2.4953344,2.3838298,-14.508271,529.60748,590.08429)"\n       spreadMethod="pad"\n       id="radialGradient398-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop396-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop397-0" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop398-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(16.999365,0,0,-17.794517,498.96738,577.06146)"\n       spreadMethod="pad"\n       id="radialGradient401-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop399-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop400-9" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop401-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.2615695,0,0,-5.5076814,497.80457,574.14221)"\n       spreadMethod="pad"\n       id="radialGradient404-9"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop402-7" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop403-2" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop404-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(3.8265915,0.6943604,0.6633328,-4.0055819,530.68677,590.08441)"\n       spreadMethod="pad"\n       id="radialGradient407-8"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop405-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop406-8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop407-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.9543619,-2.4196851,1.5915818,-3.7604923,516.02515,583.75952)"\n       spreadMethod="pad"\n       id="radialGradient410-0"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop408-7" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop409-8" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop410-1" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath410-5"><path\n         d="m 513.83501,587.57602 h 4.379 v -7.636 h -4.379 z"\n         id="path410-8" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath412-6"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path412-1" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,469.0929,479.88107)"\n       spreadMethod="pad"\n       id="radialGradient419-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop417-4" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop418-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop419-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,467.53519,476.02194)"\n       spreadMethod="pad"\n       id="radialGradient422-8"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop420-6" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop421-2" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop422-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,507.55173,439.42313)"\n       spreadMethod="pad"\n       id="radialGradient427-5"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop425-3" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop426-9" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop427-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,505.99402,435.564)"\n       spreadMethod="pad"\n       id="radialGradient430-4"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop428-6" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop429-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop430-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(19.683309,-3.5213403,3.3854092,20.473634,151.17283,330.66865)"\n       spreadMethod="pad"\n       id="radialGradient436-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop434-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop435-1" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop436-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.141743,0,0,25.111086,107.65912,349.04613)"\n       spreadMethod="pad"\n       id="radialGradient439-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop437-6" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop438-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop439-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4722473,0,0,7.7722738,114.72071,492.57303)"\n       spreadMethod="pad"\n       id="radialGradient442-5"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop440-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop441-0" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop442-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.4343553,-0.97986038,0.94203584,5.6525553,161.41859,470.0759)"\n       spreadMethod="pad"\n       id="radialGradient445-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop443-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop444-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop445-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0816238,-2.5609398,1.6952205,-3.98002,520.40234,194.37172)"\n       spreadMethod="pad"\n       id="radialGradient448-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop446-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop447-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop448-3" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath448-2"><path\n         d="m 518.07001,198.41101 h 4.663 v -8.082 h -4.663 z"\n         id="path448-5" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath450-2"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path450-5" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.094729,0,0,24.870851,124.129,154.78121)"\n       spreadMethod="pad"\n       id="radialGradient456-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop454-6" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop455-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop456-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4606481,0,0,7.7009657,113.04409,52.125052)"\n       spreadMethod="pad"\n       id="radialGradient459-2"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop457-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop458-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop459-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.150868,0,0,-18.834442,451.0462,147.19133)"\n       spreadMethod="pad"\n       id="radialGradient462-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop460-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop461-6" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop462-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(11.033005,0,0,-11.448516,478.08298,160.02737)"\n       spreadMethod="pad"\n       id="radialGradient465-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop463-2" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop464-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop465-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.9394798,3.9950626,1.6764334,-4.9237266,463.72586,152.00581)"\n       spreadMethod="pad"\n       id="radialGradient468-4"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop466-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop467-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop468-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0057206,-2.4152482,1.5934501,-3.7893524,469.14563,154.48541)"\n       spreadMethod="pad"\n       id="radialGradient471-5"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop469-0" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop470-8" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop471-3" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath471-9"><path\n         d="m 466.91501,158.32201 h 4.458 v -7.675 h -4.458 z"\n         id="path471-3" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath473-9"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path473-6" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,160.06101,208.15524)"\n       spreadMethod="pad"\n       id="radialGradient77-2-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,142.82225,469.66275)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,144.99955,477.22245)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-35"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-76" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-65" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,152.21191,481.3324)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-1-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-8-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-9-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,131.30639,474.06771)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-3-4"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-1-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-2-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.75035,480.51167)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-38-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-7-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-4-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,138.83947,489.85995)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-19-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-86-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-5-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.16259,44.279802)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-4-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-8-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,126.77659,60.979642)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,121.51251,31.391882)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-6-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-0-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-6-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,92.47587,51.717192)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-3-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-7-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-8-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,108.08659,64.332832)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-5-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-9-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-84-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,128.59835,171.61115)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-1-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-5-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-4-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,122.17299,137.41595)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-6-0" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-0-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,102.38731,35.583232)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-4-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-1-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-2-9" /></radialGradient></defs><g\n     inkscape:label="Layer 1"\n     inkscape:groupmode="layer"\n     id="layer1"\n     transform="translate(-91.443726,-115.64942)"><g\n       id="g2"\n       transform="matrix(0.50412692,0,0,0.60813557,48.358713,-1.9542538)"><path\n         id="path3"\n         d="m 120.10414,208.63228 c -18.43467,3.052 -31.585337,17.22267 -29.37067,31.65334 v 0 c 2.226666,14.43333 18.964,23.656 37.39467,20.60933 v 0 c 18.43466,-3.04667 31.58533,-17.21733 29.37066,-31.652 v 0 c -1.928,-12.576 -14.88133,-21.192 -30.37333,-21.192 v 0 c -2.29467,0 -4.644,0.18933 -7.02133,0.58133"\n         style="fill:none;stroke:#006200;stroke-width:1.33333;stroke-opacity:1" /><path\n         id="path6"\n         d="m 161.04284,200.34731 c -8.49626,6.60073 -11.435,17.24722 -6.54755,23.78136 v 0 c 4.88619,6.53415 15.72014,6.47866 24.22894,-0.12454 v 0 c 8.48371,-6.59949 11.41618,-17.24969 6.53625,-23.77644 v 0 c -2.42678,-3.24795 -6.32544,-4.86946 -10.67333,-4.86946 v 0 c -4.4056,0 -9.27297,1.66344 -13.54431,4.98908"\n         style="fill:none;stroke:#006200;stroke-width:1.24389;stroke-opacity:1" /><path\n         id="path9"\n         d="m 142.31747,220.68428 c -9.49067,4.52534 -15.964,10.93734 -14.44667,14.32267 v 0 c 0.332,0.768 1.04667,1.31733 2.06134,1.656 v 0 c 3.48266,1.14533 10.46933,-0.19867 17.81199,-3.70267 v 0 c 7.35334,-3.508 12.888,-8.14266 14.30134,-11.62933 v 0 c 0.432,-1.01333 0.484,-1.93333 0.15466,-2.69867 v 0 c -0.61866,-1.37866 -2.45733,-2.044 -5.05866,-2.044 v 0 c -3.79467,0 -9.20934,1.416 -14.824,4.096"\n         style="fill:url(#radialGradient9-5);stroke:none;stroke-width:1.33333" /><path\n         id="path75"\n         d="m 136.77022,211.23252 c -1.29334,0.06 -2.288,1.17467 -2.22667,2.488 v 0 c 0.0133,0.296 0.0813,0.57734 0.19067,0.83334 v 0 c 0.37733,0.88533 1.25866,1.48133 2.25866,1.436 v 0 c 1.00267,-0.0493 1.82267,-0.728 2.11334,-1.64534 v 0 c 0.0853,-0.264 0.124,-0.55066 0.11333,-0.84666 v 0 c -0.0613,-1.276 -1.09733,-2.26934 -2.33733,-2.26934 v 0 c -0.04,0 -0.0773,0.003 -0.112,0.004"\n         style="fill:url(#radialGradient75-6);stroke:none;stroke-width:1.33333" /><path\n         id="path77"\n         d="m 169.01715,222.43187 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-9);stroke:none;stroke-width:1.33333" /><path\n         id="path79"\n         d="m 143.46814,214.20962 c -1.28934,0.0587 -2.28667,1.172 -2.22534,2.492 v 0 c 0.0133,0.28666 0.084,0.57466 0.19467,0.832 v 0 c 0.376,0.88666 1.25733,1.48266 2.25733,1.43466 v 0 c 1.00267,-0.0467 1.82134,-0.72533 2.11334,-1.64533 v 0 c 0.0827,-0.26133 0.124,-0.552 0.11066,-0.844 v 0 c -0.0653,-1.28133 -1.09866,-2.272 -2.34666,-2.272 v 0 c -0.032,0 -0.068,0.001 -0.104,0.003"\n         style="fill:url(#radialGradient79-1);stroke:none;stroke-width:1.33333" /><path\n         id="path83"\n         d="m 163.20773,222.78654 c -1.28666,0.0573 -2.28266,1.17466 -2.22266,2.48933 v 0 c 0.0147,0.28933 0.0787,0.57733 0.192,0.83467 v 0 c 0.37466,0.88266 1.25866,1.476 2.25733,1.42933 v 0 c 1,-0.044 1.824,-0.724 2.112,-1.64267 v 0 c 0.0907,-0.26133 0.124,-0.54666 0.112,-0.84266 v 0 c -0.0613,-1.276 -1.09467,-2.27067 -2.336,-2.27067 v 0 c -0.04,0 -0.0747,0.003 -0.11467,0.003"\n         style="fill:url(#radialGradient83-0);stroke:none;stroke-width:1.33333" /><path\n         id="path85"\n         d="m 153.86147,231.77495 c -1.28667,0.056 -2.284,1.17333 -2.22534,2.488 v 0 c 0.016,0.288 0.08,0.57733 0.19334,0.836 v 0 c 0.376,0.88 1.25733,1.47467 2.25866,1.428 v 0 c 0.99867,-0.0453 1.82267,-0.72267 2.112,-1.64133 v 0 c 0.088,-0.26267 0.12267,-0.54934 0.11067,-0.84267 v 0 c -0.0613,-1.27867 -1.09467,-2.272 -2.33733,-2.272 v 0 c -0.0387,0 -0.0747,0.004 -0.112,0.004"\n         style="fill:url(#radialGradient85-7);stroke:none;stroke-width:1.33333" /><path\n         id="path87"\n         d="m 156.14439,211.48323 c -1.28667,0.0573 -2.284,1.17467 -2.224,2.488 v 0 c 0.016,0.28933 0.08,0.57733 0.192,0.83467 v 0 c 0.376,0.88266 1.25866,1.476 2.25866,1.42933 v 0 c 1,-0.044 1.82267,-0.72267 2.112,-1.64267 v 0 c 0.0893,-0.26133 0.124,-0.54666 0.112,-0.84133 v 0 c -0.0613,-1.27733 -1.09466,-2.27067 -2.33733,-2.27067 v 0 c -0.0387,0 -0.0733,0.003 -0.11333,0.003"\n         style="fill:url(#radialGradient87-8);stroke:none;stroke-width:1.33333" /></g></g></svg>\n',
    "M": '<?xml version="1.0" encoding="UTF-8" standalone="no"?>\n<!-- Created with Inkscape (http://www.inkscape.org/) -->\n\n<svg\n   width="63.709248mm"\n   height="43.70472mm"\n   viewBox="0 0 63.709246 43.704718"\n   version="1.1"\n   id="svg1"\n   xml:space="preserve"\n   inkscape:version="1.4.3 (0d15f75, 2025-12-25)"\n   sodipodi:docname="M.svg"\n   xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape"\n   xmlns:sodipodi="http://sodipodi.sourceforge.net/DTD/sodipodi-0.dtd"\n   xmlns="http://www.w3.org/2000/svg"\n   xmlns:svg="http://www.w3.org/2000/svg"><sodipodi:namedview\n     id="namedview1"\n     pagecolor="#ffffff"\n     bordercolor="#000000"\n     borderopacity="0.25"\n     inkscape:showpageshadow="2"\n     inkscape:pageopacity="0.0"\n     inkscape:pagecheckerboard="0"\n     inkscape:deskcolor="#d1d1d1"\n     inkscape:document-units="mm"\n     inkscape:zoom="2.7545174"\n     inkscape:cx="36.66704"\n     inkscape:cy="92.5752"\n     inkscape:window-width="1920"\n     inkscape:window-height="991"\n     inkscape:window-x="-9"\n     inkscape:window-y="-9"\n     inkscape:window-maximized="1"\n     inkscape:current-layer="layer1" /><defs\n     id="defs1"><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(29.944918,0,0,31.072665,124.11605,234.76371)"\n       spreadMethod="pad"\n       id="radialGradient3"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.202018,0,0,18.88752,168.72075,213.58708)"\n       spreadMethod="pad"\n       id="radialGradient6"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop4" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(14.748164,-6.5909696,2.7657445,8.1230594,145.03468,226.82087)"\n       spreadMethod="pad"\n       id="radialGradient9"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop7" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.4817462,-2.9884689,1.97163,-4.688695,495.58499,36.967648)"\n       spreadMethod="pad"\n       id="radialGradient12"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop10" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop11" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop12" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath12"><path\n         d="m 492.82501,41.714304 h 5.516 v -9.495701 h -5.516 z"\n         id="path12" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath14"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path14" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath57"><path\n         d="M -21.860801,-2.0250001e-4 H 20.816201 V -62.223204 h -42.677002 z"\n         id="path57" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath58"><path\n         d="m 620.04667,389.236 h 60.288 V 473.2 h -60.288 z"\n         id="path58" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath61"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path61" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath62"><path\n         d="M -42.421399,0 H 41.153604 V -18.934 H -42.421399 Z"\n         id="path62" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath63"><path\n         d="M 596.15867,447.332 H 708.592 v 26.24533 H 596.15867 Z"\n         id="path63" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath66"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path66" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath67"><path\n         d="M -17.009698,-5.0000001e-4 H 20.600303 V -48.434501 h -37.610001 z"\n         id="path67" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath68"><path\n         d="M 626.80133,408.33467 H 677.448 V 473.772 h -50.64667 z"\n         id="path68" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath71"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path71" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,136.87935,213.61111)"\n       spreadMethod="pad"\n       id="radialGradient75"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop74" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop75" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,169.12798,224.81072)"\n       spreadMethod="pad"\n       id="radialGradient77"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,143.58212,216.58798)"\n       spreadMethod="pad"\n       id="radialGradient79"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop78" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop79" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3591334,0,0,2.3591334,148.14712,241.9406)"\n       spreadMethod="pad"\n       id="radialGradient81"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop80" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop81" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,163.32091,225.16177)"\n       spreadMethod="pad"\n       id="radialGradient83"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop82" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop83" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,153.9732,234.14904)"\n       spreadMethod="pad"\n       id="radialGradient85"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop84" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop85" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,156.257,213.85752)"\n       spreadMethod="pad"\n       id="radialGradient87"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop86" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop87" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.088001,0,0,-17.088001,457.0029,611.94415)"\n       spreadMethod="pad"\n       id="radialGradient345"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop343" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop344" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop345" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(10.386942,0,0,-10.386942,482.45648,623.58997)"\n       spreadMethod="pad"\n       id="radialGradient348"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop346" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop347" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop348" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.416008,3.6246166,1.578266,-4.4671698,468.94009,616.31226)"\n       spreadMethod="pad"\n       id="radialGradient351"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop349" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop350" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop351" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.888271,-2.191292,1.500142,-3.4379811,474.04248,618.56189)"\n       spreadMethod="pad"\n       id="radialGradient354"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop352" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop353" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop354" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath354"><path\n         d="m 471.94301,622.04202 h 4.197 v -6.962 h -4.197 z"\n         id="path354" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath356"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path356" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.129417,0,0,-17.681181,455.73108,563.27264)"\n       spreadMethod="pad"\n       id="radialGradient360"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop358" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop359" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop360" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.3039217,0,0,-5.4747686,453.74182,559.15039)"\n       spreadMethod="pad"\n       id="radialGradient363"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop361" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop362" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop363" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(9.7889442,0,0,-9.7889442,474.04147,535.25586)"\n       spreadMethod="pad"\n       id="radialGradient365"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop364" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop365" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(13.859965,2.4953344,2.3838298,-14.508271,529.60748,590.08429)"\n       spreadMethod="pad"\n       id="radialGradient398"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop396" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop397" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop398" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(16.999365,0,0,-17.794517,498.96738,577.06146)"\n       spreadMethod="pad"\n       id="radialGradient401"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop399" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop400" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop401" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.2615695,0,0,-5.5076814,497.80457,574.14221)"\n       spreadMethod="pad"\n       id="radialGradient404"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop402" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop403" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop404" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(3.8265915,0.6943604,0.6633328,-4.0055819,530.68677,590.08441)"\n       spreadMethod="pad"\n       id="radialGradient407"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop405" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop406" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop407" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.9543619,-2.4196851,1.5915818,-3.7604923,516.02515,583.75952)"\n       spreadMethod="pad"\n       id="radialGradient410"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop408" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop409" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop410" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath410"><path\n         d="m 513.83501,587.57602 h 4.379 v -7.636 h -4.379 z"\n         id="path410" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath412"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path412" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,469.0929,479.88107)"\n       spreadMethod="pad"\n       id="radialGradient419"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop417" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop418" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop419" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,467.53519,476.02194)"\n       spreadMethod="pad"\n       id="radialGradient422"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop420" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop421" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop422" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,507.55173,439.42313)"\n       spreadMethod="pad"\n       id="radialGradient427"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop425" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop426" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop427" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,505.99402,435.564)"\n       spreadMethod="pad"\n       id="radialGradient430"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop428" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop429" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop430" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(19.683309,-3.5213403,3.3854092,20.473634,151.17283,330.66865)"\n       spreadMethod="pad"\n       id="radialGradient436"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop434" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop435" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop436" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.141743,0,0,25.111086,107.65912,349.04613)"\n       spreadMethod="pad"\n       id="radialGradient439"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop437" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop438" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop439" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4722473,0,0,7.7722738,114.72071,492.57303)"\n       spreadMethod="pad"\n       id="radialGradient442"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop440" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop441" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop442" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.4343553,-0.97986038,0.94203584,5.6525553,161.41859,470.0759)"\n       spreadMethod="pad"\n       id="radialGradient445"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop443" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop444" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop445" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0816238,-2.5609398,1.6952205,-3.98002,520.40234,194.37172)"\n       spreadMethod="pad"\n       id="radialGradient448"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop446" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop447" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop448" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath448"><path\n         d="m 518.07001,198.41101 h 4.663 v -8.082 h -4.663 z"\n         id="path448" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath450"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path450" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.094729,0,0,24.870851,124.129,154.78121)"\n       spreadMethod="pad"\n       id="radialGradient456"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop454" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop455" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop456" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4606481,0,0,7.7009657,113.04409,52.125052)"\n       spreadMethod="pad"\n       id="radialGradient459"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop457" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop458" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop459" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.150868,0,0,-18.834442,451.0462,147.19133)"\n       spreadMethod="pad"\n       id="radialGradient462"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop460" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop461" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop462" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(11.033005,0,0,-11.448516,478.08298,160.02737)"\n       spreadMethod="pad"\n       id="radialGradient465"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop463" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop464" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop465" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.9394798,3.9950626,1.6764334,-4.9237266,463.72586,152.00581)"\n       spreadMethod="pad"\n       id="radialGradient468"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop466" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop467" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop468" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0057206,-2.4152482,1.5934501,-3.7893524,469.14563,154.48541)"\n       spreadMethod="pad"\n       id="radialGradient471"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop469" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop470" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop471" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath471"><path\n         d="m 466.91501,158.32201 h 4.458 v -7.675 h -4.458 z"\n         id="path471" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath473"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path473" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,160.06101,208.15524)"\n       spreadMethod="pad"\n       id="radialGradient77-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,142.82225,469.66275)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,144.99955,477.22245)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,152.21191,481.3324)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,131.30639,474.06771)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.75035,480.51167)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-38"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-7" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,138.83947,489.85995)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-19"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-86" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.16259,44.279802)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,126.77659,60.979642)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,121.51251,31.391882)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-0" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,92.47587,51.717192)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-7" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,108.08659,64.332832)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-84" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,128.59835,171.61115)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-5" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,122.17299,137.41595)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,102.38731,35.583232)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-4"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(29.944918,0,0,31.072665,124.11605,234.76371)"\n       spreadMethod="pad"\n       id="radialGradient3-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop1-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop2-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop3-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.202018,0,0,18.88752,168.72075,213.58708)"\n       spreadMethod="pad"\n       id="radialGradient6-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop4-7" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop5-3" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop6-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(14.748164,-6.5909696,2.7657445,8.1230594,145.03468,226.82087)"\n       spreadMethod="pad"\n       id="radialGradient9-5"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop7-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop8-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop9-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.4817462,-2.9884689,1.97163,-4.688695,495.58499,36.967648)"\n       spreadMethod="pad"\n       id="radialGradient12-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop10-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop11-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop12-6" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath12-6"><path\n         d="m 492.82501,41.714304 h 5.516 v -9.495701 h -5.516 z"\n         id="path12-0" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath14-3"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path14-8" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath57-0"><path\n         d="M -21.860801,-2.0250001e-4 H 20.816201 V -62.223204 h -42.677002 z"\n         id="path57-1" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath58-2"><path\n         d="m 620.04667,389.236 h 60.288 V 473.2 h -60.288 z"\n         id="path58-5" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath61-0"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path61-9" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath62-4"><path\n         d="M -42.421399,0 H 41.153604 V -18.934 H -42.421399 Z"\n         id="path62-7" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath63-8"><path\n         d="M 596.15867,447.332 H 708.592 v 26.24533 H 596.15867 Z"\n         id="path63-3" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath66-5"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path66-1" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath67-2"><path\n         d="M -17.009698,-5.0000001e-4 H 20.600303 V -48.434501 h -37.610001 z"\n         id="path67-0" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath68-1"><path\n         d="M 626.80133,408.33467 H 677.448 V 473.772 h -50.64667 z"\n         id="path68-6" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath71-4"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path71-0" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,136.87935,213.61111)"\n       spreadMethod="pad"\n       id="radialGradient75-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop74-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop75-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,169.12798,224.81072)"\n       spreadMethod="pad"\n       id="radialGradient77-9"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,143.58212,216.58798)"\n       spreadMethod="pad"\n       id="radialGradient79-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop78-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop79-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3591334,0,0,2.3591334,148.14712,241.9406)"\n       spreadMethod="pad"\n       id="radialGradient81-9"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop80-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop81-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,163.32091,225.16177)"\n       spreadMethod="pad"\n       id="radialGradient83-0"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop82-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop83-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,153.9732,234.14904)"\n       spreadMethod="pad"\n       id="radialGradient85-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop84-8" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop85-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3601119,0,0,2.3601119,156.257,213.85752)"\n       spreadMethod="pad"\n       id="radialGradient87-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop86-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop87-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.088001,0,0,-17.088001,457.0029,611.94415)"\n       spreadMethod="pad"\n       id="radialGradient345-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop343-0" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop344-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop345-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(10.386942,0,0,-10.386942,482.45648,623.58997)"\n       spreadMethod="pad"\n       id="radialGradient348-4"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop346-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop347-6" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop348-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.416008,3.6246166,1.578266,-4.4671698,468.94009,616.31226)"\n       spreadMethod="pad"\n       id="radialGradient351-1"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop349-0" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop350-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop351-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.888271,-2.191292,1.500142,-3.4379811,474.04248,618.56189)"\n       spreadMethod="pad"\n       id="radialGradient354-6"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop352-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop353-3" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop354-4" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath354-8"><path\n         d="m 471.94301,622.04202 h 4.197 v -6.962 h -4.197 z"\n         id="path354-4" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath356-9"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path356-9" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(17.129417,0,0,-17.681181,455.73108,563.27264)"\n       spreadMethod="pad"\n       id="radialGradient360-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop358-5" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop359-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop360-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.3039217,0,0,-5.4747686,453.74182,559.15039)"\n       spreadMethod="pad"\n       id="radialGradient363-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop361-3" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop362-7" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop363-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(9.7889442,0,0,-9.7889442,474.04147,535.25586)"\n       spreadMethod="pad"\n       id="radialGradient365-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop364-8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop365-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(13.859965,2.4953344,2.3838298,-14.508271,529.60748,590.08429)"\n       spreadMethod="pad"\n       id="radialGradient398-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop396-8" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop397-0" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop398-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(16.999365,0,0,-17.794517,498.96738,577.06146)"\n       spreadMethod="pad"\n       id="radialGradient401-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop399-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop400-9" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop401-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.2615695,0,0,-5.5076814,497.80457,574.14221)"\n       spreadMethod="pad"\n       id="radialGradient404-9"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop402-7" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop403-2" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop404-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(3.8265915,0.6943604,0.6633328,-4.0055819,530.68677,590.08441)"\n       spreadMethod="pad"\n       id="radialGradient407-8"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop405-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop406-8" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop407-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(1.9543619,-2.4196851,1.5915818,-3.7604923,516.02515,583.75952)"\n       spreadMethod="pad"\n       id="radialGradient410-0"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop408-7" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop409-8" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop410-1" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath410-5"><path\n         d="m 513.83501,587.57602 h 4.379 v -7.636 h -4.379 z"\n         id="path410-8" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath412-6"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path412-1" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,469.0929,479.88107)"\n       spreadMethod="pad"\n       id="radialGradient419-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop417-4" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop418-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop419-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,467.53519,476.02194)"\n       spreadMethod="pad"\n       id="radialGradient422-8"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop420-6" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop421-2" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop422-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(22.803827,0,0,-23.53837,507.55173,439.42313)"\n       spreadMethod="pad"\n       id="radialGradient427-5"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop425-3" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop426-9" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop427-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.060936,0,0,-7.2883782,505.99402,435.564)"\n       spreadMethod="pad"\n       id="radialGradient430-4"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop428-6" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop429-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop430-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(19.683309,-3.5213403,3.3854092,20.473634,151.17283,330.66865)"\n       spreadMethod="pad"\n       id="radialGradient436-2"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop434-1" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop435-1" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop436-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.141743,0,0,25.111086,107.65912,349.04613)"\n       spreadMethod="pad"\n       id="radialGradient439-7"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop437-6" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop438-2" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop439-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4722473,0,0,7.7722738,114.72071,492.57303)"\n       spreadMethod="pad"\n       id="radialGradient442-5"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop440-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop441-0" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop442-0" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(5.4343553,-0.97986038,0.94203584,5.6525553,161.41859,470.0759)"\n       spreadMethod="pad"\n       id="radialGradient445-3"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop443-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop444-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop445-8" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0816238,-2.5609398,1.6952205,-3.98002,520.40234,194.37172)"\n       spreadMethod="pad"\n       id="radialGradient448-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop446-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop447-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop448-3" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath448-2"><path\n         d="m 518.07001,198.41101 h 4.663 v -8.082 h -4.663 z"\n         id="path448-5" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath450-2"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path450-5" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(24.094729,0,0,24.870851,124.129,154.78121)"\n       spreadMethod="pad"\n       id="radialGradient456-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop454-6" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop455-7" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop456-7" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(7.4606481,0,0,7.7009657,113.04409,52.125052)"\n       spreadMethod="pad"\n       id="radialGradient459-2"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop457-2" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop458-9" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop459-4" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(18.150868,0,0,-18.834442,451.0462,147.19133)"\n       spreadMethod="pad"\n       id="radialGradient462-1"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop460-9" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop461-6" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop462-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(11.033005,0,0,-11.448516,478.08298,160.02737)"\n       spreadMethod="pad"\n       id="radialGradient465-8"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop463-2" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop464-5" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop465-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(8.9394798,3.9950626,1.6764334,-4.9237266,463.72586,152.00581)"\n       spreadMethod="pad"\n       id="radialGradient468-4"><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0"\n         id="stop466-9" /><stop\n         style="stop-opacity:1;stop-color:#325cb2"\n         offset="0.12030684"\n         id="stop467-1" /><stop\n         style="stop-opacity:1;stop-color:#3293d7"\n         offset="1"\n         id="stop468-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.0057206,-2.4152482,1.5934501,-3.7893524,469.14563,154.48541)"\n       spreadMethod="pad"\n       id="radialGradient471-5"><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0"\n         id="stop469-0" /><stop\n         style="stop-opacity:1;stop-color:#f68712"\n         offset="0.12030684"\n         id="stop470-8" /><stop\n         style="stop-opacity:1;stop-color:#f8c760"\n         offset="1"\n         id="stop471-3" /></radialGradient><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath471-9"><path\n         d="m 466.91501,158.32201 h 4.458 v -7.675 h -4.458 z"\n         id="path471-3" /></clipPath><clipPath\n       clipPathUnits="userSpaceOnUse"\n       id="clipPath473-9"><path\n         d="M 0,0 H 958.81467 V 874.66667 H 0 Z"\n         id="path473-6" /></clipPath><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,160.06101,208.15524)"\n       spreadMethod="pad"\n       id="radialGradient77-2-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,142.82225,469.66275)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-7"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-9" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,144.99955,477.22245)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-35"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-76" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-65" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,152.21191,481.3324)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-1-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-8-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-9-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,131.30639,474.06771)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-3-4"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-1-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-2-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.75035,480.51167)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-38-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-7-1" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-4-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,138.83947,489.85995)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-19-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-86-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-5-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,137.16259,44.279802)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-4-3" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-8-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,126.77659,60.979642)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-2" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-5" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,121.51251,31.391882)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-6-3"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-0-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-6-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,92.47587,51.717192)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-3-8"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-7-6" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-8-2" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,108.08659,64.332832)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-5-1"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-9-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-84-6" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,128.59835,171.61115)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-1-2"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-5-9" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-4-1" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,122.17299,137.41595)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-2-5"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-6-0" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-0-3" /></radialGradient><radialGradient\n       fx="0"\n       fy="0"\n       cx="0"\n       cy="0"\n       r="1"\n       gradientUnits="userSpaceOnUse"\n       gradientTransform="matrix(2.3604466,0,0,2.3604466,102.38731,35.583232)"\n       spreadMethod="pad"\n       id="radialGradient77-2-5-2-0-4-6"><stop\n         style="stop-opacity:1;stop-color:#bf1b30"\n         offset="0"\n         id="stop76-4-2-2-6-1-4" /><stop\n         style="stop-opacity:1;stop-color:#d84e1e"\n         offset="1"\n         id="stop77-5-7-1-1-2-9" /></radialGradient></defs><g\n     inkscape:label="Layer 1"\n     inkscape:groupmode="layer"\n     id="layer1"\n     transform="translate(-79.436934,-180.77426)"><g\n       inkscape:label="Layer 1"\n       id="layer1-8"\n       transform="translate(39.042624,52.93468)" /><g\n       inkscape:label="Layer 1"\n       id="layer1-2"\n       transform="translate(-163.1037,62.879377)" /><g\n       id="g3"\n       transform="matrix(0.65746005,0,0,0.68830705,22.099398,-128.24301)"><path\n         id="path436"\n         d="m 155.41539,453.17536 c -11.36788,3.87855 -18.21509,13.97763 -15.29821,22.55744 v 0 c 2.92846,8.57982 14.51111,12.39528 25.87512,8.52044 v 0 c 11.36787,-3.87484 18.21766,-13.97391 15.29563,-22.55497 v 0 c -2.14522,-6.30712 -8.968,-10.03721 -16.94313,-10.03721 v 0 c -2.88216,0 -5.91351,0.48621 -8.92941,1.5143"\n         style="fill:none;stroke:#006200;stroke-width:1.2614;stroke-opacity:1" /><path\n         id="path77-7-4"\n         d="m 142.71142,467.2839 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-7);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5"\n         d="m 144.88872,474.8436 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-35);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-9"\n         d="m 152.10108,478.95355 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-1-8);stroke:none;stroke-width:1.33333" /><path\n         id="path439"\n         d="m 113.13768,467.33698 c -14.862666,2.46533 -25.465332,13.916 -23.678666,25.57866 v 0 c 1.794667,11.664 15.289336,19.11733 30.147996,16.65467 v 0 c 14.86267,-2.46134 25.464,-13.91333 23.67734,-25.57867 v 0 c -1.55334,-10.16266 -11.99334,-17.12666 -24.484,-17.12533 v 0 c -1.84934,0 -3.74667,0.15333 -5.66267,0.47067"\n         style="fill:none;stroke:#006200;stroke-width:1.33333;stroke-opacity:1" /><path\n         id="path442"\n         d="m 114.44568,486.46898 c -4.852,0.23466 -8.65867,3.156 -8.49867,6.532 v 0 c 0.0307,0.764 0.264,1.484 0.65867,2.14133 v 0 c 1.36,2.24267 4.63333,3.72 8.38667,3.536 v 0 c 3.76,-0.18533 6.88666,-1.98 8.03066,-4.35067 v 0 c 0.34134,-0.68266 0.50267,-1.42666 0.472,-2.184 v 0 c -0.15466,-3.22933 -3.86266,-5.69066 -8.41466,-5.69066 v 0 c -0.21067,0 -0.42134,0.005 -0.63467,0.016"\n         style="fill:url(#radialGradient442-5);stroke:none;stroke-width:1.33333" /><path\n         id="path445"\n         d="m 160.47768,465.67164 c -3.50133,0.808 -5.91333,3.432 -5.38667,5.868 v 0 c 0.11467,0.548 0.36934,1.04267 0.73467,1.46667 v 0 c 1.26667,1.45467 3.82533,2.10267 6.52933,1.476 v 0 c 2.71467,-0.62933 4.772,-2.34267 5.31467,-4.21867 v 0 c 0.16667,-0.54133 0.19067,-1.10266 0.0773,-1.65066 v 0 c -0.42,-1.94934 -2.588,-3.18134 -5.224,-3.18134 v 0 c -0.66,0 -1.348,0.0773 -2.04533,0.24"\n         style="fill:url(#radialGradient445-3);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-4"\n         d="m 131.19556,471.68886 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-3-4);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-7"\n         d="m 137.63952,478.13282 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-38-6);stroke:none;stroke-width:1.33333" /><path\n         id="path77-7-4-5-8"\n         d="m 138.72864,487.4811 c -1.28801,0.0587 -2.28667,1.17334 -2.22534,2.48934 v 0 c 0.0133,0.29333 0.08,0.572 0.19333,0.83333 v 0 c 0.376,0.884 1.25467,1.48133 2.25734,1.436 v 0 c 1.00133,-0.0493 1.82267,-0.72667 2.11333,-1.644 v 0 c 0.084,-0.264 0.124,-0.552 0.11067,-0.84667 v 0 c -0.06,-1.27733 -1.096,-2.27066 -2.33867,-2.27066 v 0 c -0.0373,0 -0.072,0 -0.11066,0.003"\n         style="fill:url(#radialGradient77-2-5-2-19-3);stroke:none;stroke-width:1.33333" /></g></g></svg>\n',
}


def _shade(hex_color: str, factor: float) -> str:
    """Lighten/darken a hex color by multiplying RGB values and clipping."""
    h = str(hex_color).lstrip("#")
    try:
        r, g, b = int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16)
    except Exception:
        return hex_color
    r = max(0, min(255, int(r * factor)))
    g = max(0, min(255, int(g * factor)))
    b = max(0, min(255, int(b * factor)))
    return "#{:02X}{:02X}{:02X}".format(r, g, b)


def colorize_stage_svg(svg_text: str, part_colors: dict[str, str]) -> str:
    """Recolor the uploaded SVG schematic using MorphoNet component rules."""
    cell_wall = part_colors.get("cell_wall", "#BDBDBD")
    actin = part_colors.get("actin", "#BDBDBD")
    nucleus = part_colors.get("nucleus", "#BDBDBD")

    replacements = {
        "#006200": cell_wall,
        "#bf1b30": _shade(actin, 0.85),
        "#d84e1e": actin,
        "#f68712": _shade(actin, 1.15),
        "#f8c760": _shade(actin, 1.35),
        "#325cb2": nucleus,
        "#3293d7": _shade(nucleus, 1.25),
    }
    out = svg_text
    for old, new in replacements.items():
        out = out.replace(old, new).replace(old.upper(), new)
    if "<svg" in out and "style=" not in out.split("<svg", 1)[1].split(">", 1)[0]:
        out = out.replace("<svg", '<svg style="width:100%; max-width:270px; height:auto; display:block; margin:auto;"', 1)
    return out


def render_stage_svg_html(stage: str, part_colors: dict[str, str]) -> str:
    """Return HTML containing the uploaded SVG schematic without parameter labels."""
    svg = colorize_stage_svg(SVG_TEMPLATES.get(stage, ""), part_colors)
    return f"""
    <div style="border:1px solid rgba(0,0,0,0.12); border-radius:16px; padding:10px; background:#fff;">
      <div style="text-align:center; font-weight:700; font-size:18px; margin-bottom:6px;">{html_escape(stage)}</div>
      <div style="position:relative; min-height:205px; display:flex; align-items:center; justify-content:center;">
        {svg}
      </div>
    </div>
    """


                               
                                                          
                               


def _first_existing_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
    lower_map = {str(c).strip().casefold(): c for c in df.columns}
    for c in candidates:
        key = str(c).strip().casefold()
        if key in lower_map:
            return lower_map[key]
    return None


def _selected_functional_group_gene_table(
    nodes_df: pd.DataFrame,
    info_df: pd.DataFrame,
    selected_groups: set[str],
) -> pd.DataFrame:
    columns = ["Functional group", "Systematic name", "Standard name", "SGD ID", "Name description"]
    if not selected_groups or "GO" not in nodes_df.columns:
        return pd.DataFrame(columns=columns)

    tmp_nodes = nodes_df.copy()
    tmp_nodes["GO"] = tmp_nodes["GO"].astype(str).str.strip()
    selected = {str(g).strip() for g in selected_groups}
    tmp_nodes = tmp_nodes[tmp_nodes["GO"].isin(selected)]
    if tmp_nodes.empty:
        return pd.DataFrame(columns=columns)

    std_col = _first_existing_column(info_df, ["Standard name", "Gene.name", "Gene name", "standard_name"])
    sgd_col = _first_existing_column(info_df, ["SGD ID", "SGD.ID", "SGD_ID", "sgd_id"])
    desc_col = _first_existing_column(
        info_df,
        [
            "Name description",
            "Name.Description",
            "Name.description",
            "Description",
            "Gene description",
            "Gene.Description",
            "Feature.Description",
            "Feature description",
        ],
    )

    rows = []
    for sys_name, row in tmp_nodes.iterrows():
        std = ""
        sgd = ""
        desc = ""
        if sys_name in info_df.index:
            if std_col is not None and pd.notna(info_df.loc[sys_name, std_col]):
                std = str(info_df.loc[sys_name, std_col])
            if sgd_col is not None and pd.notna(info_df.loc[sys_name, sgd_col]):
                sgd = str(info_df.loc[sys_name, sgd_col])
            if desc_col is not None and pd.notna(info_df.loc[sys_name, desc_col]):
                desc = str(info_df.loc[sys_name, desc_col])
        rows.append(
            {
                "Functional group": str(row.get("GO", "")),
                "Systematic name": str(sys_name),
                "Standard name": std,
                "SGD ID": sgd,
                "Name description": desc,
            }
        )

    out = pd.DataFrame(rows, columns=columns)
    return out.sort_values(["Functional group", "Systematic name"], kind="stable").reset_index(drop=True)


def render_functional_group_gene_export(
    nodes_df: pd.DataFrame,
    info_df: pd.DataFrame,
    selected_groups: set[str],
    prefix: str,
) -> None:
    if not selected_groups:
        return

    export_df = _selected_functional_group_gene_table(nodes_df, info_df, selected_groups)
    if export_df.empty:
        st.warning("No genes were found for the selected functional group(s).")
        return

    st.caption(f"Selected functional group gene list: {len(export_df):,} genes.")
    try:
        table_export_box = st.popover("Export table ▾", use_container_width=False)
    except Exception:
        table_export_box = st.expander("Export table", expanded=False)

    with table_export_box:
        st.download_button(
            "CSV",
            data=export_df.to_csv(index=False).encode("utf-8"),
            file_name=f"{prefix}_selected_functional_groups_genes.csv",
            mime="text/csv",
            key=f"{prefix}_selected_functional_groups_csv",
            use_container_width=True,
        )
        st.download_button(
            "TSV",
            data=export_df.to_csv(index=False, sep="\t").encode("utf-8"),
            file_name=f"{prefix}_selected_functional_groups_genes.tsv",
            mime="text/tab-separated-values",
            key=f"{prefix}_selected_functional_groups_tsv",
            use_container_width=True,
        )

with tab1:

    global_raw, global_mapped, global_dups = _panel_genes("tab1_shared", info_df)
    st.session_state["global_genes_raw"] = global_raw
    st.session_state["global_genes_mapped"] = global_mapped
    st.session_state["global_genes_dups"] = global_dups

    if global_raw:
        matching_ess = [g for g in global_mapped if g in nodes_df_Essen.index]
        matching_non = [g for g in global_mapped if g in nodes_df_NonEssen.index]
        st.info(
            f"Number of input genes: {len(global_raw)} · Found in essential: {len(matching_ess)} · Found in nonessential: {len(matching_non)}"
        )
        if global_dups:
            st.warning("Duplicate gene names detected in the input.")
    else:
        matching_ess = []
        matching_non = []

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Essential Genes Network")
        st.caption(f"Essential network contains {len(nodes_df_Essen):,} genes.")

                                                           
        go_options_ess = ESS_FUNCTIONAL_GROUPS.copy()
        go_colors_ess: dict[str, str] = {}

                                                                             
        if "GO" in nodes_df_Essen.columns and "Col" in nodes_df_Essen.columns:
            tmp = nodes_df_Essen[["GO", "Col"]].dropna()
            tmp["GO"] = tmp["GO"].astype(str).str.strip()
            tmp["Col"] = tmp["Col"].astype(str).str.strip()
            tmp = tmp[tmp["GO"].isin(go_options_ess)]
            for _, r in tmp.iterrows():
                go_colors_ess[str(r["GO"])] = str(r["Col"])

        colored_groups_ess: set[str] = set()
        if go_options_ess:
            try:
                go_box_ess = st.popover("Select functional groups ...")
            except Exception:
                go_box_ess = st.expander("Select functional groups ...", expanded=False)

            with go_box_ess:
                st.caption("Select a group to show on the map.")

                state_key = "essen_go_color_checks"
                if state_key not in st.session_state:
                    st.session_state[state_key] = {g: False for g in go_options_ess}

                checks: dict = st.session_state[state_key]
                                                
                for g in list(checks.keys()):
                    if g not in go_options_ess:
                        checks.pop(g, None)
                for g in go_options_ess:
                    checks.setdefault(g, False)

                c1, c2 = st.columns(2)
                with c1:
                    if st.button("Select all", key="essen_go_select_all", use_container_width=True):
                        for g in go_options_ess:
                            checks[g] = True
                        st.session_state[state_key] = checks
                        st.rerun()
                with c2:
                    if st.button("Clear all", key="essen_go_clear_all", use_container_width=True):
                        for g in go_options_ess:
                            checks[g] = False
                        st.session_state[state_key] = checks
                        st.rerun()

                st.divider()

                for g in go_options_ess:
                    checks[g] = st.checkbox(g, value=bool(checks.get(g, False)), key=f"essen_go_chk_{g}")

                colored_groups_ess = {g for g, v in checks.items() if v}

        render_functional_group_gene_export(
            nodes_df_Essen,
            info_df,
            colored_groups_ess,
            "essential",
        )

        show_edges_essen = st.checkbox('Show edges', value=False, key='show_edges_essen')

        graph = build_sigma_graph(
            nodes_df_Essen,
            info_df,
            highlighted=matching_ess,
            edges_df=edges_df_Essen,
            show_edges_precomputed=show_edges_essen,
            visible_groups=None,
            group_colors=go_colors_ess if go_colors_ess else None,
            group_col="GO",
            color_col="Col",
            display_name_col=None,
            color_by_group=True if go_options_ess else False,
            colored_groups=colored_groups_ess if go_options_ess else None,
            default_node_color="#00000050",
            size_scale=2.0,
        )
        sigma_viewer_component(graph, height=620, show_threshold_slider=False, title="")

    with col2:
        st.subheader("NonEssential Genes Network")
        st.caption(f"NonEssential network contains {len(nodes_df_NonEssen):,} genes.")

        go_options = NONESS_FUNCTIONAL_GROUPS.copy()
        go_colors: dict[str, str] = {}

        if "GO" in nodes_df_NonEssen.columns and "Col" in nodes_df_NonEssen.columns:
            tmp = nodes_df_NonEssen[["GO", "Col"]].dropna()
            tmp["GO"] = tmp["GO"].astype(str).str.strip()
            tmp["Col"] = tmp["Col"].astype(str).str.strip()
            tmp = tmp[tmp["GO"].isin(set(go_options))]
                                            
            for grp, sub in tmp.groupby("GO", sort=False):
                col = sub["Col"].dropna().astype(str).iloc[0] if len(sub) else None
                if col and str(col).lower() != "nan":
                    go_colors[str(grp)] = str(col)
                                                    
                                                                       
        colored_groups: set[str] = set()
        if go_options:
            try:
                go_box = st.popover("Select functional groups ...")
            except Exception:
                go_box = st.expander("Select functional groups ...", expanded=False)

            with go_box:
                st.caption("Select a group to shown on the map.")

                                                           
                state_key = "nonessen_go_color_checks"
                if state_key not in st.session_state:
                    st.session_state[state_key] = {g: False for g in go_options}

                                                                                 
                checks: dict = st.session_state[state_key]
                for g in list(checks.keys()):
                    if g not in go_options:
                        checks.pop(g, None)
                for g in go_options:
                    checks.setdefault(g, False)

                c1, c2 = st.columns(2)
                with c1:
                    if st.button("Select all", key="nonessen_go_select_all", use_container_width=True):
                        for g in go_options:
                            checks[g] = True
                        st.session_state[state_key] = checks
                        st.rerun()
                with c2:
                    if st.button("Clear all", key="nonessen_go_clear_all", use_container_width=True):
                        for g in go_options:
                            checks[g] = False
                        st.session_state[state_key] = checks
                        st.rerun()

                st.divider()

                                        
                for g in go_options:
                    checks[g] = st.checkbox(g, value=bool(checks.get(g, False)), key=f"nonessen_go_chk_{g}")

                colored_groups = {g for g, v in checks.items() if v}

        render_functional_group_gene_export(
            nodes_df_NonEssen,
            info_df,
            colored_groups,
            "nonessential",
        )

        show_edges_nonessen = st.checkbox('Show edges', value=False, key='show_edges_nonessen')

        graph = build_sigma_graph(
            nodes_df_NonEssen,
            info_df,
            highlighted=matching_non,
            edges_df=edges_df_NonEssen,
            show_edges_precomputed=show_edges_nonessen,
            visible_groups=None,
            group_colors=go_colors if go_colors else None,
            group_col="GO",
            color_col="Col",
            display_name_col=None,
            color_by_group=True if go_options else False,
            colored_groups=colored_groups if go_options else None,
            default_node_color="#00000050",
        )
        sigma_viewer_component(graph, height=620, show_threshold_slider=False, title="")


                               
                                         
                               
with tab2:
    st.header("Morphological Similarity Networkss")
    st.info("Use the search box under Morphology–Gene Function Networks tab.")
    st.info("Pearson correlation coefficients (r) were computed between morphological profiles based on CalMorph-derived Z values after projection onto wild-type (WT) principal component axes capturing 99% of the cumulative explained variance. Before PCA, morphological parameters with zero variance among WT samples were excluded, and the remaining parameters were centered and scaled")

    col1, col2 = st.columns(2)
                                                
    def render_similarity(
        title: str,
        genes_input: List[str],
        corr_df: pd.DataFrame,
        study_label: str,
        layout_seed: int,
        filename_prefix: str,
    ):
        st.subheader(title)

                                                         
        genes = [
            g for g in genes_input
            if g in info_df.index and str(info_df.loc[g, "Study"]) == study_label
        ]
        st.info(f"Input genes: {len(genes_input)} · In {study_label}: {len(genes)}")

        if len(genes) < 2:
            st.info("Provide at least 2 genes (in this study) to build a similarity network.")
            return

                                                                                                        
                                                             
        correlations_all = extract_correlations(genes, corr_df)

        r_rows = []
        for (g1, g2), r in correlations_all.items():
            std1 = None
            std2 = None
            if g1 in info_df.index and "Standard name" in info_df.columns:
                v = info_df.loc[g1, "Standard name"]
                std1 = None if pd.isna(v) else str(v)
            if g2 in info_df.index and "Standard name" in info_df.columns:
                v = info_df.loc[g2, "Standard name"]
                std2 = None if pd.isna(v) else str(v)
            r_rows.append(
                {
                    "Gene 1 (systematic)": g1,
                    "Gene 1 (standard)": std1,
                    "Gene 2 (systematic)": g2,
                    "Gene 2 (standard)": std2,
                    "r": float(r),
                }
            )
        r_table_df = pd.DataFrame(r_rows)
                                        
        r_threshold = st.slider(
            f"|r| threshold ({study_label})",
            min_value=0.0,
            max_value=1.0,
            value=0.5,
            step=0.01,
            key=f"thr_{study_label}",
        )

        show_r_values = st.checkbox(
            "Show r values",
            value=False,
            key=f"show_r_{study_label}",
        )

                                          
        correlations = extract_correlations_threshold(genes, corr_df, r_threshold)

                       
        layout = generate_layout(genes, correlations, seed=layout_seed)

              
        fig, _ = plot_similarity_network(genes, correlations, info_df, layout=layout, show_r_values=show_r_values)
        

                                                       
        try:
            export_box = st.popover("Export", use_container_width=False)
        except Exception:
                                                   
            export_box = st.expander("Export", expanded=False)

        with export_box:
            svg_buf = save_svg(fig)
            st.download_button(
                label="SVG",
                data=svg_buf,
                file_name=f"{filename_prefix}_Similarity_Network.svg",
                mime="image/svg+xml",
                key=f"dl_svg_{study_label}",
                use_container_width=True,
            )

                                                                       
            fig_tiff, _ = plot_similarity_network(genes, correlations, info_df, layout=layout, show_r_values=show_r_values)
            tiff_buf = save_tiff(fig_tiff)
            st.download_button(
                label="TIFF",
                data=tiff_buf,
                file_name=f"{filename_prefix}_Similarity_Network.tiff",
                mime="image/tiff",
                key=f"dl_tiff_{study_label}",
                use_container_width=True,
            )

            json_bytes = similarity_network_json_bytes(genes, correlations, layout, info_df)
            st.download_button(
                label="JSON",
                data=json_bytes,
                file_name=f"{filename_prefix}_Similarity_Network.json",
                mime="application/json",
                key=f"dl_json_{study_label}",
                use_container_width=True,
            )

                                                                                                               
            if r_table_df is not None and len(r_table_df) > 0:
                st.download_button(
                    label="Table (CSV)",
                    data=convert_df_to_file(r_table_df, "csv"),
                    file_name=f"{filename_prefix}_Similarity_r_values.csv",
                    mime="text/csv",
                    key=f"dl_rtable_csv_{study_label}",
                    use_container_width=True,
                )
                st.download_button(
                    label="Table (TSV)",
                    data=convert_df_to_file(r_table_df, "tsv"),
                    file_name=f"{filename_prefix}_Similarity_r_values.tsv",
                    mime="text/tab-separated-values",
                    key=f"dl_rtable_tsv_{study_label}",
                    use_container_width=True,
                )

        st.pyplot(fig, use_container_width=True)        


    with col1:
        st.subheader("Essential Genes Similarity Network")
        mapped_unique = st.session_state.get("global_genes_mapped", [])
        dups = st.session_state.get("global_genes_dups", [])
        if dups:
            st.warning("Duplicate gene names detected in the input.")
        render_similarity(
            title="",
            genes_input=mapped_unique,
            corr_df=Corr_df_Essen,
            study_label="Essential",
            layout_seed=42,
            filename_prefix="Essential",
        )

    with col2:
        st.subheader("NonEssential Genes Similarity Network")
        mapped_unique = st.session_state.get("global_genes_mapped", [])
        dups = st.session_state.get("global_genes_dups", [])
        if dups:
            st.warning("Duplicate gene names detected in the input.")
        render_similarity(
            title="",
            genes_input=mapped_unique,
            corr_df=Corr_df_NonEssen,
            study_label="NonEssential",
            layout_seed=84,
            filename_prefix="Nonessential",
        )

                               
                             
                               
with tab3:
    st.header("Statistical Profiles")

                                                
    st.info("Use the search box under Morphology–Gene Function Networks tab.")
    raw = st.session_state.get("global_genes_raw", [])
    mapped = st.session_state.get("global_genes_mapped", [])
    dups = st.session_state.get("global_genes_dups", [])

    st.info(f"Number of input genes: {len(raw)} · Found: {len(mapped)}")
    if dups:
        st.warning("There are duplicate gene names in the input")

    if mapped:
        st.markdown("### Effect size (Z values)")
        display_values_table_only(z_values_df, "Z_values")

        st.markdown("### Statistical significance (FDR-adjusted)")
        display_values_table_only(q_values_df, "q_values")




                               
                              
                               
with tab4:
    st.header("Morphological Profile")
    st.write("Search by CalMorph parameter to find affected mutants, or search by mutant to retrieve its morphological defects.")

                                 
    st.markdown("### Mutant finder")
    st.write("Enter a CalMorph parameter(s) and a q-value threshold to find mutants with significant defects.")

    mf_c1, mf_c2 = st.columns([2, 1])
    with mf_c1:
        st.markdown(
            "Parameter name(s)<br>"
            "<span style='font-size: 0.9em; color: #666;'>(comma separated)</span>",
            unsafe_allow_html=True,
        )

        mutant_finder_parameter_raw = st.text_input(
            "Parameter name(s)",
            placeholder="e.g., DCV17.1_C, C11-1_A",
            key="mutant_finder_parameter_name",
            label_visibility="collapsed",
        ).strip()
    with mf_c2:
        st.markdown(
            "q-value threshold<br>"
            "<span style='font-size: 0.9em; color: transparent;'>( )</span>",
            unsafe_allow_html=True,
        )

        mutant_finder_q_threshold_txt = st.text_input(
            "q-value threshold",
            value="0.05",
            key="mutant_finder_q_threshold",
            label_visibility="collapsed",
        ).strip()
        try:
            mutant_finder_q_threshold = float(mutant_finder_q_threshold_txt)
            if mutant_finder_q_threshold < 0 or mutant_finder_q_threshold > 1:
                st.warning("Please enter a q-value threshold between 0 and 1 for Mutant finder.")
                mutant_finder_q_threshold = 0.05
        except ValueError:
            st.warning("Please enter a numeric q-value threshold for Mutant finder. Using 0.05.")
            mutant_finder_q_threshold = 0.05

    if mutant_finder_parameter_raw:
        mutant_finder_parameters = parse_calmorph_parameter_queries(mutant_finder_parameter_raw)
        mutant_finder_df, matched_parameters, missing_parameters = find_mutants_for_calmorph_parameters(
            mutant_finder_parameters,
            mutant_finder_q_threshold,
            q_values_df,
            z_values_df,
            info_df,
        )

        if missing_parameters:
            st.warning(
                "The following parameter(s) were not found: "
                + ", ".join(missing_parameters)
            )

        if not matched_parameters:
            st.warning("No entered parameter was found.")
        elif mutant_finder_df.empty:
            st.info(
                f"No mutants passed q ≤ {format_q_value_sci2(mutant_finder_q_threshold)} for the matched parameter(s): "
                + ", ".join(display_calmorph_parameter_name(p) for p in matched_parameters)
                + "."
            )
        else:
            mutant_finder_display_df = format_q_values_for_display(mutant_finder_df)
            matched_display = ", ".join(display_calmorph_parameter_name(p) for p in matched_parameters)
            st.success(
                f"Found {len(mutant_finder_df):,} mutant(s) significant for all matched parameter(s)"
            )
            try:
                mf_export_box = st.popover("Export", use_container_width=False)
            except Exception:
                mf_export_box = st.expander("Export", expanded=False)
            with mf_export_box:
                if len(matched_parameters) == 1:
                    safe_parameter = display_calmorph_parameter_name(matched_parameters[0]).replace("/", "_").replace(" ", "_")
                else:
                    safe_parameter = "Shared_Mutants_Multiple_CalMorph_Parameters"
                st.download_button(
                    label="TSV",
                    data=mutant_finder_display_df.to_csv(sep="\t", index=False).encode("utf-8"),
                    file_name=f"{safe_parameter}_Mutant_Finder.tsv",
                    mime="text/tab-separated-values",
                    key="download_mutant_finder_table_tsv",
                    use_container_width=True,
                )
                st.download_button(
                    label="CSV",
                    data=mutant_finder_display_df.to_csv(index=False).encode("utf-8"),
                    file_name=f"{safe_parameter}_Mutant_Finder.csv",
                    mime="text/csv",
                    key="download_mutant_finder_table_csv",
                    use_container_width=True,
                )
            st.dataframe(mutant_finder_display_df, hide_index=True, use_container_width=True)

    st.divider()

                                                     
    st.markdown("### Defect finder")
    st.write("Enter a gene name to retrieve morphological defect(s), if detected.")

                                                                 
    c1, c2 = st.columns([2, 1])
    with c1:
        profile_gene_raw = st.text_input(
            "Gene name",
            placeholder="e.g., CDC25",
            key="morph_profile_gene",
        ).upper().strip()
    with c2:
        q_threshold_txt = st.text_input(
            "q-value threshold",
            value="0.05",
            key="morph_profile_q_threshold",
        ).strip()
        try:
            q_threshold = float(q_threshold_txt)
            if q_threshold < 0 or q_threshold > 1:
                st.warning("Please enter a q-value threshold between 0 and 1.")
                q_threshold = 0.05
        except ValueError:
            st.warning("Please enter a numeric q-value threshold. Using 0.05.")
            q_threshold = 0.05

    profile_df = pd.DataFrame(columns=["Parameter", "Z value", "q value", "Description"])
    mapped_gene = None

    if profile_gene_raw:
        profile_gene_query, has_multiple_profile_genes = parse_single_gene_query(profile_gene_raw)
        if has_multiple_profile_genes:
            st.warning("Please enter only one gene name. Do not include multiple genes.")
        else:
            mapped, dups = map_to_systematic_name([profile_gene_query], info_df)
            if not mapped:
                st.warning(f"Gene '{profile_gene_query}' was not found.")
            else:
                mapped_gene = mapped[0]
                sys_name, std_name, sgd_id = get_gene_display_names(mapped_gene, info_df)

                profile_df = extract_significant_morphological_profile(
                    mapped_gene,
                    q_threshold,
                    z_values_df,
                    q_values_df,
                    parameter_desc_df,
                )

                if profile_df.empty:
                    st.info(f"No CalMorph parameters passed q ≤ {format_q_value_sci2(q_threshold)} for {profile_gene_raw}.")
                else:
                    profile_display_df = format_q_values_for_display(profile_df)
                    st.success(f"Found {len(profile_df):,} significant CalMorph parameters at q ≤ {format_q_value_sci2(q_threshold)}.")
                    base_name = std_name if std_name else mapped_gene
                    try:
                        export_box = st.popover("Export", use_container_width=False)
                    except Exception:
                        export_box = st.expander("Export", expanded=False)
                    with export_box:
                        st.download_button(
                            label="TSV",
                            data=profile_display_df.to_csv(sep="	", index=False).encode("utf-8"),
                            file_name=f"{base_name}_Morphological_Profile.tsv",
                            mime="text/tab-separated-values",
                            key="download_morph_profile_table_tsv",
                            use_container_width=True,
                        )
                        st.download_button(
                            label="CSV",
                            data=profile_display_df.to_csv(index=False).encode("utf-8"),
                            file_name=f"{base_name}_Morphological_Profile.csv",
                            mime="text/csv",
                            key="download_morph_profile_table_csv",
                            use_container_width=True,
                        )
                    st.dataframe(profile_display_df, hide_index=True, use_container_width=True)

    st.divider()

                                                     
    st.markdown("### Schematic summary of morphological structures associated with significantly altered CalMorph traits")
    st.markdown(
        """
        <div style="display:flex; gap:18px; flex-wrap:wrap; align-items:center; margin: 0.25rem 0 1rem 0;">
          <div style="display:flex; align-items:center; gap:6px;"><span style="display:inline-block;width:18px;height:18px;background:#1B5E20;border:1px solid #777;"></span>Cell wall, positive Z</div>
          <div style="display:flex; align-items:center; gap:6px;"><span style="display:inline-block;width:18px;height:18px;background:#C8E6C9;border:1px solid #777;"></span>Cell wall, negative Z</div>
          <div style="display:flex; align-items:center; gap:6px;"><span style="display:inline-block;width:18px;height:18px;background:#B71C1C;border:1px solid #777;"></span>Actin, positive Z</div>
          <div style="display:flex; align-items:center; gap:6px;"><span style="display:inline-block;width:18px;height:18px;background:#FFCDD2;border:1px solid #777;"></span>Actin, negative Z</div>
          <div style="display:flex; align-items:center; gap:6px;"><span style="display:inline-block;width:18px;height:18px;background:#0D47A1;border:1px solid #777;"></span>Nucleus, positive Z</div>
          <div style="display:flex; align-items:center; gap:6px;"><span style="display:inline-block;width:18px;height:18px;background:#BBDEFB;border:1px solid #777;"></span>Nucleus, negative Z</div>
          <div style="display:flex; align-items:center; gap:6px;"><span style="display:inline-block;width:18px;height:18px;background:#BDBDBD;border:1px solid #777;"></span>No significant mapped trait</div>
        </div>
        """,
        unsafe_allow_html=True,
    )
    color_summary = summarize_profile_colors(profile_df)
    phase_cols = st.columns(3)
    for phase, col in zip(["G1", "S/G2", "M"], phase_cols):
        with col:
            st.components.v1.html(
                render_stage_svg_html(
                    phase,
                    color_summary.get(phase, {}),
                ),
                height=350,
                scrolling=False,
            )

                               
              
                               
with tab5:
    st.header("Visualizing functional relationships between S. cerevisiae genes with morphological defects")
    st.write(
        "MorphoNet is an interactive platform for exploring quantitative morphological phenotypes in Saccharomyces cerevisiae. Built on high-resolution microscopy data "
        "from essential and nonessential gene perturbations analyzed using the CalMorph pipeline, the platform enables systematic investigation of genotype–phenotype relationships. "
        "MorphoNet integrates the curated SCMD2 morphological database to construct similarity networks that reveal functional modules and shared cellular phenotypes. "
        "By organizing high-dimensional morphological traits into intuitive network visualizations, MorphoNet provides a scalable framework for linking gene function to cellular architecture. "
        "Even though morphological data are available for 1,112 essential and 4,704 nonessential mutants in this database, the constructed networks include only a subset of these genes (513 essential and 2,911 nonessential), rather than the full datasets. This restriction arises from the requirement for detectable morphological defects and the availability of functional annotations. "
        "Networks only include mutants with significant morphological phenotypes and genes annotated with GO terms. Further details on network construction can be found in the referenced studies."
    )

    st.header("References")
    st.subheader("CalMorph")
    st.write("Ohya Y, Sese J, Yukawa M, Sano F, Nakatani Y, Saito TL, et al. (2005) High-dimensional and large-scale phenotyping of yeast mutants. "
             "PNAS. DOI: https://doi.org/10.1073/pnas.0509436102")
    st.subheader("Essential Genes Network")
    st.write("Ohnuki S and Ohya Y (2018) High-dimensional single-cell phenotyping reveals extensive haploinsufficiency. "
             "PLOS Biology. DOI: https://doi.org/10.1371/journal.pbio.2005130")
    st.subheader("NonEssential Genes Network")
    st.write("Ghanegolmohammadi F, Ohnuki S, and Ohya Y (2022) Assignment of unimodal probability distribution models for quantitative morphological phenotyping. "
             "BMC Biology. DOI: https://doi.org/10.1186/s12915-022-01283-6")
    
    st.header("Contributors")
    st.write("Farzan Ghanegolmohammadi")
    st.write("Shinsuke Ohnukia")
    st.write("Yoshikazu Ohya")

    st.header("Contact")
    st.write("Yoshikazu Ohya, PhD")
    st.write("Department of Integrated Biosciences, Graduate School of Frontier Sciences, University of Tokyo")
    st.write("Department of Science and Technology Innovation, Nagaoka University of Technology")
    st.write("E-mail: ohya@vos.nagaokaut.ac.jp")
    
    st.subheader("Acknowledgments")
    st.write("Thanks to the following libraries and their developers:")
    
    col1, col2, col3, col4, col5, col6, col7 = st.columns(7)
    with col1:
        st.link_button("Streamlit", "https://streamlit.io/")
    with col2:
        st.link_button("NumPy", "https://numpy.org/")
    with col3:
        st.link_button("Pandas", "https://pandas.pydata.org/")
    with col4:
        st.link_button("NetworkX", "https://networkx.org/")
    with col5:
        st.link_button("Matplotlib", "https://matplotlib.org/")
    with col6:
        st.link_button("Sigma.js", "https://www.sigmajs.org/")
    with col7:
        st.link_button("Graphology", "https://graphology.github.io/")
    