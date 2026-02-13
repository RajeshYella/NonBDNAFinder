"""
╔══════════════════════════════════════════════════════════════════════════════╗
║                       RESULTS PAGE MODULE                                     ║
║              Analysis Results & Visualization Renderer                        ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: results.py (UI/)
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2024.1
LICENSE: MIT

DESCRIPTION:
    Renders analysis results with comprehensive visualizations including:
    - Motif track visualizations (class/subclass)
    - Distribution charts (bar, pie, KDE)
    - Density analysis plots
    - Cluster and co-occurrence analysis

VISUALIZATIONS:
    | Type                | Description                    |
    |---------------------|--------------------------------|
    | Linear Track        | Genomic position mapping       |
    | Distribution        | Class/Subclass bar charts      |
    | Nested Pie          | Hierarchical composition       |
    | Density Comparison  | Genomic vs positional          |
    | Length KDE          | Length distribution curves     |
    | Co-occurrence       | Motif interaction matrix       |

PERFORMANCE: Optimized DataFrame memory usage for large result sets
"""

# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logging
from collections import Counter
from Utilities.config.text import UI_TEXT
from Utilities.config.themes import TAB_THEMES
from Utilities.config.analysis import MAX_OVERLAP_DISPLAY
from UI.css import load_css
from UI.headers import render_section_heading
from Utilities.utilities import (
    get_basic_stats, export_results_to_dataframe, optimize_dataframe_memory,
    calculate_genomic_density, calculate_positional_density, plot_motif_distribution,
    plot_nested_pie_chart, plot_linear_motif_track, plot_linear_subclass_track,
    plot_density_comparison, plot_cluster_size_distribution, plot_motif_cooccurrence_matrix,
    plot_motif_length_kde, plot_score_distribution
)
from Utilities.visualization import NATURE_MOTIF_COLORS

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
logger = logging.getLogger(__name__)
CLUSTER_CLASSES = ['Hybrid', 'Non-B_DNA_Clusters']

def _render_section_divider(label): st.markdown(f"<div style='display:flex;align-items:center;gap:8px;padding:2px 0;margin-top:6px;'><span style='font-size:0.8rem;color:#64748b;font-weight:600;'>{label}</span><div style='flex:1;height:1px;background:linear-gradient(90deg,#a855f7 0%,transparent 100%);'></div></div>", unsafe_allow_html=True)

def _render_analysis_summary_box(cov, den, cnt, slen): st.markdown(f"<div style='display:flex;flex-wrap:wrap;gap:4px;padding:5px 10px;background:linear-gradient(135deg,#faf5ff 0%,#f3e8ff 100%);border-radius:6px;border:1px solid #e9d5ff;margin-bottom:8px;justify-content:space-around;align-items:center;'><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{cov:.2f}%</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>Coverage</span></div><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{den:.2f}</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>Motifs/kb</span></div><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{cnt:,}</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>Motifs</span></div><div style='display:flex;flex-direction:column;align-items:center;padding:1px 8px;'><span style='font-size:0.95rem;font-weight:800;background:linear-gradient(135deg,#a855f7,#8b5cf6);-webkit-background-clip:text;-webkit-text-fill-color:transparent;'>{slen:,}</span><span style='font-size:0.6rem;color:#64748b;text-transform:uppercase;'>bp</span></div></div>", unsafe_allow_html=True)

def _calculate_overlaps(motifs, by='Class'):
    overlaps = {}; sorted_m = sorted(motifs, key=lambda m: m.get('Start', 0))
    for i, m1 in enumerate(sorted_m):
        for m2 in sorted_m[i+1:]:
            if m2.get('Start', 0) >= m1.get('End', 0): break
            k1, k2 = m1.get(by, 'Unknown'), m2.get(by, 'Unknown')
            if k1 != k2: pair = tuple(sorted([k1, k2])); overlaps[pair] = overlaps.get(pair, 0) + 1
    return overlaps

def _render_overlap_matrix(overlaps, title):
    if not overlaps: return
    sorted_ovl = sorted(overlaps.items(), key=lambda x: x[1], reverse=True)[:MAX_OVERLAP_DISPLAY]
    html = f"<div style='font-size:0.8rem;font-weight:600;margin-bottom:4px;color:#334155;'>{title}</div><div style='font-size:0.75rem;color:#64748b;'>"
    for (k1, k2), c in sorted_ovl: html += f'<div style="padding:2px 0;">{k1.replace("_"," ")} ↔ {k2.replace("_"," ")}: <strong>{c}</strong></div>'
    html += "</div>"; st.markdown(html, unsafe_allow_html=True)

def render():
    load_css(TAB_THEMES.get('Results', 'genomic_purple')); render_section_heading("Analysis Results & Visualization", page="Results")
    if not st.session_state.results: st.info(UI_TEXT['status_no_results']); st.info("Run analysis first in 'Upload & Analyze' tab."); st.markdown("<div style='background:linear-gradient(135deg,#faf5ff 0%,#f3e8ff 100%);padding:1.2rem;border-radius:12px;margin-top:0.8rem;border:1px solid #e9d5ff;text-align:center;'><h3 style='color:#7c3aed;margin:0 0 0.6rem 0;'>Visualization Preview</h3><p style='color:#6b7280;margin:0;'>Upload and analyze a sequence to see motif track visualizations, distributions, density plots, and more.</p></div>", unsafe_allow_html=True); return
    with st.expander(f"Summary ({len(st.session_state.summary_df)} rows)", expanded=False): st.dataframe(st.session_state.summary_df, use_container_width=True)
    seq_idx = 0
    if len(st.session_state.seqs) > 1:
        try: sel = st.pills("Sequence:", options=list(range(len(st.session_state.seqs))), format_func=lambda i: st.session_state.names[i], selection_mode="single", default=0); seq_idx = sel or 0
        except: seq_idx = 0
    motifs = st.session_state.results[seq_idx]; slen = len(st.session_state.seqs[seq_idx])
    if not motifs: st.warning("No motifs detected."); return
    df = pd.DataFrame(motifs) if motifs else pd.DataFrame()
    if len(df) > 1000: df = optimize_dataframe_memory(df)
    stats = get_basic_stats(st.session_state.seqs[seq_idx], motifs); _render_analysis_summary_box(stats.get("Coverage%", 0), stats.get("Density", 0), len(motifs), slen)
    has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in motifs); has_hybrids = any(m.get('Class') == 'Hybrid' for m in motifs)
    viz_tabs = st.tabs(["All Motifs", "Dynamic Clusters"])
    with viz_tabs[0]:
        _render_section_divider("Track")
        try: fig = plot_linear_motif_track(motifs, slen, title="Class Track"); st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(f"Track error: {e}")
        _render_section_divider("Subclass"); sm = [m for m in motifs if m.get('Class') not in CLUSTER_CLASSES]
        if sm:
            try: fig = plot_linear_subclass_track(sm, slen, title="Subclass Track"); st.pyplot(fig); plt.close(fig)
            except Exception as e: st.error(f"Subclass track error: {e}")
        else: st.info("No non-cluster motifs.")
        _render_section_divider("Distribution"); c1, c2 = st.columns(2)
        with c1:
            try: fig = plot_motif_distribution(motifs, by='Class', title="Class"); st.pyplot(fig); plt.close(fig)
            except Exception as e: st.error(f"Class dist error: {e}")
        with c2:
            try: fig = plot_motif_distribution(motifs, by='Subclass', title="Subclass"); st.pyplot(fig); plt.close(fig)
            except Exception as e: st.error(f"Subclass dist error: {e}")
        _render_section_divider("Density")
        try:
            vk = f"seq_{seq_idx}"; cv = st.session_state.get('cached_visualizations', {}).get(vk, {}); cd = cv.get('densities', {})
            if cd: gd, pd_kbp = cd['class_genomic'], cd['class_positional']
            else: gd = calculate_genomic_density(motifs, slen, by_class=True); pd_kbp = calculate_positional_density(motifs, slen, unit='kbp', by_class=True)
            fig = plot_density_comparison(gd, pd_kbp, title="Density Analysis"); st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(f"Density error: {e}")
        _render_section_divider("Length")
        try: fig = plot_motif_length_kde(motifs, by_class=True, title="Length KDE"); st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(f"KDE error: {e}")
        _render_section_divider("Score")
        try: fig = plot_score_distribution(motifs, by_class=True, title="Score (1-3)"); st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(f"Score error: {e}")
        _render_section_divider("Composition")
        try: fig = plot_nested_pie_chart(motifs, title="Class → Subclass"); st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(f"Pie error: {e}")
    with viz_tabs[1]:
        if has_clusters or has_hybrids:
            _render_section_divider("Clusters"); chm = [m for m in motifs if m.get('Class') in CLUSTER_CLASSES]
            try: fig = plot_linear_motif_track(chm, slen, title="Hybrid & Cluster Track"); st.pyplot(fig); plt.close(fig)
            except Exception as e: st.error(f"Cluster track error: {e}")
            if has_clusters:
                _render_section_divider("Stats")
                try: fig = plot_cluster_size_distribution(motifs, title="Cluster Statistics"); st.pyplot(fig); plt.close(fig)
                except Exception as e: st.error(f"Cluster size error: {e}")
        else: st.markdown("<div style='padding:8px;background:#f0f9ff;border-radius:6px;color:#0369a1;font-size:0.8rem;text-align:center;'>No clusters or hybrids detected</div>", unsafe_allow_html=True)
        _render_section_divider("Co-occurrence")
        try: fig = plot_motif_cooccurrence_matrix(motifs, title="Co-occurrence"); st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(f"Co-occurrence error: {e}")
        _render_section_divider("Overlaps"); co, so = _calculate_overlaps(motifs, by='Class'), _calculate_overlaps(motifs, by='Subclass'); c1, c2 = st.columns(2)
        with c1:
            if co: _render_overlap_matrix(co, "Class Overlaps")
            else: st.markdown('<div style="color:#64748b;font-size:0.8rem;">No class overlaps</div>', unsafe_allow_html=True)
        with c2:
            if so: _render_overlap_matrix(so, "Subclass Overlaps")
            else: st.markdown('<div style="color:#64748b;font-size:0.8rem;">No subclass overlaps</div>', unsafe_allow_html=True)
