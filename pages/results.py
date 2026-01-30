"""
Results page for NonBDNAFinder application.
Displays analysis results with publication-quality visualizations.
"""

import streamlit as st; import pandas as pd; import matplotlib.pyplot as plt; import logging
from config.text import UI_TEXT; from config.themes import TAB_THEMES; from ui.css import load_css
from utilities import (
    get_basic_stats, export_results_to_dataframe, optimize_dataframe_memory,
    calculate_genomic_density, calculate_positional_density,
    plot_motif_distribution, plot_nested_pie_chart,
    plot_manhattan_motif_density, plot_linear_motif_track,
    plot_density_comparison, plot_cluster_size_distribution,
    plot_motif_cooccurrence_matrix, plot_motif_length_kde,
    plot_score_distribution
)
from visualization_standards import TRANSPARENCY_NOTE

logger = logging.getLogger(__name__)


def render():
    load_css(TAB_THEMES.get("Results", "genomic_purple"))

    st.markdown("""
    <div style='text-align:center;margin-bottom:2rem'>
      <h2 style='font-size:2rem;font-weight:700;
      background:linear-gradient(135deg,#a855f7,#8b5cf6);
      -webkit-background-clip:text;-webkit-text-fill-color:transparent'>
      Analysis Results & Visualizations</h2>
      <p style='color:#64748b'>Publication-quality Non-B DNA analysis</p>
    </div>""", unsafe_allow_html=True)

    if not st.session_state.results:
        st.info(UI_TEXT["status_no_results"]); st.stop()

    if st.session_state.get("performance_metrics"):
        m = st.session_state.performance_metrics
        st.markdown(f"""
        <div style='padding:1.5rem;border-radius:16px;
        background:linear-gradient(135deg,#faf5ff,#f3e8ff)'>
        <b>Performance</b><br>
        ⏱ {m['total_time']:.2f}s | 🧬 {m['total_bp']:,} bp |
        ⚡ {m['speed']:,.0f} bp/s | 🔬 {m['total_motifs']} motifs
        </div>""", unsafe_allow_html=True)

    st.markdown(f"### {UI_TEXT['heading_analysis_summary']}")
    st.dataframe(st.session_state.summary_df, use_container_width=True)

    seq_idx = 0
    if len(st.session_state.seqs) > 1:
        try:
            seq_idx = st.pills(
                "Choose Sequence",
                range(len(st.session_state.seqs)),
                format_func=lambda i: st.session_state.names[i],
                default=0
            ) or 0
        except Exception: pass

    motifs = st.session_state.results[seq_idx]
    if not motifs: st.warning("No motifs detected"); return

    seq = st.session_state.seqs[seq_idx]; L = len(seq)
    name = st.session_state.names[seq_idx]; df = pd.DataFrame(motifs)
    if len(df) > 1000: df = optimize_dataframe_memory(df)

    stats = get_basic_stats(seq, motifs); st.info(TRANSPARENCY_NOTE)

    st.markdown(f"""
    <div style='display:grid;grid-template-columns:repeat(4,1fr);gap:1rem'>
      <div><b>{stats.get("Coverage%",0):.2f}%</b><br>Coverage</div>
      <div><b>{stats.get("Density",0):.2f}</b><br>Motifs/kb</div>
      <div><b>{len(motifs)}</b><br>Total Motifs</div>
      <div><b>{L:,}</b><br>bp</div>
    </div>""", unsafe_allow_html=True)

    tabs = st.tabs(["All Motifs", "Dynamic Clusters"])

    # ========================= ALL MOTIFS =========================
    with tabs[0]:
        st.subheader("Genome-scale Localization")
        try:
            fig = (
                plot_manhattan_motif_density(motifs, L, f"Density – {name}")
                if L > 50000 else
                plot_linear_motif_track(motifs, L, f"Track – {name}")
            ); st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(e)

        st.subheader("Motif Composition")
        for fn, by in [
            (plot_nested_pie_chart, None),
            (plot_motif_distribution, "Class"),
            (plot_motif_distribution, "Subclass")
        ]:
            try:
                fig = fn(motifs, by=by, title=f"{by or 'Hierarchy'} – {name}")
                st.pyplot(fig); plt.close(fig)
            except Exception: pass

        st.subheader("Coverage & Density")
        try:
            gd = calculate_genomic_density(motifs, L, by_class=True)
            pdn = calculate_positional_density(motifs, L, unit="kbp", by_class=True)
            fig = plot_density_comparison(gd, pdn, "Density Comparison")
            st.pyplot(fig); plt.close(fig)
        except Exception as e: st.error(e)

        st.subheader("Structural Constraints")
        for fn, title in [
            (plot_motif_length_kde, "Length KDE"),
            (plot_score_distribution, "Score Distribution")
        ]:
            try:
                fig = fn(motifs, by_class=True, title=title)
                st.pyplot(fig); plt.close(fig)
            except Exception: pass

    # ======================= DYNAMIC CLUSTERS ======================
    with tabs[1]:
        st.subheader("Structural Clustering")
        if any(m.get("Class") == "Non-B_DNA_Clusters" for m in motifs):
            try:
                fig = plot_cluster_size_distribution(motifs, f"Clusters – {name}")
                st.pyplot(fig); plt.close(fig)
            except Exception: pass
        else: st.info("No clusters detected")

        st.subheader("Motif Co-occurrence")
        try:
            fig = plot_motif_cooccurrence_matrix(motifs, f"Co-occurrence – {name}")
            st.pyplot(fig); plt.close(fig)
        except Exception: pass
