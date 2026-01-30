"""
Results page for NonBDNAFinder application.
Displays analysis results with publication-quality visualizations.
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logging

from config.text import UI_TEXT
from config.themes import TAB_THEMES
from ui.css import load_css
from utilities import (
    get_basic_stats,
    optimize_dataframe_memory,
    calculate_genomic_density,
    calculate_positional_density,
    plot_motif_distribution,
    plot_nested_pie_chart,
    plot_manhattan_motif_density,
    plot_linear_motif_track,
    plot_density_comparison,
    plot_cluster_size_distribution,
    plot_motif_cooccurrence_matrix,
    plot_motif_length_kde,
    plot_score_distribution,
)
from visualization_standards import TRANSPARENCY_NOTE

logger = logging.getLogger(__name__)


# --------------------------- helpers ---------------------------

def show(fig):
    st.pyplot(fig)
    plt.close(fig)


def stat_card(label, value):
    st.markdown(
        f"<div class='stat-card'><b>{label}</b><br>{value}</div>",
        unsafe_allow_html=True
    )


# --------------------------- main ---------------------------

def render():
    load_css(TAB_THEMES.get("Results", "genomic_purple"))

    # Header
    st.markdown(
        """
        <div style='text-align:center;margin-bottom:2rem'>
            <h2 style='font-size:2rem;font-weight:700;
            background:linear-gradient(135deg,#a855f7,#8b5cf6);
            -webkit-background-clip:text;-webkit-text-fill-color:transparent'>
            Analysis Results & Visualizations</h2>
            <p style='color:#64748b'>
            Publication-quality Non-B DNA motif analysis
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # Guard
    if not st.session_state.get("results"):
        st.info(UI_TEXT["status_no_results"])
        st.stop()

    # Performance metrics
    if st.session_state.get("performance_metrics"):
        m = st.session_state.performance_metrics
        st.markdown("<div class='stats-grid'>", unsafe_allow_html=True)
        stat_card("Time", f"{m['total_time']:.2f}s")
        stat_card("Base Pairs", f"{m['total_bp']:,}")
        stat_card("Speed", f"{m['speed']:,.0f} bp/s")
        stat_card("Detectors", m.get("detector_count", 9))
        stat_card("Sequences", m["sequences"])
        stat_card("Motifs", m["total_motifs"])
        st.markdown("</div>", unsafe_allow_html=True)

    # Summary table
    st.markdown(f"### {UI_TEXT['heading_analysis_summary']}")
    st.dataframe(st.session_state.summary_df, use_container_width=True)

    # Sequence selector
    seq_idx = 0
    if len(st.session_state.seqs) > 1:
        seq_idx = st.pills(
            "Choose sequence:",
            list(range(len(st.session_state.seqs))),
            format_func=lambda i: st.session_state.names[i],
            default=0,
        ) or 0

    seq = st.session_state.seqs[seq_idx]
    motifs = st.session_state.results[seq_idx]
    name = st.session_state.names[seq_idx]
    length = len(seq)

    if not motifs:
        st.warning("No motifs detected.")
        return

    df = pd.DataFrame(motifs)
    if len(df) > 1000:
        df = optimize_dataframe_memory(df)

    hybrids = [m for m in motifs if m.get("Class") == "Hybrid"]
    clusters = [m for m in motifs if m.get("Class") == "Non-B_DNA_Clusters"]

    stats = get_basic_stats(seq, motifs)

    # Overview
    st.markdown("<div class='stats-grid'>", unsafe_allow_html=True)
    stat_card("Coverage", f"{stats.get('Coverage%',0):.2f}%")
    stat_card("Density", f"{stats.get('Density',0):.2f} motifs/kb")
    stat_card("Motifs", len(motifs))
    stat_card("Length", f"{length:,} bp")
    st.markdown("</div>", unsafe_allow_html=True)

    st.info(TRANSPARENCY_NOTE)

    tabs = st.tabs(["All Motifs", "Dynamic Clusters"])

    # ===================== ALL MOTIFS =====================
    with tabs[0]:
        st.markdown("## Global Non-B DNA Landscape")

        st.markdown("### Motif Composition")
        show(plot_nested_pie_chart(motifs, title=name))

        c1, c2 = st.columns(2)
        with c1:
            show(plot_motif_distribution(motifs, by="Class"))
        with c2:
            show(plot_motif_distribution(motifs, by="Subclass"))

        st.markdown("### Structural Constraints")
        show(plot_motif_length_kde(motifs, by_class=True))
        show(plot_score_distribution(motifs, by_class=True))

        st.markdown("### Genome Localization")
        if length <= 50000:
            show(plot_linear_motif_track(motifs, length))
        else:
            show(plot_manhattan_motif_density(motifs, length))

        st.markdown("### Density & Coverage")
        gd = calculate_genomic_density(motifs, length, by_class=True)
        pdn = calculate_positional_density(motifs, length, by_class=True)
        show(plot_density_comparison(gd, pdn))

    # ===================== DYNAMIC CLUSTERS =====================
    with tabs[1]:
        st.markdown("## Hybrid & Cluster Motifs")

        st.markdown("### Hybrid Motifs")
        if hybrids:
            st.metric("Hybrid Regions", len(hybrids))
            show(plot_motif_cooccurrence_matrix(hybrids))
        else:
            st.info("No hybrid motifs detected.")

        st.markdown("### Non-B DNA Clusters")
        if clusters:
            st.metric("Clusters", len(clusters))
            show(plot_cluster_size_distribution(clusters))
        else:
            st.info("No clusters detected.")

        st.markdown("### Co-occurrence (Global)")
        show(plot_motif_cooccurrence_matrix(motifs))
