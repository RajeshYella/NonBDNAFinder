"""
Results page for NonBDNAFinder application.
Displays analysis results with publication-quality visualizations.
"""

import streamlit as st; import pandas as pd; import matplotlib.pyplot as plt; import logging
from collections import Counter; from config.text import UI_TEXT; from config.themes import TAB_THEMES
from ui.css import load_css; from utilities import (get_basic_stats, export_results_to_dataframe, optimize_dataframe_memory,
    calculate_genomic_density, calculate_positional_density, plot_motif_distribution, plot_nested_pie_chart,
    plot_manhattan_motif_density, plot_linear_motif_track, plot_density_comparison, plot_cluster_size_distribution,
    plot_motif_cooccurrence_matrix, plot_motif_length_kde, plot_score_distribution)
from nonbscanner import get_motif_info as get_motif_classification_info; from visualization_standards import (
    NATURE_MOTIF_COLORS, PlotDominance, FigurePanel, MetricFilter, LabelPolicy, UILayout, TRANSPARENCY_NOTE, 
    SUPPLEMENTARY_NOTE, should_show_plot, get_nature_style_params)

logger = logging.getLogger(__name__)

def render():
    load_css(TAB_THEMES.get('Results', 'genomic_purple'))
    st.markdown("""
    <div style='text-align: center; margin-bottom: 2rem;'>
        <h2 style='margin: 0; font-size: 2rem; background: linear-gradient(135deg, #a855f7, #8b5cf6);
                   -webkit-background-clip: text; -webkit-text-fill-color: transparent; font-weight: 700;'>
            Analysis Results & Visualizations
        </h2>
        <p style='margin: 0.5rem 0 0 0; color: #64748b; font-size: 1rem;'>
            Comprehensive Non-B DNA motif detection results with publication-quality visualizations
        </p>
    </div>
    """, unsafe_allow_html=True)

    if not st.session_state.results: st.info(UI_TEXT['status_no_results']); st.info("Run analysis first in the 'Upload & Analyze' tab"); st.stop()

    if st.session_state.get('performance_metrics'):
        metrics = st.session_state.performance_metrics; st.markdown(f"""
        <div style='background: linear-gradient(135deg, #faf5ff 0%, #f3e8ff 100%); padding: 1.5rem; border-radius: 16px;'>
            <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 1rem;'>
                {[f'''<div style='background: white; padding: 1rem; border-radius: 12px; text-align: center;'>
                    <div style='font-size: 1.8rem; font-weight: 700;'>
                        {value}
                    </div><div style='color: #64748b; font-size: 0.85rem; margin-top: 0.3rem;'>
                        {label}</div></div>''' for value, label in zip(
                            [f"{metrics['total_time']:.2f}s", f"{metrics['total_bp']:,}", f"{metrics['speed']:,}",
                            metrics.get('detector_count', 9), metrics['sequences'], metrics['total_motifs']],
                            ["Processing Time", "Base Pairs", "bp/second", "Detector Processes", "Sequences", "Total Motifs"])]).join('')}
            </div>
        </div>
        """, unsafe_allow_html=True)

    seq_idx = st.session_state.get('selected_seq', 0)
    if len(motifs := st.session_state.results[seq_idx] if st.session_state.results else []) > 0:
        df = pd.DataFrame(motifs); df = optimize_dataframe_memory(df) if len(df) > 1000 else df; stats = get_basic_stats(st.session_state.seqs[seq_idx], df);
        st.markdown(f"""
        <div class='progress-panel progress-panel--results'>
            <h3 class='progress-panel__title progress-panel__title--large'>
                NBDScanner Analysis Results</h3>
            <div class="stats"> {[f''' \<Stat>9]}")
