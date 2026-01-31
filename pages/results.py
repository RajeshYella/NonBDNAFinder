"""
Results page for NonBDNAFinder application.
Displays analysis results with publication-quality visualizations.

Visual Compression Strategy:
- Inline section bars with emoji + thin dividers (not headings)
- Single horizontal metrics ribbon at top
- 2-panel grid layouts for related plots
- Toggle controls instead of repeated captions
- Collapsed preview for data tables
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logging
from collections import Counter

from config.text import UI_TEXT
from config.themes import TAB_THEMES
from ui.css import load_css
from utilities import (
    get_basic_stats, 
    export_results_to_dataframe, 
    optimize_dataframe_memory,
    calculate_genomic_density, 
    calculate_positional_density,
    # Visualization functions
    plot_motif_distribution, 
    plot_nested_pie_chart,
    plot_manhattan_motif_density,
    plot_linear_motif_track,
    plot_density_comparison,
    plot_cluster_size_distribution,
    plot_motif_cooccurrence_matrix,
    plot_motif_length_kde,
    plot_score_distribution
)
from nonbscanner import get_motif_info as get_motif_classification_info
from visualization_standards import (
    NATURE_MOTIF_COLORS, 
    PlotDominance, 
    FigurePanel, 
    MetricFilter,
    LabelPolicy, 
    UILayout, 
    TRANSPARENCY_NOTE, 
    SUPPLEMENTARY_NOTE,
    should_show_plot, 
    get_nature_style_params
)


def _render_section_bar(emoji: str, title: str, info_icon: bool = False) -> None:
    """Render a compact inline section bar with emoji and thin divider."""
    info_html = '<span style="margin-left: 8px; cursor: help;" title="More information available" aria-label="Information">ℹ️</span>' if info_icon else ''
    st.markdown(f"""
    <div class="section-bar">
        <span class="section-bar__emoji" aria-hidden="true">{emoji}</span>
        <span class="section-bar__title">{title}</span>
        {info_html}
    </div>
    <div class="section-divider" aria-hidden="true"></div>
    """, unsafe_allow_html=True)


def _render_metrics_ribbon(metrics: dict) -> None:
    """Render a compact horizontal metrics ribbon."""
    # Build metrics HTML
    items = []
    
    if 'total_time' in metrics:
        items.append(f'<span class="ribbon-item">⏱ {metrics["total_time"]:.1f}s</span>')
    
    if 'total_bp' in metrics:
        # Format as Mb or kb depending on size
        bp = metrics['total_bp']
        if bp >= 1_000_000:
            bp_str = f'{bp/1_000_000:.1f} Mb'
        else:
            bp_str = f'{bp/1_000:.1f} kb'
        items.append(f'<span class="ribbon-item">🧬 {bp_str}</span>')
    
    if 'speed' in metrics:
        speed = metrics['speed']
        if speed >= 1_000_000:
            speed_str = f'{speed/1_000_000:.1f} Mb/s'
        elif speed >= 1_000:
            speed_str = f'{speed/1_000:.0f} kb/s'
        else:
            speed_str = f'{speed:.0f} bp/s'
        items.append(f'<span class="ribbon-item">⚡ {speed_str}</span>')
    
    if 'total_motifs' in metrics:
        items.append(f'<span class="ribbon-item">🔬 {metrics["total_motifs"]:,} motifs</span>')
    
    if 'detector_count' in metrics or 'sequences' in metrics:
        det_count = metrics.get('detector_count', 9)
        items.append(f'<span class="ribbon-item">🧪 {det_count} detectors</span>')
    
    ribbon_html = ' '.join(items)
    st.markdown(f"""
    <div class="metrics-ribbon">
        {ribbon_html}
    </div>
    """, unsafe_allow_html=True)

logger = logging.getLogger(__name__)


def render():
    """Render the Results tab content."""
    # Apply Results tab theme based on configuration
    load_css(TAB_THEMES.get('Results', 'genomic_purple'))
    
    # Compact header for results page
    st.markdown("""
    <div style='text-align: center; margin-bottom: 1rem;'>
        <h2 style='margin: 0; font-size: 1.8rem; background: linear-gradient(135deg, #a855f7, #8b5cf6);
                   -webkit-background-clip: text; -webkit-text-fill-color: transparent;
                   font-weight: 700;'>
            📊 Analysis Results
        </h2>
    </div>
    """, unsafe_allow_html=True)
    
    # Deterministic Results Page: Only render, never compute
    # If results are missing, show info and stop
    if not st.session_state.results:
        st.info(UI_TEXT['status_no_results'])
        st.info("Run analysis first in the 'Upload & Analyze' tab")
        st.stop()  # Explicit stop to prevent any further execution
    
    # Compact horizontal metrics ribbon (replacing large performance metrics block)
    if st.session_state.get('performance_metrics'):
        metrics = st.session_state.performance_metrics
        _render_metrics_ribbon(metrics)
    
    # Compact Analysis Summary section bar
    _render_section_bar("📋", "Analysis Summary")
    
    # Show selected class/subclass filter information (compact)
    selected_classes_used = st.session_state.get('selected_classes_used', [])
    selected_subclasses_used = st.session_state.get('selected_subclasses_used', [])
    analysis_mode_used = st.session_state.get('analysis_mode_used', 'Motif Level')
    
    if selected_classes_used:
        st.markdown(f"""
        <div class="filter-badge">
            <strong>Filter:</strong> {len(selected_classes_used)} classes / {len(selected_subclasses_used)} subclasses ({analysis_mode_used})
        </div>
        """, unsafe_allow_html=True)
    
    # Collapsed data table preview
    with st.expander(f"View summary table ({len(st.session_state.summary_df)} rows)", expanded=False):
        st.dataframe(st.session_state.summary_df, use_container_width=True)
    
    # Sequence selection for detailed analysis using pills for better UX
    seq_idx = 0
    if len(st.session_state.seqs) > 1:
        # Use pills for sequence selection - a more modern and visual alternative to dropdown
        try:
            selected_seq = st.pills(
                "Choose Sequence for Details:",
                options=list(range(len(st.session_state.seqs))),
                format_func=lambda i: st.session_state.names[i],
                selection_mode="single",
                default=0,
                help="Select a sequence to view detailed analysis results"
            )
            seq_idx = selected_seq or 0
        except Exception:
            seq_idx = 0
    
    motifs = st.session_state.results[seq_idx]
    sequence_length = len(st.session_state.seqs[seq_idx])
    sequence_name = st.session_state.names[seq_idx]
    
    if not motifs:
        st.warning("No motifs detected for this sequence.")
    else:
        # Show all motifs including hybrid/cluster motifs
        # No filtering is applied - all results are displayed
        filtered_motifs = motifs
        hybrid_cluster_motifs = [m for m in motifs if m.get('Class') in ['Hybrid', 'Non-B_DNA_Clusters']]
        
        # Create enhanced motifs DataFrame
        df = pd.DataFrame(filtered_motifs) if filtered_motifs else pd.DataFrame()
        
        # Memory optimization: Optimize DataFrame for large result sets
        if len(df) > 1000:
            df = optimize_dataframe_memory(df)
            logger.debug(f"Optimized DataFrame memory for {len(df)} motifs")
        
        # Calculate and display enhanced coverage statistics (using filtered motifs)
        stats = get_basic_stats(st.session_state.seqs[seq_idx], filtered_motifs)
        
        motif_count = len(filtered_motifs)
        hybrid_cluster_count = len(hybrid_cluster_motifs)
        coverage_pct = stats.get("Motif Coverage %", 0)
        non_b_density = (motif_count / sequence_length * 1000) if sequence_length > 0 else 0
        
        # Compact summary card (horizontal metrics strip style)
        st.markdown(f"""
        <div class="stats-ribbon">
            <div class="stats-ribbon__item">
                <span class="stats-ribbon__value">{stats.get("Coverage%", 0):.2f}%</span>
                <span class="stats-ribbon__label">Coverage</span>
            </div>
            <div class="stats-ribbon__item">
                <span class="stats-ribbon__value">{stats.get("Density", 0):.2f}</span>
                <span class="stats-ribbon__label">motifs/kb</span>
            </div>
            <div class="stats-ribbon__item">
                <span class="stats-ribbon__value">{motif_count}</span>
                <span class="stats-ribbon__label">Motifs</span>
            </div>
            <div class="stats-ribbon__item">
                <span class="stats-ribbon__value">{sequence_length:,}</span>
                <span class="stats-ribbon__label">bp</span>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Compact hybrid/cluster info badge
        if hybrid_cluster_count > 0:
            st.markdown(f'<div class="info-badge">🔗 {hybrid_cluster_count} Hybrid/Cluster motifs → Dynamic Clusters tab</div>', unsafe_allow_html=True)
        
        # Show cached visualization summary if available (compact)
        viz_cache_key = f"seq_{seq_idx}"
        cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
        if cached_viz.get('summary'):
            viz_summary = cached_viz['summary']
            st.markdown(f'<div class="success-badge">✅ {viz_summary["unique_classes"]} classes, {viz_summary["unique_subclasses"]} subclasses ready</div>', unsafe_allow_html=True)
        
        # Visualization section bar (replaces heading)
        _render_section_bar("📈", "Visual Landscape", info_icon=True)
        
        # Single scientific transparency note (compact badge, not info box)
        st.markdown(f'<div class="transparency-badge">📊 Only non-redundant metrics shown. Full data in exports.</div>', unsafe_allow_html=True)
        
        # Create simplified visualization tabs (text-first for accessibility)
        viz_tabs = st.tabs(["All Motifs", "Dynamic Clusters"])
        
        # Check if clusters exist
        has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in filtered_motifs)
        
        # =================================================================
        # ALL MOTIFS TAB: Global Non-B DNA Landscape & Structural Constraints
        # =================================================================
        with viz_tabs[0]:
            # Compact section: Motif Composition (2-panel grid)
            _render_section_bar("🧬", "Motif Composition")
            
            # 2-panel grid: Donut + Class Bar side by side
            col_donut, col_bars = st.columns([1, 1])
            
            with col_donut:
                try:
                    fig_composition = plot_nested_pie_chart(
                        filtered_motifs, 
                        title="Class → Subclass"
                    )
                    st.pyplot(fig_composition)
                    plt.close(fig_composition)
                except Exception as e:
                    st.error(f"Error: {e}")
            
            with col_bars:
                # Stacked bar plots in right column
                try:
                    fig_class_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Class',
                        title="Class Distribution"
                    )
                    st.pyplot(fig_class_bar)
                    plt.close(fig_class_bar)
                except Exception as e:
                    st.error(f"Error: {e}")
                
                try:
                    fig_subclass_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Subclass',
                        title="Subclass Distribution"
                    )
                    st.pyplot(fig_subclass_bar)
                    plt.close(fig_subclass_bar)
                except Exception as e:
                    st.error(f"Error: {e}")
            
            # Genome Localization: ONE plot with toggle
            _render_section_bar("📍", "Genome Localization")
            
            # Toggle for Manhattan vs Linear (user controls representation)
            plot_type = st.radio(
                "Genome visualization view mode",
                options=["Manhattan", "Linear"],
                index=0 if sequence_length > 50000 else 1,
                horizontal=True,
                help="Manhattan: density hotspots view for large genomes. Linear: track view for small sequences."
            )
            
            try:
                if plot_type == "Manhattan":
                    fig_position = plot_manhattan_motif_density(
                        filtered_motifs, sequence_length,
                        title=f"Density Hotspots"
                    )
                else:
                    fig_position = plot_linear_motif_track(
                        filtered_motifs, sequence_length,
                        title=f"Motif Track"
                    )
                st.pyplot(fig_position)
                plt.close(fig_position)
            except Exception as e:
                st.error(f"Error: {e}")
            
            # Structural Constraints: ONE row with 3 plots (grid layout)
            _render_section_bar("📐", "Structural Constraints")
            
            # 3-column grid for peer diagnostics
            col_length, col_score, col_density = st.columns(3)
            
            with col_length:
                try:
                    fig_length = plot_motif_length_kde(
                        filtered_motifs,
                        by_class=True,
                        title="Length KDE"
                    )
                    st.pyplot(fig_length)
                    plt.close(fig_length)
                except Exception as e:
                    st.error(f"Error: {e}")
            
            with col_score:
                try:
                    fig_score = plot_score_distribution(
                        filtered_motifs, by_class=True,
                        title="Score (1-3)"
                    )
                    st.pyplot(fig_score)
                    plt.close(fig_score)
                except Exception as e:
                    st.error(f"Error: {e}")
            
            with col_density:
                try:
                    # Calculate or retrieve density metrics
                    viz_cache_key = f"seq_{seq_idx}"
                    cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
                    cached_densities = cached_viz.get('densities', {})
                    
                    if cached_densities:
                        genomic_density = cached_densities['class_genomic']
                        positional_density_kbp = cached_densities['class_positional']
                    else:
                        genomic_density = calculate_genomic_density(filtered_motifs, sequence_length, by_class=True)
                        positional_density_kbp = calculate_positional_density(filtered_motifs, sequence_length, unit='kbp', by_class=True)
                    
                    fig_density = plot_density_comparison(
                        genomic_density, positional_density_kbp,
                        title="Density"
                    )
                    st.pyplot(fig_density)
                    plt.close(fig_density)
                except Exception as e:
                    st.error(f"Error: {e}")
        
        # =================================================================
        # DYNAMIC CLUSTERS TAB: Compact layout, no intro text
        # =================================================================
        with viz_tabs[1]:
            # No intro paragraph - just plots
            
            # Cluster Size Distribution (conditional)
            if has_clusters:
                _render_section_bar("📊", "Cluster Size Distribution")
                try:
                    fig_cluster = plot_cluster_size_distribution(
                        filtered_motifs,
                        title="Cluster Statistics"
                    )
                    st.pyplot(fig_cluster)
                    plt.close(fig_cluster)
                except Exception as e:
                    st.error(f"Error: {e}")
            else:
                st.markdown('<div class="info-badge">ℹ️ No clusters detected</div>', unsafe_allow_html=True)
            
            # Co-occurrence Matrix
            _render_section_bar("🔗", "Co-occurrence Matrix")
            try:
                fig_cooccur = plot_motif_cooccurrence_matrix(
                    filtered_motifs,
                    title="Class Co-occurrence"
                )
                st.pyplot(fig_cooccur)
                plt.close(fig_cooccur)
            except Exception as e:
                st.error(f"Error: {e}")
