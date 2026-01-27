"""
Results page for NonBDNAFinder application.
Displays analysis results with publication-quality visualizations.
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

logger = logging.getLogger(__name__)


def render():
    """Render the Results tab content."""
    # Apply Results tab theme based on configuration
    load_css(TAB_THEMES.get('Results', 'genomic_purple'))
    
    # Modern header for results page
    st.markdown("""
    <div style='text-align: center; margin-bottom: 2rem;'>
        <h2 style='margin: 0; font-size: 2rem; background: linear-gradient(135deg, #a855f7, #8b5cf6);
                   -webkit-background-clip: text; -webkit-text-fill-color: transparent;
                   font-weight: 700;'>
            Analysis Results & Visualizations
        </h2>
        <p style='margin: 0.5rem 0 0 0; color: #64748b; font-size: 1rem;'>
            Comprehensive Non-B DNA motif detection results with publication-quality visualizations
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Deterministic Results Page: Only render, never compute
    # If results are missing, show info and stop
    if not st.session_state.results:
        st.info(UI_TEXT['status_no_results'])
        st.info("Run analysis first in the 'Upload & Analyze' tab")
        st.stop()  # Explicit stop to prevent any further execution
    
    # Performance metrics display if available
    if st.session_state.get('performance_metrics'):
        metrics = st.session_state.performance_metrics
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, #faf5ff 0%, #f3e8ff 100%); 
                    padding: 1.5rem; border-radius: 16px; margin-bottom: 2rem;
                    border: 1px solid #e9d5ff; box-shadow: 0 4px 16px rgba(168, 85, 247, 0.1);'>
            <h3 style='margin: 0 0 1.2rem 0; color: #7c3aed; font-size: 1.3rem; font-weight: 600;'>
                Performance Metrics
            </h3>
            <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 1rem;'>
                <div style='background: white; padding: 1rem; border-radius: 12px; text-align: center;
                            border: 1px solid #e9d5ff; box-shadow: 0 2px 8px rgba(168, 85, 247, 0.08);'>
                    <div style='font-size: 1.8rem; font-weight: 700; 
                                background: linear-gradient(135deg, #a855f7, #8b5cf6);
                                -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                        {metrics['total_time']:.2f}s
                    </div>
                    <div style='color: #64748b; font-size: 0.85rem; margin-top: 0.3rem; font-weight: 500;'>
                        Processing Time
                    </div>
                </div>
                <div style='background: white; padding: 1rem; border-radius: 12px; text-align: center;
                            border: 1px solid #e9d5ff; box-shadow: 0 2px 8px rgba(168, 85, 247, 0.08);'>
                    <div style='font-size: 1.8rem; font-weight: 700; 
                                background: linear-gradient(135deg, #a855f7, #8b5cf6);
                                -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                        {metrics['total_bp']:,}
                    </div>
                    <div style='color: #64748b; font-size: 0.85rem; margin-top: 0.3rem; font-weight: 500;'>
                        Base Pairs
                    </div>
                </div>
                <div style='background: white; padding: 1rem; border-radius: 12px; text-align: center;
                            border: 1px solid #e9d5ff; box-shadow: 0 2px 8px rgba(168, 85, 247, 0.08);'>
                    <div style='font-size: 1.8rem; font-weight: 700; 
                                background: linear-gradient(135deg, #a855f7, #8b5cf6);
                                -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                        {metrics['speed']:,.0f}
                    </div>
                    <div style='color: #64748b; font-size: 0.85rem; margin-top: 0.3rem; font-weight: 500;'>
                        bp/second
                    </div>
                </div>
                <div style='background: white; padding: 1rem; border-radius: 12px; text-align: center;
                            border: 1px solid #e9d5ff; box-shadow: 0 2px 8px rgba(168, 85, 247, 0.08);'>
                    <div style='font-size: 1.8rem; font-weight: 700; 
                                background: linear-gradient(135deg, #a855f7, #8b5cf6);
                                -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                        {metrics.get('detector_count', 9)}
                    </div>
                    <div style='color: #64748b; font-size: 0.85rem; margin-top: 0.3rem; font-weight: 500;'>
                        Detector Processes
                    </div>
                </div>
                <div style='background: white; padding: 1rem; border-radius: 12px; text-align: center;
                            border: 1px solid #e9d5ff; box-shadow: 0 2px 8px rgba(168, 85, 247, 0.08);'>
                    <div style='font-size: 1.8rem; font-weight: 700; 
                                background: linear-gradient(135deg, #a855f7, #8b5cf6);
                                -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                        {metrics['sequences']}
                    </div>
                    <div style='color: #64748b; font-size: 0.85rem; margin-top: 0.3rem; font-weight: 500;'>
                        Sequences
                    </div>
                </div>
                <div style='background: white; padding: 1rem; border-radius: 12px; text-align: center;
                            border: 1px solid #e9d5ff; box-shadow: 0 2px 8px rgba(168, 85, 247, 0.08);'>
                    <div style='font-size: 1.8rem; font-weight: 700; 
                                background: linear-gradient(135deg, #a855f7, #8b5cf6);
                                -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                        {metrics['total_motifs']}
                    </div>
                    <div style='color: #64748b; font-size: 0.85rem; margin-top: 0.3rem; font-weight: 500;'>
                        Total Motifs
                    </div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Enhanced summary display
    st.markdown(f"### {UI_TEXT['heading_analysis_summary']}")
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
        
        # Enhanced summary card with modern research-quality styling
        st.markdown(f"""
        <div class='progress-panel progress-panel--results'>
            <h3 class='progress-panel__title progress-panel__title--large'>
                NBDScanner Analysis Results
            </h3>
            <div class='stats-grid stats-grid--extra-wide'>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat_card_value_large'>
                        {stats.get("Coverage%", 0):.2f}%
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Sequence Coverage
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {stats.get("Density", 0):.2f}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Motif Density<br>(motifs/kb)
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {motif_count}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Total Motifs
                    </p>
                </div>
                <div class='stat-card stat-card--large'>
                    <h2 class='stat-card__value stat-card__value--large'>
                        {sequence_length:,}
                    </h2>
                    <p class='stat-card__label stat-card__label--large'>
                        Sequence Length (bp)
                    </p>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Add info about hybrid/cluster motifs being shown separately
        if hybrid_cluster_count > 0:
            st.info(f"{hybrid_cluster_count} Hybrid/Cluster motifs detected. View them in the 'Dynamic Clusters' tab below.")
        
        # Show cached visualization summary if available
        viz_cache_key = f"seq_{seq_idx}"
        cached_viz = st.session_state.get('cached_visualizations', {}).get(viz_cache_key, {})
        if cached_viz.get('summary'):
            viz_summary = cached_viz['summary']
            st.success(f"""**Pre-generated Analysis Ready:** 
            {viz_summary['unique_classes']} unique classes, 
            {viz_summary['unique_subclasses']} unique subclasses analyzed
            """)
        
        # NATURE-READY VISUALIZATION SUITE
        st.markdown(f'<h3>{UI_TEXT["heading_results_viz"]}</h3>', unsafe_allow_html=True)
        
        # Scientific Transparency Badge
        st.info(TRANSPARENCY_NOTE)
        
        # Create simplified visualization tabs: "All Motifs" and "Dynamic Clusters"
        viz_tabs = st.tabs([
            "All Motifs", 
            "Dynamic Clusters"
        ])
        
        # Check if clusters exist
        has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in filtered_motifs)
        
        # =================================================================
        # ALL MOTIFS TAB: Global Non-B DNA Landscape & Structural Constraints
        # =================================================================
        with viz_tabs[0]:
            st.markdown("#### Global Non-B DNA Landscape")
            st.caption("*Purpose: What structures exist, where are they located, and what are their physical constraints?*")
            
            # Panel A: Motif Composition (nested donut + bar plots)
            st.markdown("##### Motif Composition (Class → Subclass)")
            
            # Nested pie chart for hierarchical view
            try:
                fig_composition = plot_nested_pie_chart(
                    filtered_motifs, 
                    title=f"Motif Composition - {sequence_name}"
                )
                st.pyplot(fig_composition)
                plt.close(fig_composition)
            except Exception as e:
                st.error(f"Error generating composition plot: {e}")
            
            # Bar plots for class and subclass distributions
            col1, col2 = st.columns(2)
            with col1:
                try:
                    fig_class_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Class',
                        title=f"Motif Class Distribution - {sequence_name}"
                    )
                    st.pyplot(fig_class_bar)
                    plt.close(fig_class_bar)
                except Exception as e:
                    st.error(f"Error generating class bar plot: {e}")
            
            with col2:
                try:
                    fig_subclass_bar = plot_motif_distribution(
                        filtered_motifs,
                        by='Subclass',
                        title=f"Motif Subclass Distribution - {sequence_name}"
                    )
                    st.pyplot(fig_subclass_bar)
                    plt.close(fig_subclass_bar)
                except Exception as e:
                    st.error(f"Error generating subclass bar plot: {e}")
            
            # Panel B: Genome-Scale Localization (size-dependent: Manhattan OR Linear)
            st.markdown("##### Genome-Scale Localization")
            try:
                if sequence_length > 50000:
                    # Large sequences: Manhattan plot
                    st.caption("*Manhattan plot for large genome (>50kb)*")
                    fig_position = plot_manhattan_motif_density(
                        filtered_motifs, sequence_length,
                        title=f"Motif Density Hotspots - {sequence_name}"
                    )
                    st.pyplot(fig_position)
                    plt.close(fig_position)
                else:
                    # Small sequences: Linear track
                    st.caption("*Linear track for short sequence (≤50kb)*")
                    fig_position = plot_linear_motif_track(
                        filtered_motifs, sequence_length,
                        title=f"Motif Track - {sequence_name}"
                    )
                    st.pyplot(fig_position)
                    plt.close(fig_position)
            except Exception as e:
                st.error(f"Error generating localization plot: {e}")
            
            # Panel C: Genome Coverage (compact bar chart)
            st.markdown("##### Genome Coverage (% per motif class)")
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
                
                # Compact visualization: density comparison
                fig_density = plot_density_comparison(
                    genomic_density, positional_density_kbp,
                    title="Motif Density Analysis"
                )
                st.pyplot(fig_density)
                plt.close(fig_density)
                    
            except Exception as e:
                st.error(f"Error calculating density metrics: {e}")
            
            # Structural Constraints Section
            st.markdown("---")
            st.markdown("##### Structural Constraints")
            st.caption("*Purpose: What are the physical constraints shaping each motif class?*")
            
            # Length Distribution by Class (KDE)
            st.markdown("**Length Distributions by Class**")
            try:
                fig_length = plot_motif_length_kde(
                    filtered_motifs,
                    by_class=True,
                    title=f"Length Distribution (KDE) - {sequence_name}"
                )
                st.pyplot(fig_length)
                plt.close(fig_length)
            except Exception as e:
                st.error(f"Error generating length KDE: {e}")
            
            # Score distribution
            st.markdown("**Score Distribution by Class**")
            try:
                fig_score = plot_score_distribution(
                    filtered_motifs, by_class=True,
                    title="Score Distribution (1-3 Scale)"
                )
                st.pyplot(fig_score)
                plt.close(fig_score)
            except Exception as e:
                st.error(f"Error generating score distribution: {e}")
        
        # =================================================================
        # DYNAMIC CLUSTERS TAB: Structural Clustering & Co-Occurrence
        # =================================================================
        with viz_tabs[1]:
            st.markdown("#### Structural Clustering & Co-Occurrence")
            st.caption("*Purpose: Which structures co-localize and form regulatory hotspots?*")
            
            # Cluster Size Distribution (conditional on cluster existence)
            if has_clusters:
                st.markdown("##### Cluster Size Distribution")
                try:
                    fig_cluster = plot_cluster_size_distribution(
                        filtered_motifs,
                        title=f"Cluster Statistics - {sequence_name}"
                    )
                    st.pyplot(fig_cluster)
                    plt.close(fig_cluster)
                except Exception as e:
                    st.error(f"Error generating cluster plot: {e}")
            else:
                st.info("No clusters detected (requires multiple motifs in close proximity)")
            
            # Motif Co-occurrence Matrix (always shown)
            st.markdown("##### Motif Co-occurrence Matrix")
            st.caption("*Shows which motif classes tend to appear together (overlapping or within 1bp)*")
            try:
                fig_cooccur = plot_motif_cooccurrence_matrix(
                    filtered_motifs,
                    title=f"Co-occurrence Matrix - {sequence_name}"
                )
                st.pyplot(fig_cooccur)
                plt.close(fig_cooccur)
            except Exception as e:
                st.error(f"Error generating co-occurrence matrix: {e}")
