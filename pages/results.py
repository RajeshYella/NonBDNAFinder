"""
Results page for NonBDNAFinder application.
Displays analysis results with publication-quality visualizations.

Visual-First Scientific Dashboard:
- Compact horizontal analysis summary (icons only)
- Plot order: Track → Subclass Track → Distributions → Density → KDE → Score → Pie
- NATURE_MOTIF_COLORS used consistently across all visualizations
- No redundant text, titles, or explanatory prose
"""

import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import logging
from collections import Counter

from config.text import UI_TEXT
from config.themes import TAB_THEMES
from config.analysis import MAX_OVERLAP_DISPLAY
from ui.css import load_css
from ui.headers import render_section_heading
from utilities import (
    get_basic_stats, 
    export_results_to_dataframe, 
    optimize_dataframe_memory,
    calculate_genomic_density, 
    calculate_positional_density,
    # Visualization functions
    plot_motif_distribution, 
    plot_nested_pie_chart,
    plot_linear_motif_track,
    plot_linear_subclass_track,
    plot_density_comparison,
    plot_cluster_size_distribution,
    plot_motif_cooccurrence_matrix,
    plot_motif_length_kde,
    plot_score_distribution
)
from visualization import NATURE_MOTIF_COLORS

logger = logging.getLogger(__name__)

# Cluster class names excluded from certain visualizations
CLUSTER_CLASSES = ['Hybrid', 'Non-B_DNA_Clusters']


def _render_section_divider(label: str) -> None:
    """Render a minimal section divider with a label.
    
    Args:
        label: Required text label for the section (e.g., 'Track', 'Density')
    """
    st.markdown(f"""
    <div style="display: flex; align-items: center; gap: 8px; padding: 2px 0; margin-top: 6px;">
        <span style="font-size: 0.8rem; color: #64748b; font-weight: 600;">{label}</span>
        <div style="flex: 1; height: 1px; background: linear-gradient(90deg, #a855f7 0%, transparent 100%);"></div>
    </div>
    """, unsafe_allow_html=True)


def _render_analysis_summary_box(coverage_pct: float, density: float, motif_count: int, seq_length: int) -> None:
    """Render compact horizontal analysis summary box (no icons/emojis)."""
    st.markdown(f"""
    <div style="display: flex; flex-wrap: wrap; gap: 6px; padding: 8px 12px; 
                background: linear-gradient(135deg, #faf5ff 0%, #f3e8ff 100%);
                border-radius: 6px; border: 1px solid #e9d5ff; margin-bottom: 10px;
                justify-content: space-around; align-items: center;">
        <div style="display: flex; flex-direction: column; align-items: center; padding: 2px 10px;">
            <span style="font-size: 1.1rem; font-weight: 800; background: linear-gradient(135deg, #a855f7, #8b5cf6);
                         -webkit-background-clip: text; -webkit-text-fill-color: transparent;">{coverage_pct:.2f}%</span>
            <span style="font-size: 0.65rem; color: #64748b; text-transform: uppercase;">Coverage</span>
        </div>
        <div style="display: flex; flex-direction: column; align-items: center; padding: 2px 10px;">
            <span style="font-size: 1.1rem; font-weight: 800; background: linear-gradient(135deg, #a855f7, #8b5cf6);
                         -webkit-background-clip: text; -webkit-text-fill-color: transparent;">{density:.2f}</span>
            <span style="font-size: 0.65rem; color: #64748b; text-transform: uppercase;">Motifs/kb</span>
        </div>
        <div style="display: flex; flex-direction: column; align-items: center; padding: 2px 10px;">
            <span style="font-size: 1.1rem; font-weight: 800; background: linear-gradient(135deg, #a855f7, #8b5cf6);
                         -webkit-background-clip: text; -webkit-text-fill-color: transparent;">{motif_count:,}</span>
            <span style="font-size: 0.65rem; color: #64748b; text-transform: uppercase;">Motifs</span>
        </div>
        <div style="display: flex; flex-direction: column; align-items: center; padding: 2px 10px;">
            <span style="font-size: 1.1rem; font-weight: 800; background: linear-gradient(135deg, #a855f7, #8b5cf6);
                         -webkit-background-clip: text; -webkit-text-fill-color: transparent;">{seq_length:,}</span>
            <span style="font-size: 0.65rem; color: #64748b; text-transform: uppercase;">bp</span>
        </div>
    </div>
    """, unsafe_allow_html=True)


def render():
    """Render the Results tab content."""
    # Apply Results tab theme
    load_css(TAB_THEMES.get('Results', 'genomic_purple'))
    
    # Uniform section heading (thin blue box with white glowing text)
    render_section_heading("Analysis Results & Visualization")
    
    # No results - show info message but continue rendering page structure
    if not st.session_state.results:
        st.info(UI_TEXT['status_no_results'])
        st.info("Run analysis first in the 'Upload & Analyze' tab to see results and visualizations.")
        
        # Show placeholder content when no results
        st.markdown("""
        <div style='background: linear-gradient(135deg, #faf5ff 0%, #f3e8ff 100%); 
                    padding: 2rem; border-radius: 12px; margin-top: 1rem;
                    border: 1px solid #e9d5ff; text-align: center;'>
            <h3 style='color: #7c3aed; margin: 0 0 1rem 0;'>Visualization Preview</h3>
            <p style='color: #6b7280; margin: 0;'>
                Once you upload and analyze a sequence, this page will display:
            </p>
            <ul style='color: #6b7280; text-align: left; margin: 1rem auto; max-width: 500px;'>
                <li>Motif track visualizations</li>
                <li>Class and subclass distributions</li>
                <li>Density analysis plots</li>
                <li>Length and score distributions</li>
                <li>Co-occurrence matrices</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
        return
    
    # Collapsed data table preview (no emoji)
    with st.expander(f"Data table ({len(st.session_state.summary_df)} rows)", expanded=False):
        st.dataframe(st.session_state.summary_df, use_container_width=True)
    
    # Sequence selection for multi-sequence files
    seq_idx = 0
    if len(st.session_state.seqs) > 1:
        try:
            selected_seq = st.pills(
                "Sequence:",
                options=list(range(len(st.session_state.seqs))),
                format_func=lambda i: st.session_state.names[i],
                selection_mode="single",
                default=0
            )
            seq_idx = selected_seq or 0
        except Exception:
            seq_idx = 0
    
    motifs = st.session_state.results[seq_idx]
    sequence_length = len(st.session_state.seqs[seq_idx])
    
    if not motifs:
        st.warning("No motifs detected for this sequence.")
        return
    
    # All motifs shown (no filtering)
    filtered_motifs = motifs
    
    # Create DataFrame for large sets
    df = pd.DataFrame(filtered_motifs) if filtered_motifs else pd.DataFrame()
    if len(df) > 1000:
        df = optimize_dataframe_memory(df)
    
    # Calculate stats
    stats = get_basic_stats(st.session_state.seqs[seq_idx], filtered_motifs)
    motif_count = len(filtered_motifs)
    coverage_pct = stats.get("Coverage%", 0)
    density = stats.get("Density", 0)
    
    # Render compact analysis summary box (icons only)
    _render_analysis_summary_box(coverage_pct, density, motif_count, sequence_length)
    
    # Check for clusters/hybrids
    has_clusters = any(m.get('Class') == 'Non-B_DNA_Clusters' for m in filtered_motifs)
    has_hybrids = any(m.get('Class') == 'Hybrid' for m in filtered_motifs)
    
    # Create visualization tabs (no emojis)
    viz_tabs = st.tabs(["All Motifs", "Dynamic Clusters"])
    
    # =================================================================
    # TAB 1: ALL MOTIFS
    # Order: Track → Subclass Track → Distributions → Density → KDE → Score → Pie
    # =================================================================
    with viz_tabs[0]:
        
        # MOTIF TRACK — Class Level
        _render_section_divider("Track")
        
        try:
            fig_track = plot_linear_motif_track(
                filtered_motifs, sequence_length,
                title="Class Track"
            )
            st.pyplot(fig_track)
            plt.close(fig_track)
        except Exception as e:
            st.error(f"Track error: {e}")
        
        # SUBCLASS TRACK (same coordinates, parent class colors, no clusters)
        _render_section_divider("Subclass")
        
        try:
            # Filter out clusters from subclass track
            subclass_motifs = [m for m in filtered_motifs 
                               if m.get('Class') not in CLUSTER_CLASSES]
            
            if subclass_motifs:
                fig_subtrack = plot_linear_subclass_track(
                    subclass_motifs, sequence_length,
                    title="Subclass Track"
                )
                st.pyplot(fig_subtrack)
                plt.close(fig_subtrack)
            else:
                st.info("No non-cluster motifs found for subclass track")
        except Exception as e:
            st.error(f"Subclass track error: {e}")
        
        # DISTRIBUTION PLOTS (Side-by-Side)
        _render_section_divider("Distribution")
        
        col_class, col_subclass = st.columns(2)
        
        with col_class:
            try:
                fig_class_dist = plot_motif_distribution(
                    filtered_motifs,
                    by='Class',
                    title="Class"
                )
                st.pyplot(fig_class_dist)
                plt.close(fig_class_dist)
            except Exception as e:
                st.error(f"Class dist error: {e}")
        
        with col_subclass:
            try:
                fig_subclass_dist = plot_motif_distribution(
                    filtered_motifs,
                    by='Subclass',
                    title="Subclass"
                )
                st.pyplot(fig_subclass_dist)
                plt.close(fig_subclass_dist)
            except Exception as e:
                st.error(f"Subclass dist error: {e}")
        
        # MOTIF DENSITY ANALYSIS
        _render_section_divider("Density")
        
        try:
            # Get cached or calculate density
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
                title="Density Analysis"
            )
            st.pyplot(fig_density)
            plt.close(fig_density)
        except Exception as e:
            st.error(f"Density error: {e}")
        
        # MOTIF LENGTH KDE
        _render_section_divider("Length")
        
        try:
            fig_kde = plot_motif_length_kde(
                filtered_motifs,
                by_class=True,
                title="Length KDE"
            )
            st.pyplot(fig_kde)
            plt.close(fig_kde)
        except Exception as e:
            st.error(f"KDE error: {e}")
        
        # MOTIF SCORE DISTRIBUTION
        _render_section_divider("Score")
        
        try:
            fig_score = plot_score_distribution(
                filtered_motifs, 
                by_class=True,
                title="Score (1-3)"
            )
            st.pyplot(fig_score)
            plt.close(fig_score)
        except Exception as e:
            st.error(f"Score error: {e}")
        
        # PIE / DONUT CHARTS (End)
        _render_section_divider("Composition")
        
        try:
            fig_pie = plot_nested_pie_chart(
                filtered_motifs, 
                title="Class → Subclass"
            )
            st.pyplot(fig_pie)
            plt.close(fig_pie)
        except Exception as e:
            st.error(f"Pie chart error: {e}")
    
    # =================================================================
    # TAB 2: DYNAMIC CLUSTERS
    # Order: Cluster Track → Co-occurrence Matrix → Overlap Analysis
    # =================================================================
    with viz_tabs[1]:
        
        # HYBRID & CLUSTER TRACK
        if has_clusters or has_hybrids:
            _render_section_divider("Clusters")
            
            # Filter to just hybrid/cluster motifs for the track
            cluster_hybrid_motifs = [m for m in filtered_motifs 
                                    if m.get('Class') in CLUSTER_CLASSES]
            
            try:
                fig_cluster_track = plot_linear_motif_track(
                    cluster_hybrid_motifs, sequence_length,
                    title="Hybrid & Cluster Track"
                )
                st.pyplot(fig_cluster_track)
                plt.close(fig_cluster_track)
            except Exception as e:
                st.error(f"Cluster track error: {e}")
            
            # Cluster Size Distribution
            if has_clusters:
                _render_section_divider("Stats")
                try:
                    fig_cluster_size = plot_cluster_size_distribution(
                        filtered_motifs,
                        title="Cluster Statistics"
                    )
                    st.pyplot(fig_cluster_size)
                    plt.close(fig_cluster_size)
                except Exception as e:
                    st.error(f"Cluster size error: {e}")
        else:
            st.markdown("""
            <div style="padding: 10px; background: #f0f9ff; border-radius: 6px; 
                        color: #0369a1; font-size: 0.85rem; text-align: center;">
                No clusters or hybrids detected in this sequence
            </div>
            """, unsafe_allow_html=True)
        
        # CO-OCCURRENCE MATRIX
        _render_section_divider("Co-occurrence")
        
        try:
            fig_cooccur = plot_motif_cooccurrence_matrix(
                filtered_motifs,
                title="Co-occurrence"
            )
            st.pyplot(fig_cooccur)
            plt.close(fig_cooccur)
        except Exception as e:
            st.error(f"Co-occurrence error: {e}")
        
        # OVERLAP ANALYSIS (Class × Class, Subclass × Subclass)
        _render_section_divider("Overlaps")
        
        # Calculate overlap stats
        class_overlaps = _calculate_overlaps(filtered_motifs, by='Class')
        subclass_overlaps = _calculate_overlaps(filtered_motifs, by='Subclass')
        
        col_class_overlap, col_subclass_overlap = st.columns(2)
        
        with col_class_overlap:
            if class_overlaps:
                _render_overlap_matrix(class_overlaps, "Class Overlaps")
            else:
                st.markdown('<div style="color: #64748b; font-size: 0.8rem;">No class overlaps</div>', 
                           unsafe_allow_html=True)
        
        with col_subclass_overlap:
            if subclass_overlaps:
                _render_overlap_matrix(subclass_overlaps, "Subclass Overlaps")
            else:
                st.markdown('<div style="color: #64748b; font-size: 0.8rem;">No subclass overlaps</div>', 
                           unsafe_allow_html=True)


def _calculate_overlaps(motifs: list, by: str = 'Class') -> dict:
    """Calculate overlapping motifs by Class or Subclass.
    
    Identifies pairs of motifs that have overlapping genomic coordinates
    and belong to different classes/subclasses. This helps visualize
    co-localization patterns between different motif types.
    
    Args:
        motifs: List of motif dictionaries, each containing at minimum:
            - 'Start': Start position (1-based)
            - 'End': End position (inclusive)
            - 'Class' or 'Subclass': Category for grouping
        by: Field to group motifs by ('Class' or 'Subclass')
        
    Returns:
        Dictionary mapping (class1, class2) tuples to overlap counts.
        Pairs are sorted alphabetically to ensure consistent keys.
        Only includes overlaps between DIFFERENT classes/subclasses.
        
    Example:
        >>> overlaps = _calculate_overlaps(motifs, by='Class')
        >>> overlaps
        {('G-Quadruplex', 'Z-DNA'): 5, ('Cruciform', 'R-Loop'): 2}
    """
    overlaps = {}
    
    # Sort motifs by start position
    sorted_motifs = sorted(motifs, key=lambda m: m.get('Start', 0))
    
    for i, m1 in enumerate(sorted_motifs):
        for m2 in sorted_motifs[i+1:]:
            # Check overlap
            if m2.get('Start', 0) >= m1.get('End', 0):
                break  # No more overlaps possible
            
            key1 = m1.get(by, 'Unknown')
            key2 = m2.get(by, 'Unknown')
            
            if key1 != key2:  # Only count different classes/subclasses
                pair = tuple(sorted([key1, key2]))
                overlaps[pair] = overlaps.get(pair, 0) + 1
    
    return overlaps


def _render_overlap_matrix(overlaps: dict, title: str) -> None:
    """Render overlap data as a compact matrix/table.
    
    Displays the top overlap pairs sorted by count in descending order.
    Limited to MAX_OVERLAP_DISPLAY entries for UI compactness.
    
    Args:
        overlaps: Dictionary from _calculate_overlaps(), mapping 
            (class1, class2) tuples to overlap counts.
        title: Title to display above the overlap list.
        
    Renders:
        Streamlit HTML with title and list of overlap pairs formatted as:
        "Class1 ↔ Class2: <count>"
    """
    if not overlaps:
        return
    
    # Sort by count descending, limited to MAX_OVERLAP_DISPLAY
    sorted_overlaps = sorted(overlaps.items(), key=lambda x: x[1], reverse=True)[:MAX_OVERLAP_DISPLAY]
    
    html = f"""
    <div style="font-size: 0.8rem; font-weight: 600; margin-bottom: 4px; color: #334155;">{title}</div>
    <div style="font-size: 0.75rem; color: #64748b;">
    """
    
    for (k1, k2), count in sorted_overlaps:
        html += f'<div style="padding: 2px 0;">{k1.replace("_", " ")} ↔ {k2.replace("_", " ")}: <strong>{count}</strong></div>'
    
    html += "</div>"
    st.markdown(html, unsafe_allow_html=True)
