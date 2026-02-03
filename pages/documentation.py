"""
Documentation page for NonBDNAFinder.
Contains motif class documentation and references.
"""

import streamlit as st
import pandas as pd
from config.text import UI_TEXT
from config.typography import FONT_CONFIG
from config.themes import TAB_THEMES
from ui.css import load_css
from ui.headers import render_section_heading

# Configuration availability flag
# This is set to False by default as configuration is not available
CONFIG_AVAILABLE = False

# These would be imported from config if available
MOTIF_LENGTH_LIMITS = {}
SCORING_METHODS = {}


def render():
    """Render the Documentation page content."""
    # Apply Documentation tab theme based on configuration
    load_css(TAB_THEMES.get('Documentation', 'midnight'))
    
    # Uniform section heading (no caption)
    render_section_heading("Scientific Documentation & References")
    
    # Motif classes documentation with elegant dark theme styling
    st.markdown("""
    <div style='background: linear-gradient(135deg, #1e1b4b 0%, #312e81 100%); 
                border-radius: 16px; padding: 24px; font-size: 1rem; 
                box-shadow: 0 10px 40px rgba(30, 27, 75, 0.4);
                border: 1px solid rgba(139, 92, 246, 0.2);'>
        <div style='display: flex; align-items: center; gap: 10px; margin-bottom: 20px;'>
            <div style='width: 4px; height: 28px; background: linear-gradient(180deg, #a855f7, #ec4899); border-radius: 2px;'></div>
            <h3 style='color: #e9d5ff; font-size: 1.4rem; margin: 0; font-weight: 700;'>
                Motif Classes Detected
            </h3>
        </div>
        <ul style='color: #e2e8f0; line-height: 1.9; padding-left: 20px; margin: 0;'>
            <li style='margin-bottom: 12px;'><b style='color: #f59e0b;'>Curved DNA</b>: Identifies phased poly(A) or poly(T) tracts using regex and spacing rules, reflecting intrinsic curvature. Scoring is based on tract length/grouping.</li>
            <li style='margin-bottom: 12px;'><b style='color: #6366f1;'>Z-DNA</b>: Detects alternating purine-pyrimidine patterns, GC-rich segments. Uses windowed scoring; regex finds dinucleotide repeats.</li>
            <li style='margin-bottom: 12px;'><b style='color: #8b5cf6;'>eGZ-motif (Extruded-G Z-DNA)</b>: Searches for long (CGG)<sub>n</sub> runs via regex. Scored by repeat count.</li>
            <li style='margin-bottom: 12px;'><b style='color: #3b82f6;'>Slipped DNA</b>: Recognizes direct/tandem repeats (10–50 nt repeat unit, 0 nt spacer) by repeat-unit matching. Scoring by length and unit copies.</li>
            <li style='margin-bottom: 12px;'><b style='color: #10b981;'>R-Loop</b>: Finds G-rich regions for stable RNA-DNA hybrids; RLFS model and regex. Thermodynamic scoring for hybrid stability.</li>
            <li style='margin-bottom: 12px;'><b style='color: #ec4899;'>Cruciform</b>: Finds palindromic inverted repeats (10–100 nt arms, 0–3 nt spacer) using reverse complement matching. Scoring by arm length and A/T content.</li>
            <li style='margin-bottom: 12px;'><b style='color: #eab308;'>Triplex DNA / Mirror Repeat</b>: Detects purine/pyrimidine mirror repeats (10–100 nt arms, 0–8 nt spacer, ≥90% purine or pyrimidine). Scoring by composition/purity.</li>
            <li style='margin-bottom: 12px;'><b style='color: #ef4444;'>Sticky DNA</b>: Searches extended GAA/TTC repeats. Scoring by repeat count.</li>
            <li style='margin-bottom: 12px;'><b style='color: #a855f7;'>G-Triplex</b>: Finds three consecutive guanine runs by regex and loop length. Scoring by G-run sum and loop penalty.</li>
            <li style='margin-bottom: 12px;'><b style='color: #22c55e;'>G4 (G-Quadruplex) and Variants</b>: Detects canonical/variant G4 motifs by G-run/loop regex. G4Hunter scoring for content/structure.</li>
            <li style='margin-bottom: 12px;'><b style='color: #f97316;'>i-Motif</b>: C-rich sequences for i-motif under acid. Regex for C runs/loops; scoring by run count and content.</li>
            <li style='margin-bottom: 12px;'><b style='color: #06b6d4;'>AC-Motif</b>: Alternating A-rich/C-rich consensus regions by regex. Scoring by pattern presence.</li>
            <li style='margin-bottom: 12px;'><b style='color: #14b8a6;'>A-philic DNA</b>: Uses tetranucleotide log2 odds scoring to identify A-tract-favoring sequences with high protein-binding affinity. Classified as high-confidence or moderate A-philic based on score.</li>
            <li style='margin-bottom: 12px;'><b style='color: #64748b;'>Hybrid Motif</b>: Regions where motif classes overlap; found by interval intersection, scored on diversity/size.</li>
            <li style='margin-bottom: 12px;'><b style='color: #475569;'>Non-B DNA Clusters</b>: Hotspots with multiple motifs in a window; sliding algorithm, scored by motif count/diversity.</li>
        </ul>
        
        <div style='margin-top: 24px; padding-top: 20px; border-top: 1px solid rgba(139, 92, 246, 0.2);'>
            <div style='display: flex; align-items: center; gap: 10px; margin-bottom: 16px;'>
                <div style='width: 4px; height: 24px; background: linear-gradient(180deg, #06b6d4, #10b981); border-radius: 2px;'></div>
                <h4 style='color: #e9d5ff; font-size: 1.2rem; margin: 0; font-weight: 700;'>
                    References
                </h4>
            </div>
            <ul style='color: #cbd5e1; line-height: 1.8; padding-left: 20px; margin: 0; font-size: 0.9rem;'>
                <li>Bedrat et al., 2016 Nucleic Acids Research</li>
                <li>Ho et al., 2010 Nature Chemical Biology</li>
                <li>Kim et al., 2018 Nucleic Acids Research</li>
                <li>Zeraati et al., 2018 Nature Chemistry</li>
                <li>Bacolla et al., 2006 Nucleic Acids Research</li>
                <li>Mirkin & Frank-Kamenetskii, 1994 Annual Review of Biophysics</li>
                <li>Vinogradov, 2003 Bioinformatics (A-philic DNA tetranucleotide analysis)</li>
                <li>Bolshoy et al., 1991 PNAS (A-tract structural properties)</li>
                <li>Rohs et al., 2009 Nature (Protein-DNA interactions, A-philic binding)</li>
                <li>New et al., 2020 Journal of DNA Structure</li>
            </ul>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Add configuration information if available
    if CONFIG_AVAILABLE:
        st.markdown("""
        <div style='background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%); 
                    border-radius: 16px; padding: 24px; margin-top: 24px;
                    border: 1px solid rgba(99, 102, 241, 0.2);'>
            <div style='display: flex; align-items: center; gap: 10px; margin-bottom: 16px;'>
                <div style='width: 4px; height: 24px; background: linear-gradient(180deg, #6366f1, #8b5cf6); border-radius: 2px;'></div>
                <h3 style='color: #e2e8f0; font-size: 1.2rem; margin: 0; font-weight: 700;'>
                    Scoring Configuration Details
                </h3>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        st.markdown("### Motif Length Constraints")
        
        config_df = pd.DataFrame([
            {
                "Motif Class": motif_class,
                "Min Length (bp)": limits["S_min"],
                "Max Length (bp)": limits["S_max"],
                "Biological Basis": f"Based on {SCORING_METHODS.get(motif_class, {}).get('reference', 'literature survey')}"
            }
            for motif_class, limits in MOTIF_LENGTH_LIMITS.items()
        ])
        
        st.dataframe(config_df, use_container_width=True)
        
        st.markdown("### Scoring Methods")
        scoring_df = pd.DataFrame([
            {
                "Motif Class": motif_class,
                "Method": method_info.get("method", ""),
                "Description": method_info.get("description", ""),
                "Reference": method_info.get("reference", "")
            }
            for motif_class, method_info in SCORING_METHODS.items()
        ])
        
        st.dataframe(scoring_df, use_container_width=True)

    st.markdown(f"""
<div style='margin-top: 36px; padding: 24px; background: linear-gradient(135deg, #312e81 0%, #1e1b4b 100%);
            border-radius: 16px; border: 1px solid rgba(139, 92, 246, 0.2);'>
    <div style='font-size: 1.1rem; color: #e2e8f0; text-align: left; font-family:{FONT_CONFIG['primary_font']};'>
        <div style='display: flex; align-items: center; gap: 10px; margin-bottom: 12px;'>
            <span style='font-size: 1.5rem;'>👨‍🔬</span>
            <b style='color: #a5b4fc;'>Developed by</b>
        </div>
        <div style='color: #f1f5f9; font-weight: 600;'>{UI_TEXT['author']}</div>
        <div style='margin-top: 8px;'>
            <a href='mailto:{UI_TEXT['author_email']}' style='color: #a78bfa; text-decoration: none;'>📧 {UI_TEXT['author_email']}</a>
            <span style='color: #64748b; margin: 0 12px;'>|</span>
            <a href='https://github.com/VRYella' target='_blank' style='color: #a78bfa; text-decoration: none;'>🔗 GitHub: VRYella</a>
        </div>
    </div>
</div>
""", unsafe_allow_html=True)
