import streamlit as st
import os
from Utilities.config.text import UI_TEXT
from Utilities.config.typography import FONT_CONFIG
from Utilities.config.themes import TAB_THEMES
from Utilities.config.colors import (
    SEMANTIC_COLORS, 
    UNIFIED_MOTIF_COLORS, 
    MOTIF_CARD_COLORS,
    MOTIF_CLASS_INFO,
    get_motif_card_style
)
from UI.css import load_css, get_page_colors
from UI.headers import render_section_heading


def _build_motif_class_card(class_info: dict) -> str:
    """
    Build HTML for a single motif class card using centralized colors.
    
    Args:
        class_info: Dictionary with 'key', 'name', 'subtitle', 'num' keys
        
    Returns:
        HTML string for the card
    """
    style = get_motif_card_style(class_info['key'])
    return f"""<div style='padding: 0.5rem; background: {style['background']}; 
                border-radius: 8px; border-left: 4px solid {style['border']};'>
    <div style='font-weight: 600; color: {style['text']}; font-size: 0.8rem;'>{class_info['num']}. {class_info['name']}</div>
    <div style='color: {style['text']}; font-size: 0.7rem; margin-top: 0.15rem;'>{class_info['subtitle']}</div>
</div>"""


def render():
    # Apply Home theme based on configuration
    load_css(TAB_THEMES.get('Home', 'scientific_blue'))
    
    # Get page-specific colors from centralized token system for inline HTML styles
    # NOTE: Inline HTML styles cannot use CSS variables, so we inject literal values
    # from the centralized color token system defined at the top of this file.
    # This ensures all colors are still managed centrally even in inline styles.
    colors = get_page_colors('Home')
    
    # ========== SECTION HEADING (using uniform component) ==========
    render_section_heading("Non-B DNA Finder: Systematic Detection of Non-B DNA motifs", page="Home")
    

    
    # ========== TWO-COLUMN SECTION: NBD CIRCLE (COL 1) & SCIENTIFIC FOUNDATION (COL 2) ==========
    col_nbd, col_science = st.columns([1, 1], gap="large")
    
    with col_nbd:
        # ========== NBD CIRCLE IMAGE ==========
        try:
            # Display NBD Circle logo
            possible_paths = ["nbdcircle.JPG", "archive/nbdcircle.JPG", "./nbdcircle.JPG"]
            image_found = False
            for img_path in possible_paths:
                if os.path.exists(img_path):
                    st.image(img_path, caption=UI_TEXT['home_image_caption'], use_container_width=True)
                    image_found = True
                    break
            if not image_found:
                raise FileNotFoundError("Image not found")  # Intentional: triggers fallback
        except Exception:
            # Placeholder if image not found - using page colors
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                        border-radius: 15px; padding: 40px; text-align: center; color: {colors['white']}; margin-bottom: 1rem;'>
                <h2 style='margin: 0; color: {colors['white']}; font-size: 2rem;'>{UI_TEXT['home_image_fallback_title']}</h2>
                <h3 style='margin: 10px 0 0 0; color: {colors['white']};'>{UI_TEXT['home_image_fallback_subtitle']}</h3>
                <p style='margin: 5px 0 0 0; color: rgba(255,255,255,0.9);'>{UI_TEXT['home_image_fallback_caption']}</p>
            </div>
            """, unsafe_allow_html=True)
    
    with col_science:
        # ========== SCIENTIFIC FOUNDATION ==========
        st.markdown(f"""
        <div style='background: {colors['white']}; padding: 1.2rem; border-radius: 16px; 
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; height: 100%;'>
            <h2 style='color: {colors['text']}; font-size: 1.4rem; margin: 0 0 0.6rem 0; font-weight: 600;'>
                Scientific Foundation
            </h2>
            <p style='color: {colors['neutral_700']}; font-size: 0.9rem; line-height: 1.6; margin-bottom: 0.8rem;'>
                <b style='color: {colors['primary']};'>Non-canonical DNA structures</b> are critical regulatory elements 
                implicated in genome stability, transcriptional regulation, replication, and disease mechanisms. 
                These structures deviate from the canonical B-form DNA helix and play essential roles in:
            </p>
            <ul style='color: {colors['neutral_700']}; font-size: 0.85rem; line-height: 1.5; padding-left: 1.5rem;'>
                <li><b>Genome Instability:</b> Hotspots for mutations and chromosomal rearrangements</li>
                <li><b>Gene Regulation:</b> Promoter and enhancer activity modulation</li>
                <li><b>DNA Replication:</b> Origins of replication and fork progression</li>
                <li><b>Disease Association:</b> Cancer, neurological disorders, and aging</li>
            </ul>
            
        </div>
        """, unsafe_allow_html=True)
    
    # ========== TWO-COLUMN SECTION: MOTIF CLASSES & CALL TO ACTION ==========
    left, right = st.columns([1, 1], gap="large")
    
    with left:
        # Motif Classes visualization using UNIFIED colors from centralized config
        # Build cards dynamically from MOTIF_CLASS_INFO for consistency
        cards_html = ''.join(_build_motif_class_card(info) for info in MOTIF_CLASS_INFO)
        
        st.markdown(f"""
        <div style='background: {colors['white']}; padding: 1.2rem; border-radius: 16px; 
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; height: 100%;'>
            <h2 style='color: {colors['text']}; font-size: 1.4rem; margin: 0 0 0.6rem 0; font-weight: 600;'>
                Detected Motif Classes
            </h2>
            <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; margin-top: 0.6rem;'>
                {cards_html}
            </div>
        </div>
        """, unsafe_allow_html=True)
        
    with right:
        # Call to Action - Using page-specific colors
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                    padding: 1rem; border-radius: 12px; text-align: center; 
                    box-shadow: 0 4px 12px {colors['shadow']}; height: 100%;
                    display: flex; flex-direction: column; justify-content: center;'>
            <h3 style='color: {colors['white']}; margin: 0 0 0.3rem 0; font-size: 1.1rem;'>
                {UI_TEXT['home_call_to_action_title']}
            </h3>
            <p style='color: rgba(255,255,255,0.95); margin: 0 0 0.6rem 0; font-size: 0.85rem;'>
                {UI_TEXT['home_call_to_action_text']}
            </p>
            <div style='background: {colors['white']}; color: {colors['primary']}; padding: 0.5rem 1.2rem; 
                        border-radius: 8px; display: inline-block; font-weight: 600; font-size: 0.9rem; margin: 0 auto;'>
                {UI_TEXT['home_call_to_action_button']}
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # ========== HOW TO CITE SECTION ==========
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"""
    <div style='background: {colors['white']}; padding: 1.2rem; border-radius: 16px; 
                box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; margin-top: 1.2rem;'>
        <h2 style='color: {colors['primary']}; font-size: 1.4rem; margin: 0 0 0.6rem 0; font-weight: 600;'>
            How to Cite
        </h2>
        <div style='background: {colors['neutral_50']}; padding: 0.8rem; border-radius: 8px; border-left: 4px solid {colors['primary']}; 
                    font-family: "Courier New", monospace; font-size: 0.8rem; line-height: 1.5; color: {colors['neutral_700']};'>
            <b>NonBFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs</b><br>
            Dr. Venkata Rajesh Yella<br>
            GitHub: <a href="https://github.com/VRYella/NonBFinder" style="color: {colors['primary']};">https://github.com/VRYella/NonBFinder</a><br>
            Email: yvrajesh_bt@kluniversity.in
        </div>
        <p style='color: {colors['neutral_600']}; font-size: 0.8rem; margin-top: 0.6rem; line-height: 1.4;'>
            If you use NonBFinder in your research, please cite this resource. 
            For methodology references, see the <b>Documentation</b> tab.
        </p>
    </div>
    """, unsafe_allow_html=True)
    

