import streamlit as st
import os
from config.text import UI_TEXT
from config.typography import FONT_CONFIG
from config.themes import TAB_THEMES
from config.colors import SEMANTIC_COLORS
from ui.css import load_css, get_page_colors


def render():
    # Apply Home theme based on configuration
    load_css(TAB_THEMES.get('Home', 'scientific_blue'))
    
    # Get page-specific colors from centralized token system for inline HTML styles
    # NOTE: Inline HTML styles cannot use CSS variables, so we inject literal values
    # from the centralized color token system defined at the top of this file.
    # This ensures all colors are still managed centrally even in inline styles.
    colors = get_page_colors('Home')
    
    # ========== THIN BLUE BOX WITH WHITE GLOWING TEXT ==========
    st.markdown(f"""
    <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                padding: 1rem 1.5rem; border-radius: 8px; margin-bottom: 1.5rem; 
                border: 2px solid {colors['primary']}; text-align: center;
                box-shadow: 0 0 15px rgba(30, 64, 175, 0.4);'>
        <p style='color: {colors['white']}; font-size: 1.1rem; font-weight: 600; 
                  margin: 0; text-shadow: 0 0 10px rgba(255,255,255,0.8), 0 0 20px rgba(255,255,255,0.6);'>
            Non-B DNA Motif Detection System: Comprehensive Analysis of Non-B DNA Structures in Genomic Sequences
        </p>
    </div>
    """, unsafe_allow_html=True)
    

    
    # ========== NBD CIRCLE IMAGE (ABOVE SCIENTIFIC FOUNDATION) ==========
    try:
        # Display NBD Circle logo
        possible_paths = ["nbdcircle.JPG", "archive/nbdcircle.JPG", "./nbdcircle.JPG"]
        image_found = False
        for img_path in possible_paths:
            if os.path.exists(img_path):
                col_img1, col_img2, col_img3 = st.columns([1, 2, 1])
                with col_img2:
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
    
    # ========== MAIN CONTENT GRID ==========
    left, right = st.columns([1, 1], gap="large")
    
    with left:
        st.markdown(f"""
        <div style='background: {colors['white']}; padding: 2rem; border-radius: 16px; 
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; height: 100%;'>
            <h2 style='color: {colors['text']}; font-size: 1.6rem; margin: 0 0 1rem 0; font-weight: 600;'>
                Scientific Foundation
            </h2>
            <p style='color: {colors['neutral_700']}; font-size: 1rem; line-height: 1.8; margin-bottom: 1.2rem;'>
                <b style='color: {colors['primary']};'>Non-canonical DNA structures</b> are critical regulatory elements 
                implicated in genome stability, transcriptional regulation, replication, and disease mechanisms. 
                These structures deviate from the canonical B-form DNA helix and play essential roles in:
            </p>
            <ul style='color: {colors['neutral_700']}; font-size: 0.95rem; line-height: 1.7; padding-left: 1.5rem;'>
                <li><b>Genome Instability:</b> Hotspots for mutations and chromosomal rearrangements</li>
                <li><b>Gene Regulation:</b> Promoter and enhancer activity modulation</li>
                <li><b>DNA Replication:</b> Origins of replication and fork progression</li>
                <li><b>Disease Association:</b> Cancer, neurological disorders, and aging</li>
            </ul>
            
        </div>
        """, unsafe_allow_html=True)
    
    with right:
        # NOTE: Compact motif display with abbreviated labels organized by category
        # Each tag has a tooltip (title attribute) with full description
        # Colors are coordinated with the centralized token system
        # Using inline styles to ensure proper rendering in Streamlit
        tag_base_style = "display: inline-block; padding: 0.3rem 0.6rem; border-radius: 6px; font-size: 0.78rem; font-weight: 500; margin: 0.15rem; cursor: help;"
        category_style = f"font-weight: 600; font-size: 0.85rem; margin: 0.6rem 0 0.3rem 0; padding-bottom: 0.2rem; border-bottom: 1px solid {colors['neutral_200']};"
        
        st.markdown(f"""
        <div style='background: {colors['white']}; padding: 1.5rem; border-radius: 16px; 
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; margin-bottom: 1.5rem;'>
            <h2 style='color: {colors['text']}; font-size: 1.4rem; margin: 0 0 0.8rem 0; font-weight: 600;'>
                Detected Motif Classes
            </h2>
            <div style="{category_style} color: {SEMANTIC_COLORS['warning_dark']};">Curvature &amp; Repeats</div>
            <div style='display: flex; flex-wrap: wrap; gap: 0.2rem;'>
                <span style="{tag_base_style} background: linear-gradient(135deg, {SEMANTIC_COLORS['warning_light']} 0%, {SEMANTIC_COLORS['warning_border']} 100%); color: {SEMANTIC_COLORS['warning_dark']}; border: 1px solid {SEMANTIC_COLORS['warning']};" title="Global Curvature: Overall DNA bending across the entire sequence">Global Curv</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, {SEMANTIC_COLORS['warning_light']} 0%, {SEMANTIC_COLORS['warning_border']} 100%); color: {SEMANTIC_COLORS['warning_dark']}; border: 1px solid {SEMANTIC_COLORS['warning']};" title="Local Curvature: Region-specific DNA bending patterns">Local Curv</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, {SEMANTIC_COLORS['info_light']} 0%, {SEMANTIC_COLORS['info_border']} 100%); color: {SEMANTIC_COLORS['info_dark']}; border: 1px solid {SEMANTIC_COLORS['info']};" title="Direct Repeat: Tandem repeated sequences in the same orientation">Dir. Repeat</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, {SEMANTIC_COLORS['info_light']} 0%, {SEMANTIC_COLORS['info_border']} 100%); color: {SEMANTIC_COLORS['info_dark']}; border: 1px solid {SEMANTIC_COLORS['info']};" title="Short Tandem Repeat: Microsatellite sequences with repeated units">STR</span>
            </div>
            <div style="{category_style} color: #9f1239;">Structural Motifs</div>
            <div style='display: flex; flex-wrap: wrap; gap: 0.2rem;'>
                <span style="{tag_base_style} background: linear-gradient(135deg, #fce7f3 0%, #fbcfe8 100%); color: #9f1239; border: 1px solid #ec4899;" title="Cruciform forming Inverted Repeats: Palindromic sequences that can form cross-shaped structures">Cruciform IR</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%); color: #065f46; border: 1px solid #10b981;" title="R-loop formation sites: Regions prone to RNA-DNA hybrid formation with displaced single-stranded DNA">R-loop Sites</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #fef9c3 0%, #fef08a 100%); color: #713f12; border: 1px solid #eab308;" title="Triplex: Triple-stranded DNA structures formed by mirror repeats">Triplex</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #fef9c3 0%, #fef08a 100%); color: #713f12; border: 1px solid #eab308;" title="Sticky DNA: Triplex-mediated structures that cause DNA to self-associate">Sticky DNA</span>
            </div>
            <div style="{category_style} color: #5b21b6;">G-Quadruplex-Related</div>
            <div style='display: flex; flex-wrap: wrap; gap: 0.2rem;'>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Telomeric G4: G-quadruplex structures found in telomeric repeat sequences">Telo G4</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Stacked canonical G4s: Multiple G-quadruplex units stacked without intervening sequences">Stacked G4</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Stacked G4s with linker: Multiple G-quadruplex units connected by linker sequences">G4 + Linker</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Canonical intramolecular G4: Single-strand G-quadruplex with standard G3+ runs and loops">Intra G4</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Extended-loop canonical G4: G-quadruplex with longer loop sequences">Ext. Loop G4</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Higher-order G4 array / G4-wire: Large-scale G-quadruplex assemblies">G4 Array</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Intramolecular G-triplex: Three G-tract structures that may fold into partial G4">Intra G-triplex</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); color: #5b21b6; border: 1px solid #8b5cf6;" title="Two-tetrad weak PQS: Putative G-quadruplex with only two G-tetrads">Weak PQS</span>
            </div>
            <div style="{category_style} color: #7c2d12;">Other Non-B Motifs</div>
            <div style='display: flex; flex-wrap: wrap; gap: 0.2rem;'>
                <span style="{tag_base_style} background: linear-gradient(135deg, #fed7aa 0%, #fdba74 100%); color: #7c2d12; border: 1px solid #f97316;" title="Canonical i-motif: C-rich structures forming intercalated C-C+ base pairs at acidic pH">i-Motif</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #fed7aa 0%, #fdba74 100%); color: #7c2d12; border: 1px solid #f97316;" title="Relaxed i-motif: i-motif variants with longer loops or fewer C-tracts">Relaxed iM</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #fed7aa 0%, #fdba74 100%); color: #7c2d12; border: 1px solid #f97316;" title="AC-motif: Alternating adenine-cytosine rich sequences">AC Motif</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #e0e7ff 0%, #c7d2fe 100%); color: #3730a3; border: 1px solid #6366f1;" title="Z-DNA: Left-handed double helix formed by alternating purine-pyrimidine sequences">Z-DNA</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #e0e7ff 0%, #c7d2fe 100%); color: #3730a3; border: 1px solid #6366f1;" title="Extended Z-DNA (eGZ): Z-DNA with extended recognition patterns">eGZ</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #ccfbf1 0%, #99f6e4 100%); color: #134e4a; border: 1px solid #14b8a6;" title="A-philic DNA: A/T-rich sequences favoring A-form DNA conformation">A-DNA</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #f3e8ff 0%, #e9d5ff 100%); color: #6b21a8; border: 1px solid #a855f7;" title="Dynamic overlaps: Regions where multiple motif types overlap">Dyn. Overlaps</span>
                <span style="{tag_base_style} background: linear-gradient(135deg, #e5e7eb 0%, #d1d5db 100%); color: #1f2937; border: 1px solid #6b7280;" title="Dynamic clusters: Genomic hotspots with high density of non-B DNA motifs">Dyn. Clusters</span>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Call to Action - Using page-specific colors
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                    padding: 1.5rem; border-radius: 12px; text-align: center; 
                    box-shadow: 0 4px 12px {colors['shadow']};'>
            <h3 style='color: {colors['white']}; margin: 0 0 0.5rem 0; font-size: 1.2rem;'>
                {UI_TEXT['home_call_to_action_title']}
            </h3>
            <p style='color: rgba(255,255,255,0.95); margin: 0 0 1rem 0; font-size: 0.95rem;'>
                {UI_TEXT['home_call_to_action_text']}
            </p>
            <div style='background: {colors['white']}; color: {colors['primary']}; padding: 0.7rem 1.5rem; 
                        border-radius: 8px; display: inline-block; font-weight: 600; font-size: 1rem;'>
                {UI_TEXT['home_call_to_action_button']}
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # ========== HOW TO CITE SECTION ==========
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"""
    <div style='background: {colors['white']}; padding: 2rem; border-radius: 16px; 
                box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; margin-top: 2rem;'>
        <h2 style='color: {colors['primary']}; font-size: 1.6rem; margin: 0 0 1rem 0; font-weight: 600;'>
            How to Cite
        </h2>
        <div style='background: {colors['neutral_50']}; padding: 1.2rem; border-radius: 8px; border-left: 4px solid {colors['primary']}; 
                    font-family: "Courier New", monospace; font-size: 0.9rem; line-height: 1.7; color: {colors['neutral_700']};'>
            <b>NonBFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs</b><br>
            Dr. Venkata Rajesh Yella<br>
            GitHub: <a href="https://github.com/VRYella/NonBFinder" style="color: {colors['primary']};">https://github.com/VRYella/NonBFinder</a><br>
            Email: yvrajesh_bt@kluniversity.in
        </div>
        <p style='color: {colors['neutral_600']}; font-size: 0.9rem; margin-top: 1rem; line-height: 1.6;'>
            If you use NonBFinder in your research, please cite this resource. 
            For methodology references, see the <b>Documentation</b> tab.
        </p>
    </div>
    """, unsafe_allow_html=True)
    

