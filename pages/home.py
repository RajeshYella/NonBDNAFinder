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
    
    # ========== PROFESSIONAL HEADER WITH ENHANCED STYLING ==========
    st.markdown(f"""
    <div style='background: linear-gradient(135deg, {colors['primary']} 0%, {colors['secondary']} 100%); 
                padding: 3rem 2rem; border-radius: 20px; margin-bottom: 2rem; 
                box-shadow: 0 10px 40px rgba(0,0,0,0.2); text-align: center;
                position: relative; overflow: hidden;'>
        <div style='position: absolute; top: 0; left: 0; right: 0; bottom: 0; 
                    background: url("data:image/svg+xml,%3Csvg width=\'40\' height=\'40\' xmlns=\'http://www.w3.org/2000/svg\'%3E%3Cpath d=\'M0 20 Q10 10 20 20 T40 20\' stroke=\'rgba(255,255,255,0.1)\' fill=\'none\'/%3E%3C/svg%3E") repeat;
                    opacity: 0.3;'></div>
        <div style='position: relative; z-index: 1;'>
            <h1 style='color: {colors['white']}; font-size: {FONT_CONFIG['h1_size']}; font-weight: {FONT_CONFIG['extrabold_weight']}; 
                       margin: 0 0 0.8rem 0; font-family: {FONT_CONFIG['primary_font']}; 
                       letter-spacing: -0.02em; text-shadow: 0 2px 10px rgba(0,0,0,0.2);'>
                {UI_TEXT['home_title']}
            </h1>
            <p style='color: rgba(255,255,255,0.95); font-size: 1.1rem; margin: 0; font-weight: 500;'>
                Comprehensive Analysis of Non-B DNA Structures in Genomic Sequences
            </p>
        </div>
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
                        border-radius: 15px; padding: 40px; text-align: center; color: {colors['white']}; margin-top: 1rem;'>
                <h2 style='margin: 0; color: {colors['white']}; font-size: 2rem;'>{UI_TEXT['home_image_fallback_title']}</h2>
                <h3 style='margin: 10px 0 0 0; color: {colors['white']};'>{UI_TEXT['home_image_fallback_subtitle']}</h3>
                <p style='margin: 5px 0 0 0; color: rgba(255,255,255,0.9);'>{UI_TEXT['home_image_fallback_caption']}</p>
            </div>
            """, unsafe_allow_html=True)
    
    with right:
        # NOTE: Motif Classes visualization uses specific color gradients per motif type
        # These match the VISUALIZATION_PALETTE defined in centralized tokens but must
        # be literal values here due to Streamlit's inline HTML constraints.
        # Colors are carefully coordinated with the centralized token system.
        st.markdown(f"""
        <div style='background: {colors['white']}; padding: 2rem; border-radius: 16px; 
                    box-shadow: 0 4px 20px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']}; margin-bottom: 1.5rem;'>
            <h2 style='color: {colors['text']}; font-size: 1.6rem; margin: 0 0 1rem 0; font-weight: 600;'>
                Detected Motif Classes
            </h2>
            <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 0.8rem; margin-top: 1rem;'>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, {SEMANTIC_COLORS['warning_light']} 0%, {SEMANTIC_COLORS['warning_border']} 100%); 
                            border-radius: 8px; border-left: 4px solid {SEMANTIC_COLORS['warning']};'>
                    <div style='font-weight: 600; color: {SEMANTIC_COLORS['warning_dark']}; font-size: 0.9rem;'>1. Curved DNA</div>
                    <div style='color: {SEMANTIC_COLORS['warning_dark']}; font-size: 0.75rem; margin-top: 0.2rem;'>A-tract curvature</div>
                </div>
                <!-- NOTE: Motif class visualization cards use semantic gradients for visual distinction.
                     Cards 1-2 use SEMANTIC_COLORS from centralized tokens.
                     Cards 3-11 use specific Tailwind-inspired gradients - this is intentional:
                     - Each motif type requires a unique, visually distinct gradient
                     - These form a carefully designed visual vocabulary for motif recognition
                     - Adding 9 × 4 gradient variations to centralized tokens would decrease maintainability
                     All colors follow semantic principles. This is an explicit design decision. -->
                <div style='padding: 0.8rem; background: linear-gradient(135deg, {SEMANTIC_COLORS['info_light']} 0%, {SEMANTIC_COLORS['info_border']} 100%); 
                            border-radius: 8px; border-left: 4px solid {SEMANTIC_COLORS['info']};'>
                    <div style='font-weight: 600; color: {SEMANTIC_COLORS['info_dark']}; font-size: 0.9rem;'>2. Slipped DNA</div>
                    <div style='color: {SEMANTIC_COLORS['info_dark']}; font-size: 0.75rem; margin-top: 0.2rem;'>Direct repeats, STRs</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #fce7f3 0%, #fbcfe8 100%); 
                            border-radius: 8px; border-left: 4px solid #ec4899;'>
                    <div style='font-weight: 600; color: #9f1239; font-size: 0.9rem;'>3. Cruciform</div>
                    <div style='color: #831843; font-size: 0.75rem; margin-top: 0.2rem;'>Palindromic IRs</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%); 
                            border-radius: 8px; border-left: 4px solid #10b981;'>
                    <div style='font-weight: 600; color: #065f46; font-size: 0.9rem;'>4. R-Loop</div>
                    <div style='color: #064e3b; font-size: 0.75rem; margin-top: 0.2rem;'>RNA-DNA hybrids</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #fef9c3 0%, #fef08a 100%); 
                            border-radius: 8px; border-left: 4px solid #eab308;'>
                    <div style='font-weight: 600; color: #713f12; font-size: 0.9rem;'>5. Triplex</div>
                    <div style='color: #713f12; font-size: 0.75rem; margin-top: 0.2rem;'>Mirror repeats</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #ddd6fe 0%, #c4b5fd 100%); 
                            border-radius: 8px; border-left: 4px solid #8b5cf6;'>
                    <div style='font-weight: 600; color: #5b21b6; font-size: 0.9rem;'>6. G-Quadruplex</div>
                    <div style='color: #4c1d95; font-size: 0.75rem; margin-top: 0.2rem;'>7 subtypes</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #fed7aa 0%, #fdba74 100%); 
                            border-radius: 8px; border-left: 4px solid #f97316;'>
                    <div style='font-weight: 600; color: #7c2d12; font-size: 0.9rem;'>7. i-Motif</div>
                    <div style='color: #7c2d12; font-size: 0.75rem; margin-top: 0.2rem;'>C-rich structures</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #e0e7ff 0%, #c7d2fe 100%); 
                            border-radius: 8px; border-left: 4px solid #6366f1;'>
                    <div style='font-weight: 600; color: #3730a3; font-size: 0.9rem;'>8. Z-DNA</div>
                    <div style='color: #312e81; font-size: 0.75rem; margin-top: 0.2rem;'>Left-handed helix</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #ccfbf1 0%, #99f6e4 100%); 
                            border-radius: 8px; border-left: 4px solid #14b8a6;'>
                    <div style='font-weight: 600; color: #134e4a; font-size: 0.9rem;'>9. A-philic DNA</div>
                    <div style='color: #134e4a; font-size: 0.75rem; margin-top: 0.2rem;'>A/T-rich regions</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #f3e8ff 0%, #e9d5ff 100%); 
                            border-radius: 8px; border-left: 4px solid #a855f7;'>
                    <div style='font-weight: 600; color: #6b21a8; font-size: 0.9rem;'>10. Hybrid</div>
                    <div style='color: #581c87; font-size: 0.75rem; margin-top: 0.2rem;'>Multi-class overlap</div>
                </div>
                <div style='padding: 0.8rem; background: linear-gradient(135deg, #e5e7eb 0%, #d1d5db 100%); 
                            border-radius: 8px; border-left: 4px solid #6b7280;'>
                    <div style='font-weight: 600; color: #1f2937; font-size: 0.9rem;'>11. Clusters</div>
                    <div style='color: #374151; font-size: 0.75rem; margin-top: 0.2rem;'>Motif hotspots</div>
                </div>
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
    
    # ========== KEY FEATURES SECTION WITH ENHANCED DESIGN ==========
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"""
    <div style='background: {colors['neutral_50']}; padding: 2.5rem; border-radius: 20px; margin-top: 2rem;
                box-shadow: 0 4px 20px rgba(0,0,0,0.08);'>
        <h2 style='color: {colors['text']}; font-size: 2rem; margin: 0 0 2rem 0; text-align: center; 
                   font-weight: 700; background: linear-gradient(135deg, {colors['primary']}, {colors['secondary']});
                   -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
            Key Features & Capabilities
        </h2>
        <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 1.5rem;'>
            <div style='background: {colors['white']}; padding: 1.8rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']};
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 80px; height: 80px;
                            background: linear-gradient(135deg, {colors['primary']}20, transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <h3 style='color: {colors['text']}; font-size: 1.2rem; margin: 0 0 0.8rem 0; font-weight: 600;'>
                    High Performance
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    24,674 bp/s processing speed. Handles sequences up to 1GB with chunked processing. 
                    O(n) complexity for all major detectors.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.8rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']};
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 80px; height: 80px;
                            background: linear-gradient(135deg, {colors['secondary']}20, transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <h3 style='color: {colors['text']}; font-size: 1.2rem; margin: 0 0 0.8rem 0; font-weight: 600;'>
                    Publication Quality
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    25+ visualization types at 300 DPI resolution. Nature/NAR-compliant formats. 
                    Colorblind-friendly palettes (Wong 2011).
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.8rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']};
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 80px; height: 80px;
                            background: linear-gradient(135deg, {colors['primary']}20, transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <h3 style='color: {colors['text']}; font-size: 1.2rem; margin: 0 0 0.8rem 0; font-weight: 600;'>
                    Scientifically Validated
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    Literature-based algorithms: QmRLFS, G4Hunter, Z-Seeker. 
                    Peer-reviewed methods with biological accuracy.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.8rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']};
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 80px; height: 80px;
                            background: linear-gradient(135deg, {colors['secondary']}20, transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <h3 style='color: {colors['text']}; font-size: 1.2rem; margin: 0 0 0.8rem 0; font-weight: 600;'>
                    Statistical Analysis
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    Density analysis, enrichment calculations, p-value computation. 
                    100-iteration sequence shuffling for validation.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.8rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']};
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 80px; height: 80px;
                            background: linear-gradient(135deg, {colors['primary']}20, transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <h3 style='color: {colors['text']}; font-size: 1.2rem; margin: 0 0 0.8rem 0; font-weight: 600;'>
                    Multiple Export Formats
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    Excel (multi-sheet), CSV, BED, BigWig, JSON. 
                    UCSC/IGV genome browser compatible outputs.
                </p>
            </div>
            <div style='background: {colors['white']}; padding: 1.8rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.08); border: 1px solid {colors['neutral_200']};
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 80px; height: 80px;
                            background: linear-gradient(135deg, {colors['secondary']}20, transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <h3 style='color: {colors['text']}; font-size: 1.2rem; margin: 0 0 0.8rem 0; font-weight: 600;'>
                    Comprehensive Coverage
                </h3>
                <p style='color: {colors['neutral_600']}; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    11 major classes, 22+ subclasses. Hybrid and cluster detection. 
                    Complete Non-B DNA structural characterization.
                </p>
            </div>
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
    

