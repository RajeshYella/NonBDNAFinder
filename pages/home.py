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
    
    # ========== PROFESSIONAL HEADER WITH ENHANCED GRADIENT STYLING ==========
    st.markdown(f"""
    <div style='background: linear-gradient(135deg, #6366f1 0%, #8b5cf6 35%, #a855f7 65%, #ec4899 100%); 
                padding: 3rem 2.5rem; border-radius: 24px; margin-bottom: 2rem; 
                box-shadow: 0 20px 60px rgba(99, 102, 241, 0.3), 0 8px 24px rgba(168, 85, 247, 0.2);
                text-align: center;
                position: relative; overflow: hidden;'>
        <div style='position: absolute; top: 0; left: 0; right: 0; bottom: 0; 
                    background: url("data:image/svg+xml,%3Csvg width=\'60\' height=\'60\' xmlns=\'http://www.w3.org/2000/svg\'%3E%3Cpath d=\'M0 30 Q15 10 30 30 T60 30\' stroke=\'rgba(255,255,255,0.15)\' fill=\'none\' stroke-width=\'2\'/%3E%3C/svg%3E") repeat;
                    opacity: 0.5;'></div>
        <div style='position: relative; z-index: 1;'>
            <h1 style='color: white; font-size: {FONT_CONFIG['h1_size']}; font-weight: {FONT_CONFIG['extrabold_weight']}; 
                       margin: 0 0 1rem 0; font-family: {FONT_CONFIG['primary_font']}; 
                       letter-spacing: -0.03em; text-shadow: 0 4px 20px rgba(0,0,0,0.2);
                       background: none; -webkit-text-fill-color: white;'>
                {UI_TEXT['home_title']}
            </h1>
            <p style='color: rgba(255,255,255,0.95); font-size: 1.2rem; margin: 0; font-weight: 500;
                      max-width: 700px; margin: 0 auto; line-height: 1.6;'>
                Comprehensive Analysis of Non-B DNA Structures in Genomic Sequences
            </p>
            <div style='display: flex; justify-content: center; gap: 12px; margin-top: 1.5rem;'>
                <span style='background: rgba(255,255,255,0.2); padding: 8px 16px; border-radius: 20px;
                            font-size: 0.85rem; color: white; font-weight: 600;'>🧬 11 Motif Classes</span>
                <span style='background: rgba(255,255,255,0.2); padding: 8px 16px; border-radius: 20px;
                            font-size: 0.85rem; color: white; font-weight: 600;'>⚡ High Performance</span>
                <span style='background: rgba(255,255,255,0.2); padding: 8px 16px; border-radius: 20px;
                            font-size: 0.85rem; color: white; font-weight: 600;'>📊 Publication Ready</span>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
        
    # ========== MAIN CONTENT GRID ==========
    left, right = st.columns([1, 1], gap="large")
    
    with left:
        st.markdown(f"""
        <div style='background: linear-gradient(180deg, #ffffff 0%, #fafbfc 100%); 
                    padding: 2rem; border-radius: 20px; 
                    box-shadow: 0 4px 24px rgba(0,0,0,0.06), 0 1px 2px rgba(0,0,0,0.04); 
                    border: 1px solid #e2e8f0; height: 100%;
                    transition: all 0.3s ease;'>
            <div style='display: flex; align-items: center; gap: 10px; margin-bottom: 1.2rem;'>
                <div style='width: 4px; height: 28px; background: linear-gradient(180deg, #6366f1, #a855f7); border-radius: 2px;'></div>
                <h2 style='color: #1e293b; font-size: 1.5rem; margin: 0; font-weight: 700; letter-spacing: -0.02em;'>
                    Scientific Foundation
                </h2>
            </div>
            <p style='color: #475569; font-size: 1rem; line-height: 1.8; margin-bottom: 1.2rem;'>
                <b style='color: #6366f1;'>Non-canonical DNA structures</b> are critical regulatory elements 
                implicated in genome stability, transcriptional regulation, replication, and disease mechanisms. 
                These structures deviate from the canonical B-form DNA helix and play essential roles in:
            </p>
            <ul style='color: #475569; font-size: 0.95rem; line-height: 1.8; padding-left: 1.2rem; margin: 0;'>
                <li style='margin-bottom: 0.5rem;'><b style='color: #ef4444;'>Genome Instability:</b> Hotspots for mutations and chromosomal rearrangements</li>
                <li style='margin-bottom: 0.5rem;'><b style='color: #f59e0b;'>Gene Regulation:</b> Promoter and enhancer activity modulation</li>
                <li style='margin-bottom: 0.5rem;'><b style='color: #10b981;'>DNA Replication:</b> Origins of replication and fork progression</li>
                <li style='margin-bottom: 0.5rem;'><b style='color: #8b5cf6;'>Disease Association:</b> Cancer, neurological disorders, and aging</li>
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
            # Placeholder if image not found - using elegant gradient
            st.markdown(f"""
            <div style='background: linear-gradient(135deg, #6366f1 0%, #8b5cf6 50%, #a855f7 100%); 
                        border-radius: 20px; padding: 50px; text-align: center; color: white; margin-top: 1.5rem;
                        box-shadow: 0 10px 40px rgba(99, 102, 241, 0.3);'>
                <h2 style='margin: 0; color: white; font-size: 2rem; font-weight: 700;'>{UI_TEXT['home_image_fallback_title']}</h2>
                <h3 style='margin: 12px 0 0 0; color: rgba(255,255,255,0.9); font-weight: 500;'>{UI_TEXT['home_image_fallback_subtitle']}</h3>
                <p style='margin: 8px 0 0 0; color: rgba(255,255,255,0.8);'>{UI_TEXT['home_image_fallback_caption']}</p>
            </div>
            """, unsafe_allow_html=True)
    
    with right:
        # NOTE: Motif Classes visualization uses specific color gradients per motif type
        # These match the VISUALIZATION_PALETTE defined in centralized tokens but must
        # be literal values here due to Streamlit's inline HTML constraints.
        # Colors are carefully coordinated with the centralized token system.
        st.markdown(f"""
        <div style='background: linear-gradient(180deg, #ffffff 0%, #fafbfc 100%); 
                    padding: 2rem; border-radius: 20px; 
                    box-shadow: 0 4px 24px rgba(0,0,0,0.06), 0 1px 2px rgba(0,0,0,0.04);
                    border: 1px solid #e2e8f0; margin-bottom: 1.5rem;'>
            <div style='display: flex; align-items: center; gap: 10px; margin-bottom: 1.2rem;'>
                <div style='width: 4px; height: 28px; background: linear-gradient(180deg, #10b981, #06b6d4); border-radius: 2px;'></div>
                <h2 style='color: #1e293b; font-size: 1.5rem; margin: 0; font-weight: 700; letter-spacing: -0.02em;'>
                    Detected Motif Classes
                </h2>
            </div>
            <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 10px; margin-top: 1rem;'>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%); 
                            border-radius: 12px; border-left: 4px solid #f59e0b;
                            box-shadow: 0 2px 8px rgba(245, 158, 11, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #92400e; font-size: 0.95rem;'>1. Curved DNA</div>
                    <div style='color: #a16207; font-size: 0.8rem; margin-top: 4px;'>A-tract curvature</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%); 
                            border-radius: 12px; border-left: 4px solid #3b82f6;
                            box-shadow: 0 2px 8px rgba(59, 130, 246, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #1e40af; font-size: 0.95rem;'>2. Slipped DNA</div>
                    <div style='color: #1d4ed8; font-size: 0.8rem; margin-top: 4px;'>Direct repeats, STRs</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #fce7f3 0%, #fbcfe8 100%); 
                            border-radius: 12px; border-left: 4px solid #ec4899;
                            box-shadow: 0 2px 8px rgba(236, 72, 153, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #9f1239; font-size: 0.95rem;'>3. Cruciform</div>
                    <div style='color: #be185d; font-size: 0.8rem; margin-top: 4px;'>Palindromic IRs</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%); 
                            border-radius: 12px; border-left: 4px solid #10b981;
                            box-shadow: 0 2px 8px rgba(16, 185, 129, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #065f46; font-size: 0.95rem;'>4. R-Loop</div>
                    <div style='color: #047857; font-size: 0.8rem; margin-top: 4px;'>RNA-DNA hybrids</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #fef9c3 0%, #fef08a 100%); 
                            border-radius: 12px; border-left: 4px solid #eab308;
                            box-shadow: 0 2px 8px rgba(234, 179, 8, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #713f12; font-size: 0.95rem;'>5. Triplex</div>
                    <div style='color: #854d0e; font-size: 0.8rem; margin-top: 4px;'>Mirror repeats</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #ede9fe 0%, #ddd6fe 100%); 
                            border-radius: 12px; border-left: 4px solid #8b5cf6;
                            box-shadow: 0 2px 8px rgba(139, 92, 246, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #5b21b6; font-size: 0.95rem;'>6. G-Quadruplex</div>
                    <div style='color: #6d28d9; font-size: 0.8rem; margin-top: 4px;'>7 subtypes</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #ffedd5 0%, #fed7aa 100%); 
                            border-radius: 12px; border-left: 4px solid #f97316;
                            box-shadow: 0 2px 8px rgba(249, 115, 22, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #7c2d12; font-size: 0.95rem;'>7. i-Motif</div>
                    <div style='color: #9a3412; font-size: 0.8rem; margin-top: 4px;'>C-rich structures</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #e0e7ff 0%, #c7d2fe 100%); 
                            border-radius: 12px; border-left: 4px solid #6366f1;
                            box-shadow: 0 2px 8px rgba(99, 102, 241, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #3730a3; font-size: 0.95rem;'>8. Z-DNA</div>
                    <div style='color: #4338ca; font-size: 0.8rem; margin-top: 4px;'>Left-handed helix</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #ccfbf1 0%, #99f6e4 100%); 
                            border-radius: 12px; border-left: 4px solid #14b8a6;
                            box-shadow: 0 2px 8px rgba(20, 184, 166, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #134e4a; font-size: 0.95rem;'>9. A-philic DNA</div>
                    <div style='color: #115e59; font-size: 0.8rem; margin-top: 4px;'>A/T-rich regions</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #f3e8ff 0%, #e9d5ff 100%); 
                            border-radius: 12px; border-left: 4px solid #a855f7;
                            box-shadow: 0 2px 8px rgba(168, 85, 247, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #6b21a8; font-size: 0.95rem;'>10. Hybrid</div>
                    <div style='color: #7c3aed; font-size: 0.8rem; margin-top: 4px;'>Multi-class overlap</div>
                </div>
                <div style='padding: 12px 14px; background: linear-gradient(135deg, #f1f5f9 0%, #e2e8f0 100%); 
                            border-radius: 12px; border-left: 4px solid #64748b;
                            box-shadow: 0 2px 8px rgba(100, 116, 139, 0.15); transition: all 0.2s ease;'>
                    <div style='font-weight: 700; color: #334155; font-size: 0.95rem;'>11. Clusters</div>
                    <div style='color: #475569; font-size: 0.8rem; margin-top: 4px;'>Motif hotspots</div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        
        # Call to Action - Elegant gradient design
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, #6366f1 0%, #8b5cf6 50%, #a855f7 100%); 
                    padding: 1.75rem 2rem; border-radius: 16px; text-align: center; 
                    box-shadow: 0 10px 40px rgba(99, 102, 241, 0.3);'>
            <h3 style='color: white; margin: 0 0 0.75rem 0; font-size: 1.3rem; font-weight: 700;'>
                {UI_TEXT['home_call_to_action_title']}
            </h3>
            <p style='color: rgba(255,255,255,0.9); margin: 0 0 1.25rem 0; font-size: 1rem; line-height: 1.5;'>
                {UI_TEXT['home_call_to_action_text']}
            </p>
            <div style='background: white; color: #6366f1; padding: 12px 28px; 
                        border-radius: 12px; display: inline-block; font-weight: 700; font-size: 1rem;
                        box-shadow: 0 4px 12px rgba(0,0,0,0.15);'>
                {UI_TEXT['home_call_to_action_button']}
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # ========== KEY FEATURES SECTION WITH ENHANCED DESIGN ==========
    st.markdown("<br>", unsafe_allow_html=True)
    st.markdown(f"""
    <div style='background: linear-gradient(180deg, #f8fafc 0%, #f1f5f9 100%); 
                padding: 2.5rem; border-radius: 24px; margin-top: 2rem;
                box-shadow: 0 4px 24px rgba(0,0,0,0.04); border: 1px solid #e2e8f0;'>
        <div style='text-align: center; margin-bottom: 2rem;'>
            <h2 style='color: #1e293b; font-size: 2rem; margin: 0; 
                       font-weight: 800; letter-spacing: -0.02em;
                       background: linear-gradient(135deg, #6366f1, #8b5cf6, #a855f7);
                       -webkit-background-clip: text; -webkit-text-fill-color: transparent;'>
                Key Features & Capabilities
            </h2>
        </div>
        <div style='display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 1.5rem;'>
            <div style='background: white; padding: 1.75rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.06); border: 1px solid #e2e8f0;
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 100px; height: 100px;
                            background: linear-gradient(135deg, rgba(99, 102, 241, 0.1), transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <div style='font-size: 2rem; margin-bottom: 0.75rem;'>⚡</div>
                <h3 style='color: #1e293b; font-size: 1.2rem; margin: 0 0 0.75rem 0; font-weight: 700;'>
                    High Performance
                </h3>
                <p style='color: #64748b; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    24,674 bp/s processing speed. Handles sequences up to 1GB with chunked processing. 
                    O(n) complexity for all major detectors.
                </p>
            </div>
            <div style='background: white; padding: 1.75rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.06); border: 1px solid #e2e8f0;
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 100px; height: 100px;
                            background: linear-gradient(135deg, rgba(139, 92, 246, 0.1), transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <div style='font-size: 2rem; margin-bottom: 0.75rem;'>📊</div>
                <h3 style='color: #1e293b; font-size: 1.2rem; margin: 0 0 0.75rem 0; font-weight: 700;'>
                    Publication Quality
                </h3>
                <p style='color: #64748b; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    25+ visualization types at 300 DPI resolution. Nature/NAR-compliant formats. 
                    Colorblind-friendly palettes (Wong 2011).
                </p>
            </div>
            <div style='background: white; padding: 1.75rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.06); border: 1px solid #e2e8f0;
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 100px; height: 100px;
                            background: linear-gradient(135deg, rgba(168, 85, 247, 0.1), transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <div style='font-size: 2rem; margin-bottom: 0.75rem;'>🔬</div>
                <h3 style='color: #1e293b; font-size: 1.2rem; margin: 0 0 0.75rem 0; font-weight: 700;'>
                    Scientifically Validated
                </h3>
                <p style='color: #64748b; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    Literature-based algorithms: QmRLFS, G4Hunter, Z-Seeker. 
                    Peer-reviewed methods with biological accuracy.
                </p>
            </div>
            <div style='background: white; padding: 1.75rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.06); border: 1px solid #e2e8f0;
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 100px; height: 100px;
                            background: linear-gradient(135deg, rgba(236, 72, 153, 0.1), transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <div style='font-size: 2rem; margin-bottom: 0.75rem;'>📈</div>
                <h3 style='color: #1e293b; font-size: 1.2rem; margin: 0 0 0.75rem 0; font-weight: 700;'>
                    Statistical Analysis
                </h3>
                <p style='color: #64748b; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    Density analysis, enrichment calculations, p-value computation. 
                    100-iteration sequence shuffling for validation.
                </p>
            </div>
            <div style='background: white; padding: 1.75rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.06); border: 1px solid #e2e8f0;
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 100px; height: 100px;
                            background: linear-gradient(135deg, rgba(16, 185, 129, 0.1), transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <div style='font-size: 2rem; margin-bottom: 0.75rem;'>💾</div>
                <h3 style='color: #1e293b; font-size: 1.2rem; margin: 0 0 0.75rem 0; font-weight: 700;'>
                    Multiple Export Formats
                </h3>
                <p style='color: #64748b; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
                    Excel (multi-sheet), CSV, BED, BigWig, JSON. 
                    UCSC/IGV genome browser compatible outputs.
                </p>
            </div>
            <div style='background: white; padding: 1.75rem; border-radius: 16px; 
                        box-shadow: 0 4px 16px rgba(0,0,0,0.06); border: 1px solid #e2e8f0;
                        transition: all 0.3s ease; position: relative; overflow: hidden;'>
                <div style='position: absolute; top: 0; right: 0; width: 100px; height: 100px;
                            background: linear-gradient(135deg, rgba(6, 182, 212, 0.1), transparent);
                            border-radius: 0 16px 0 100%;'></div>
                <div style='font-size: 2rem; margin-bottom: 0.75rem;'>🧬</div>
                <h3 style='color: #1e293b; font-size: 1.2rem; margin: 0 0 0.75rem 0; font-weight: 700;'>
                    Comprehensive Coverage
                </h3>
                <p style='color: #64748b; font-size: 0.95rem; line-height: 1.7; margin: 0;'>
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
    <div style='background: linear-gradient(180deg, #ffffff 0%, #fafbfc 100%); 
                padding: 2rem; border-radius: 20px; 
                box-shadow: 0 4px 24px rgba(0,0,0,0.06); border: 1px solid #e2e8f0; margin-top: 1rem;'>
        <div style='display: flex; align-items: center; gap: 10px; margin-bottom: 1.2rem;'>
            <div style='width: 4px; height: 28px; background: linear-gradient(180deg, #6366f1, #a855f7); border-radius: 2px;'></div>
            <h2 style='color: #6366f1; font-size: 1.5rem; margin: 0; font-weight: 700;'>
                How to Cite
            </h2>
        </div>
        <div style='background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%); 
                    padding: 1.25rem; border-radius: 12px; 
                    border-left: 4px solid #6366f1; 
                    font-family: "JetBrains Mono", "Fira Code", monospace; 
                    font-size: 0.9rem; line-height: 1.8; color: #475569;'>
            <b style='color: #1e293b;'>NonBFinder: Comprehensive Detection and Analysis of Non-B DNA Motifs</b><br>
            Dr. Venkata Rajesh Yella<br>
            GitHub: <a href="https://github.com/VRYella/NonBFinder" style="color: #6366f1; text-decoration: none; font-weight: 600;">https://github.com/VRYella/NonBFinder</a><br>
            Email: <a href="mailto:yvrajesh_bt@kluniversity.in" style="color: #6366f1; text-decoration: none;">yvrajesh_bt@kluniversity.in</a>
        </div>
        <p style='color: #64748b; font-size: 0.9rem; margin-top: 1rem; line-height: 1.6;'>
            If you use NonBFinder in your research, please cite this resource. 
            For methodology references, see the <b style='color: #6366f1;'>Documentation</b> tab.
        </p>
    </div>
    """, unsafe_allow_html=True)

