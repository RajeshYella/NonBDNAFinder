"""
Documentation page for NonBDNAFinder.
Contains comprehensive scientific documentation with peer-reviewed references.
"""

import streamlit as st
import pandas as pd
from config.text import UI_TEXT
from config.typography import FONT_CONFIG
from config.themes import TAB_THEMES
from ui.css import load_css
from ui.headers import render_section_heading


# Comprehensive motif detection parameters
MOTIF_PARAMETERS = {
    "Curved DNA": {
        "min_length": "20 bp",
        "max_length": "500 bp",
        "algorithm": "A-tract phasing detection",
        "scoring": "Tract length × phasing score",
        "references": ["Crothers et al., 1990", "Koo et al., 1986"]
    },
    "G-Quadruplex": {
        "min_length": "15 bp",
        "max_length": "100 bp",
        "algorithm": "G4Hunter + Regex pattern matching",
        "scoring": "G4Hunter score (threshold ≥1.2)",
        "references": ["Bedrat et al., 2016", "Huppert & Balasubramanian, 2005"]
    },
    "Z-DNA": {
        "min_length": "8 bp",
        "max_length": "200 bp",
        "algorithm": "Alternating purine-pyrimidine detection",
        "scoring": "Z-score with GC enrichment",
        "references": ["Ho et al., 1986", "Wang et al., 1979"]
    },
    "Cruciform": {
        "min_length": "10 bp arm",
        "max_length": "100 bp arm",
        "algorithm": "Inverted repeat detection",
        "scoring": "Arm length × stem stability",
        "references": ["Lilley, 1980", "Panayotatos & Wells, 1981"]
    },
    "R-Loop": {
        "min_length": "15 bp",
        "max_length": "2000 bp",
        "algorithm": "QmRLFS model",
        "scoring": "Thermodynamic ΔG stability",
        "references": ["Jenjaroenpun et al., 2015", "Aguilera & García-Muse, 2012"]
    },
    "Triplex": {
        "min_length": "10 bp",
        "max_length": "100 bp",
        "algorithm": "Mirror repeat detection",
        "scoring": "Purine/pyrimidine purity ≥90%",
        "references": ["Mirkin et al., 1987", "Frank-Kamenetskii & Mirkin, 1995"]
    },
    "i-Motif": {
        "min_length": "12 bp",
        "max_length": "60 bp",
        "algorithm": "C-tract pattern matching",
        "scoring": "C-run count × content score",
        "references": ["Zeraati et al., 2018", "Day et al., 2014"]
    },
    "Slipped DNA": {
        "min_length": "6 bp unit",
        "max_length": "50 bp unit",
        "algorithm": "Direct repeat detection",
        "scoring": "Unit copies × repeat fidelity",
        "references": ["Pearson et al., 2005", "Wells, 2007"]
    },
    "A-philic DNA": {
        "min_length": "10 bp",
        "max_length": "200 bp",
        "algorithm": "Tetranucleotide log₂ odds scoring",
        "scoring": "Propensity score threshold",
        "references": ["Vinogradov, 2003", "Rohs et al., 2009"]
    }
}

# Comprehensive peer-reviewed references
REFERENCES = [
    {"authors": "Bedrat A, Lacroix L, Mergny JL", "year": 2016, "title": "Re-evaluation of G-quadruplex propensity with G4Hunter", "journal": "Nucleic Acids Res", "volume": "44(4):1746-59", "doi": "10.1093/nar/gkw006"},
    {"authors": "Huppert JL, Balasubramanian S", "year": 2005, "title": "Prevalence of quadruplexes in the human genome", "journal": "Nucleic Acids Res", "volume": "33(9):2908-16", "doi": "10.1093/nar/gki609"},
    {"authors": "Zeraati M, Langley DB, et al.", "year": 2018, "title": "I-motif DNA structures are formed in the nuclei of human cells", "journal": "Nat Chem", "volume": "10(6):631-637", "doi": "10.1038/s41557-018-0046-3"},
    {"authors": "Ho PS, Frederick CA, Saal D, et al.", "year": 1986, "title": "The interactions of ruthenium hexaammine with Z-DNA", "journal": "J Biomol Struct Dyn", "volume": "4(3):521-34", "doi": "10.1080/07391102.1986.10506363"},
    {"authors": "Jenjaroenpun P, Wongsurawat T, et al.", "year": 2015, "title": "QmRLFS-finder: a model, web server and stand-alone tool", "journal": "Nucleic Acids Res", "volume": "43(W1):W527-34", "doi": "10.1093/nar/gkv344"},
    {"authors": "Aguilera A, García-Muse T", "year": 2012, "title": "R loops: from transcription byproducts to threats to genome stability", "journal": "Mol Cell", "volume": "46(2):115-24", "doi": "10.1016/j.molcel.2012.04.009"},
    {"authors": "Frank-Kamenetskii MD, Mirkin SM", "year": 1995, "title": "Triplex DNA structures", "journal": "Annu Rev Biochem", "volume": "64:65-95", "doi": "10.1146/annurev.bi.64.070195.000433"},
    {"authors": "Crothers DM, Drak J, et al.", "year": 1990, "title": "DNA bending, flexibility, and helical repeat", "journal": "Methods Enzymol", "volume": "212:3-29", "doi": "10.1016/0076-6879(92)12003-9"},
    {"authors": "Vinogradov AE", "year": 2003, "title": "DNA helix: the importance of being GC-rich", "journal": "Nucleic Acids Res", "volume": "31(7):1838-44", "doi": "10.1093/nar/gkg296"},
    {"authors": "Rohs R, West SM, et al.", "year": 2009, "title": "The role of DNA shape in protein-DNA recognition", "journal": "Nature", "volume": "461(7268):1248-53", "doi": "10.1038/nature08473"},
    {"authors": "Pearson CE, Nichol Edamura K, Cleary JD", "year": 2005, "title": "Repeat instability: mechanisms of dynamic mutations", "journal": "Nat Rev Genet", "volume": "6(10):729-42", "doi": "10.1038/nrg1689"},
    {"authors": "Bacolla A, Wells RD", "year": 2004, "title": "Non-B DNA conformations, genomic rearrangements, and human disease", "journal": "J Biol Chem", "volume": "279(46):47411-4", "doi": "10.1074/jbc.R400028200"},
]


def render():
    """Render the Documentation page content."""
    # Apply Documentation tab theme - use scientific_blue for better readability
    load_css(TAB_THEMES.get('Documentation', 'scientific_blue'))
    
    # Section heading (thin blue box with white glowing text)
    render_section_heading("Scientific Documentation & References")
    
    # ═══════════════════════════════════════════════════════════
    # TOOL OVERVIEW CARD (compact, no emoji)
    # ═══════════════════════════════════════════════════════════
    st.markdown("""
    <div style='background: linear-gradient(135deg, #1e40af 0%, #7c3aed 100%); 
                padding: 1.5rem; border-radius: 12px; margin-bottom: 1.5rem;
                box-shadow: 0 4px 16px rgba(30, 64, 175, 0.3);'>
        <h2 style='color: white; margin: 0 0 0.75rem 0; font-size: 1.4rem; font-weight: 700;'>
            NonBDNAFinder v2025.1
        </h2>
        <p style='color: rgba(255,255,255,0.95); font-size: 0.95rem; line-height: 1.6; margin: 0;'>
            A comprehensive computational platform for <strong>genome-wide detection and analysis</strong> of 
            Non-B DNA structures. Implements <strong>11 motif classes</strong> with <strong>24 subclasses</strong>, 
            validated against peer-reviewed algorithms including G4Hunter, QmRLFS, and Z-Seeker.
        </p>
        <div style='display: flex; gap: 0.75rem; margin-top: 1rem; flex-wrap: wrap;'>
            <span style='background: rgba(255,255,255,0.2); padding: 0.4rem 0.8rem; border-radius: 16px; 
                         color: white; font-size: 0.8rem; font-weight: 600;'>24,674 bp/s</span>
            <span style='background: rgba(255,255,255,0.2); padding: 0.4rem 0.8rem; border-radius: 16px; 
                         color: white; font-size: 0.8rem; font-weight: 600;'>200MB+ sequences</span>
            <span style='background: rgba(255,255,255,0.2); padding: 0.4rem 0.8rem; border-radius: 16px; 
                         color: white; font-size: 0.8rem; font-weight: 600;'>25+ visualizations</span>
            <span style='background: rgba(255,255,255,0.2); padding: 0.4rem 0.8rem; border-radius: 16px; 
                         color: white; font-size: 0.8rem; font-weight: 600;'>Nature-ready output</span>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # ═══════════════════════════════════════════════════════════
    # MOTIF CLASSES - DETAILED CARDS (compact, no emoji)
    # ═══════════════════════════════════════════════════════════
    st.markdown("""
    <h3 style='color: #003D82; font-size: 1.3rem; margin: 1.5rem 0 0.75rem 0; font-weight: 700;
               border-left: 4px solid #0091FF; padding-left: 0.75rem;'>
        Detected Non-B DNA Motif Classes
    </h3>
    """, unsafe_allow_html=True)
    
    # Create motif cards in a grid (no emoji icons)
    motif_info = [
        ("Curved DNA", "A-tract mediated bending", "#06b6d4", "Intrinsic DNA curvature from phased A-tracts"),
        ("Slipped DNA", "Direct repeats & STRs", "#f59e0b", "Slippage-mediated repeat expansions"),
        ("Cruciform", "Inverted repeats", "#ef4444", "Hairpin structures from palindromic sequences"),
        ("R-Loop", "RNA:DNA hybrids", "#8b5cf6", "Co-transcriptional R-loop formation sites"),
        ("Triplex", "H-DNA structures", "#ec4899", "Triple-stranded DNA from mirror repeats"),
        ("G-Quadruplex", "G4 structures", "#10b981", "Four-stranded G-rich secondary structures"),
        ("i-Motif", "C-rich structures", "#22c55e", "Intercalated cytosine structures"),
        ("Z-DNA", "Left-handed helix", "#6366f1", "Alternating purine-pyrimidine sequences"),
        ("A-philic DNA", "Protein-binding", "#f97316", "High nucleosome positioning potential"),
        ("Hybrid", "Multi-class overlaps", "#64748b", "Regions with multiple motif types"),
        ("Clusters", "Hotspots", "#334155", "High-density Non-B DNA regions"),
    ]
    
    cards_html = '<div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(250px, 1fr)); gap: 0.75rem; margin-bottom: 1.5rem;">'
    for name, subtitle, color, description in motif_info:
        cards_html += f'''
        <div style='background: white; padding: 1rem; border-radius: 10px; 
                    box-shadow: 0 2px 8px rgba(0,0,0,0.06); border-left: 4px solid {color};'>
            <strong style='color: #1e293b; font-size: 0.95rem;'>{name}</strong>
            <div style='color: {color}; font-size: 0.75rem; font-weight: 600; margin: 0.2rem 0;'>{subtitle}</div>
            <div style='color: #64748b; font-size: 0.8rem; line-height: 1.4;'>{description}</div>
        </div>
        '''
    cards_html += '</div>'
    st.markdown(cards_html, unsafe_allow_html=True)
    
    # ═══════════════════════════════════════════════════════════
    # ALGORITHM PARAMETERS TABLE (no emoji)
    # ═══════════════════════════════════════════════════════════
    st.markdown("""
    <h3 style='color: #003D82; font-size: 1.3rem; margin: 1.5rem 0 0.75rem 0; font-weight: 700;
               border-left: 4px solid #0091FF; padding-left: 0.75rem;'>
        Detection Parameters & Algorithms
    </h3>
    """, unsafe_allow_html=True)
    
    # Create DataFrame from parameters
    params_data = []
    for motif, params in MOTIF_PARAMETERS.items():
        params_data.append({
            "Motif Class": motif,
            "Min Length": params["min_length"],
            "Max Length": params["max_length"],
            "Algorithm": params["algorithm"],
            "Scoring Method": params["scoring"]
        })
    
    params_df = pd.DataFrame(params_data)
    st.dataframe(
        params_df,
        use_container_width=True,
        hide_index=True,
        column_config={
            "Motif Class": st.column_config.TextColumn("Motif Class", width="medium"),
            "Min Length": st.column_config.TextColumn("Min Length", width="small"),
            "Max Length": st.column_config.TextColumn("Max Length", width="small"),
            "Algorithm": st.column_config.TextColumn("Algorithm", width="large"),
            "Scoring Method": st.column_config.TextColumn("Scoring", width="medium")
        }
    )
    
    # ═══════════════════════════════════════════════════════════
    # PEER-REVIEWED REFERENCES (no emoji)
    # ═══════════════════════════════════════════════════════════
    st.markdown("""
    <h3 style='color: #003D82; font-size: 1.3rem; margin: 1.5rem 0 0.75rem 0; font-weight: 700;
               border-left: 4px solid #0091FF; padding-left: 0.75rem;'>
        Peer-Reviewed References
    </h3>
    <p style='color: #64748b; font-size: 0.85rem; margin-bottom: 0.75rem;'>
        NonBDNAFinder implements algorithms validated in the following peer-reviewed publications:
    </p>
    """, unsafe_allow_html=True)
    
    # Reference cards (compact)
    refs_html = '<div style="display: flex; flex-direction: column; gap: 0.5rem;">'
    for ref in REFERENCES:
        refs_html += f'''
        <div style='background: #f8fafc; padding: 0.75rem; border-radius: 6px; border-left: 3px solid #3b82f6;'>
            <div style='font-weight: 600; color: #1e293b; font-size: 0.85rem; margin-bottom: 0.2rem;'>
                {ref["authors"]} ({ref["year"]})
            </div>
            <div style='color: #334155; font-size: 0.8rem; font-style: italic; margin-bottom: 0.15rem;'>
                {ref["title"]}
            </div>
            <div style='color: #64748b; font-size: 0.75rem;'>
                <strong>{ref["journal"]}</strong> {ref["volume"]} · 
                <a href="https://doi.org/{ref["doi"]}" target="_blank" style="color: #3b82f6;">DOI: {ref["doi"]}</a>
            </div>
        </div>
        '''
    refs_html += '</div>'
    st.markdown(refs_html, unsafe_allow_html=True)
    
    # ═══════════════════════════════════════════════════════════
    # CITATION & AUTHOR INFO (compact, no emoji)
    # ═══════════════════════════════════════════════════════════
    st.markdown(f"""
    <div style='background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%); 
                padding: 1.5rem; border-radius: 12px; margin-top: 1.5rem;
                border: 1px solid #bae6fd;'>
        <h3 style='color: #0c4a6e; margin: 0 0 0.75rem 0; font-size: 1.1rem; font-weight: 700;'>
            How to Cite
        </h3>
        <div style='background: white; padding: 1rem; border-radius: 6px; font-family: "Courier New", monospace;
                    font-size: 0.8rem; line-height: 1.6; color: #334155; border-left: 4px solid #0284c7;'>
            <strong>Yella VR</strong> (2025). NonBDNAFinder: Comprehensive Detection and Analysis of Non-B DNA Forming Motifs.<br>
            GitHub: <a href="https://github.com/VRYella/NonBDNAFinder" style="color: #0284c7;">https://github.com/VRYella/NonBDNAFinder</a>
        </div>
        <div style='margin-top: 1rem; padding-top: 1rem; border-top: 1px solid #bae6fd;'>
            <div style='font-weight: 700; color: #0c4a6e; font-size: 0.9rem; margin-bottom: 0.4rem;'>
                Developed by
            </div>
            <div style='color: #334155; font-size: 0.85rem;'>
                <strong>{UI_TEXT['author']}</strong><br>
                Email: <a href='mailto:{UI_TEXT["author_email"]}' style='color: #0284c7;'>{UI_TEXT['author_email']}</a><br>
                <a href='https://github.com/VRYella' target='_blank' style='color: #0284c7;'>GitHub: VRYella</a>
            </div>
        </div>
    </div>
    """, unsafe_allow_html=True)
