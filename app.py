"""
NBDScanner - Non-B DNA Motif Finder

ARCHITECTURE:
    Input → Detection → Scoring → Overlap Resolution → Visualization → Export

This is the main application entry point. All configuration, UI utilities, and page
logic have been modularized into separate packages for better maintainability.
"""

import streamlit as st
import pandas as pd
import os
import sys
import logging

# Ensure the current directory is in the Python path for module imports
# This is needed for Streamlit Cloud deployment to find local modules
_current_dir = os.path.dirname(os.path.abspath(__file__))
if _current_dir not in sys.path:
    sys.path.insert(0, _current_dir)

# Import configuration modules
from config.text import UI_TEXT
from config.layout import LAYOUT_CONFIG
from config.themes import TAB_THEMES

# Import UI utilities
from ui.css import load_css

# Import page modules
from pages import home, upload, results, download, documentation

# Import core scanner functionality
from nonbscanner import get_motif_info as get_motif_classification_info

# Optional imports with fallbacks
try:
    from Bio import Entrez, SeqIO
    BIO_AVAILABLE = True
except ImportError:
    BIO_AVAILABLE = False

# Try to import Hyperscan (optional for acceleration)
logger = logging.getLogger(__name__)

HYPERSCAN_AVAILABLE = False
HYPERSCAN_VERSION = None
HYPERSCAN_ERROR = None

try:
    import hyperscan
    HYPERSCAN_AVAILABLE = True
    try:
        HYPERSCAN_VERSION = hyperscan.__version__
    except AttributeError:
        HYPERSCAN_VERSION = 'unknown'
    logger.info(f"Hyperscan loaded successfully (version: {HYPERSCAN_VERSION})")
except ImportError as e:
    HYPERSCAN_ERROR = f"Hyperscan Python bindings not installed: {e}"
    logger.info(f"Hyperscan not available - using pure Python fallback. {HYPERSCAN_ERROR}")
except Exception as e:
    HYPERSCAN_ERROR = f"Hyperscan initialization failed: {e}"
    logger.warning(f"Hyperscan not available - using pure Python fallback. {HYPERSCAN_ERROR}")

# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title=f"{UI_TEXT['app_title']} - Non-B DNA Motif Finder",
    layout=LAYOUT_CONFIG['layout_mode'],
    page_icon=None,
    menu_items={'About': f"NBDScanner | Developed by {UI_TEXT['author']}"}
)

# Get motif classification info
CLASSIFICATION_INFO = get_motif_classification_info()

# Configure Entrez if available
if BIO_AVAILABLE:
    Entrez.email = "raazbiochem@gmail.com"
    Entrez.api_key = None

# Define pages
PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Sequence Upload and Motif Analysis",
    "Results": "Analysis Results and Visualization",
    "Download": "Export Data",
    "Documentation": "Scientific Documentation & References"
}

# Initialize theme state in session state
if 'theme_mode' not in st.session_state:
    st.session_state.theme_mode = 'light'
if 'table_density' not in st.session_state:
    st.session_state.table_density = 'relaxed'
if 'color_theme' not in st.session_state:
    st.session_state.color_theme = 'scientific_blue'

# Streamlined session state initialization
for k, v in {
    'seqs': [],
    'names': [],
    'results': [],
    'summary_df': pd.DataFrame(),
    'analysis_status': "Ready",
    'selected_classes': [],  # Initialize empty list for motif class selection
    'current_job_id': None  # Current job ID for result delivery
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

# ---- APPLICATION HEADER ----
st.markdown("""
<div style='text-align: center; padding: 1.5rem 0 1rem 0; margin-bottom: 1.5rem;
            border-bottom: 2px solid #e2e8f0;'>
    <h1 style='margin: 0; font-size: 2.5rem; background: linear-gradient(135deg, #4A90E2, #764ba2);
               -webkit-background-clip: text; -webkit-text-fill-color: transparent;
               font-weight: 800; letter-spacing: -0.02em;'>
        🧬 NonBDNAFinder
    </h1>
    <p style='margin: 0.5rem 0 0 0; color: #64748b; font-size: 1rem;'>
        Advanced Non-B DNA Motif Detection & Analysis Platform
    </p>
</div>
""", unsafe_allow_html=True)

# ---- TABS WITH ICONS ----
tab_icons = ["🏠", "📤", "📊", "💾", "📚"]
tab_labels = [f"{icon} {label}" for icon, label in zip(tab_icons, PAGES.keys())]
tabs = st.tabs(tab_labels)
tab_pages = dict(zip(PAGES.keys(), tabs))

# Render each page in its respective tab
with tab_pages["Home"]:
    home.render()

with tab_pages["Upload & Analyze"]:
    upload.render()

with tab_pages["Results"]:
    results.render()

with tab_pages["Download"]:
    download.render()

with tab_pages["Documentation"]:
    documentation.render()
