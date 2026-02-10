"""
CSS and styling utilities for NBDScanner.

This module contains functions for:
- Loading CSS with theme variables
- Color conversion utilities
- SVG pattern generation
- Page color retrieval
"""

import os
import streamlit as st

from Utilities.config.colors import (
    GLOBAL_COLORS,
    HOME_COLORS,
    INPUT_COLORS,
    ANALYSIS_COLORS,
    RESULTS_COLORS,
    VISUALIZATION_COLORS,
    DOWNLOAD_COLORS,
    DOCUMENTATION_COLORS,
    SEMANTIC_COLORS,
)
from Utilities.config.themes import COLOR_THEMES
from Utilities.config.typography import FONT_CONFIG
from Utilities.config.layout import LAYOUT_CONFIG
from Utilities.config.animation import ANIMATION_CONFIG


# ─────────────────────────────────────────────────────────────
# Utilities
# ─────────────────────────────────────────────────────────────
def hex_to_rgb(hex_color: str) -> tuple:
    """Convert hex color to RGB tuple for CSS rgba() usage."""
    hex_color = hex_color.lstrip("#")
    return tuple(int(hex_color[i:i + 2], 16) for i in (0, 2, 4))


def get_dna_pattern_svg(stroke_color: str) -> str:
    """Generate subtle DNA helix SVG pattern for background."""
    svg = (
        "%3Csvg xmlns='http://www.w3.org/2000/svg' width='60' height='60' viewBox='0 0 60 60'%3E"
        f"%3Cg fill='none' stroke='%23{stroke_color}' stroke-width='0.6' opacity='0.12'%3E"
        "%3Cpath d='M10 30 C 18 12, 42 12, 50 30'/%3E"
        "%3Cpath d='M10 30 C 18 48, 42 48, 50 30'/%3E"
        "%3C/g%3E%3C/svg%3E"
    )
    return f"url(\"data:image/svg+xml,{svg}\")"


# ─────────────────────────────────────────────────────────────
# CSS Loader
# ─────────────────────────────────────────────────────────────
def load_css(theme_name=None):
    """
    Load external CSS file and inject dynamic theme variables.
    Typography is compact, bold, and vibrant by default.
    """
    theme_to_use = COLOR_THEMES.get(
        theme_name,
        COLOR_THEMES.get(st.session_state.color_theme, COLOR_THEMES["scientific_blue"]),
    )
    is_dark = st.session_state.theme_mode == "dark"

    css_file_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "styles.css")
    try:
        with open(css_file_path, "r") as f:
            css_content = f.read()
    except Exception:
        css_content = ""

    p_rgb = hex_to_rgb(theme_to_use["primary"])
    s_rgb = hex_to_rgb(theme_to_use["secondary"])

    dna_svg_local = get_dna_pattern_svg("1e3a5f" if is_dark else "bbdefb")

    theme_vars = f"""
    <style>
    :root {{
        /* Core colors */
        --primary-color: {theme_to_use['primary']};
        --secondary-color: {theme_to_use['secondary']};
        --accent-color: {theme_to_use['accent']};
        --bg-light: {theme_to_use['bg_light']};
        --bg-card: {theme_to_use['bg_card']};
        --text-color: {theme_to_use['text']};
        --shadow-color: {theme_to_use['shadow']};

        --primary-rgb: {p_rgb[0]}, {p_rgb[1]}, {p_rgb[2]};
        --secondary-rgb: {s_rgb[0]}, {s_rgb[1]}, {s_rgb[2]};

        /* Typography */
        --font-primary: {FONT_CONFIG['primary_font']};
        --font-monospace: {FONT_CONFIG['monospace_font']};

        /* Compact sizing */
        --font-h1: calc({FONT_CONFIG['h1_size']} * 0.95);
        --font-h2: calc({FONT_CONFIG['h2_size']} * 0.95);
        --font-h3: calc({FONT_CONFIG['h3_size']} * 0.96);
        --font-h4: calc({FONT_CONFIG['h4_size']} * 0.97);
        --font-body: calc({FONT_CONFIG['body_size']} * 0.97);
        --font-small: calc({FONT_CONFIG['small_size']} * 0.95);

        /* Heavy weights */
        --font-weight-light: 400;
        --font-weight-normal: 600;
        --font-weight-medium: 700;
        --font-weight-semibold: 800;
        --font-weight-bold: 900;
        --font-weight-extrabold: 950;

        /* Density */
        --line-height-tight: 1.15;
        --line-height-normal: 1.25;
        --letter-spacing-tight: -0.015em;
        --letter-spacing-normal: -0.01em;

        /* Vibrancy */
        --text-glow: 0 0 6px rgba(var(--primary-rgb), 0.25);
        --accent-glow: 0 0 8px rgba(var(--secondary-rgb), 0.35);

        /* Layout */
        --border-radius-sm: {LAYOUT_CONFIG['border_radius']['small']};
        --border-radius-md: {LAYOUT_CONFIG['border_radius']['medium']};
        --border-radius-lg: {LAYOUT_CONFIG['border_radius']['large']};
        --spacing-sm: {LAYOUT_CONFIG['padding']['small']};
        --spacing-md: {LAYOUT_CONFIG['padding']['medium']};
        --spacing-lg: {LAYOUT_CONFIG['padding']['large']};

        /* Motion */
        --transition-fast: {ANIMATION_CONFIG['transition_fast']} {ANIMATION_CONFIG['easing_smooth']};
        --transition-normal: {ANIMATION_CONFIG['transition_normal']} {ANIMATION_CONFIG['easing_smooth']};

        /* Theme state */
        --dark-mode: {1 if is_dark else 0};
        --dna-pattern: {dna_svg_local};
    }}

    /* ───────── Global Typography Overrides ───────── */

    body, .stApp {{
        font-family: var(--font-primary);
        font-weight: var(--font-weight-normal);
        line-height: var(--line-height-normal);
        letter-spacing: var(--letter-spacing-normal);
        color: var(--text-color);
    }}

    h1, h2, h3, h4 {{
        font-weight: var(--font-weight-extrabold);
        letter-spacing: var(--letter-spacing-tight);
        line-height: var(--line-height-tight);
        text-shadow: var(--text-glow);
        margin-bottom: 0.35em;
    }}

    h1 {{ font-size: var(--font-h1); }}
    h2 {{ font-size: var(--font-h2); }}
    h3 {{ font-size: var(--font-h3); }}
    h4 {{ font-size: var(--font-h4); }}

    p, span, label {{
        font-weight: var(--font-weight-medium);
        line-height: var(--line-height-normal);
    }}

    /* Metric cards */
    .stat-card__value {{
        font-weight: var(--font-weight-extrabold);
        letter-spacing: -0.02em;
        text-shadow: var(--accent-glow);
    }}

    /* Tables */
    .stDataFrame td, .stDataFrame th {{
        font-weight: 600;
        font-size: 0.82rem;
        line-height: 1.2;
    }}

    /* Tabs */
    button[data-baseweb="tab"] {{
        font-weight: 800;
        letter-spacing: -0.01em;
    }}

    {css_content}
    </style>
    """

    st.markdown(theme_vars, unsafe_allow_html=True)


# ─────────────────────────────────────────────────────────────
# Page Color Resolver
# ─────────────────────────────────────────────────────────────
def get_page_colors(page_name="Home"):
    """Return literal color values for inline HTML styles."""
    page_color_map = {
        "Home": HOME_COLORS,
        "Upload & Analyze": INPUT_COLORS,
        "Analysis": ANALYSIS_COLORS,
        "Results": RESULTS_COLORS,
        "Visualization": VISUALIZATION_COLORS,
        "Download": DOWNLOAD_COLORS,
        "Documentation": DOCUMENTATION_COLORS,
    }

    page_palette = page_color_map.get(page_name, HOME_COLORS)

    return {
        **page_palette,
        "white": GLOBAL_COLORS["white"],
        "neutral_50": GLOBAL_COLORS["neutral_50"],
        "neutral_100": GLOBAL_COLORS["neutral_100"],
        "neutral_200": GLOBAL_COLORS["neutral_200"],
        "neutral_500": GLOBAL_COLORS["neutral_500"],
        "neutral_600": GLOBAL_COLORS["neutral_600"],
        "neutral_700": GLOBAL_COLORS["neutral_700"],
        "success": SEMANTIC_COLORS["success"],
        "warning": SEMANTIC_COLORS["warning"],
        "error": SEMANTIC_COLORS["error"],
        "info": SEMANTIC_COLORS["info"],
    }
