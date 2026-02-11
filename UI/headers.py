"""
Uniform section heading component for NonBDNAFinder.

This module provides a single, consistent heading style used across all pages
to ensure visual consistency following scientific dashboard design principles.
"""

import streamlit as st
import html

# Import centralized color configuration
from Utilities.config.colors import HOME_COLORS, GLOBAL_COLORS

# Page-specific background color patterns for headings
# Note: 'primary' is used for border color, 'secondary' is the middle gradient stop (for reference),
# and 'pattern' defines the full gradient including a third color for the final stop.
PAGE_HEADING_COLORS = {
    'Home': {'primary': '#0091FF', 'pattern': 'linear-gradient(135deg, #0091FF 0%, #00B4FF 50%, #00D4FF 100%)'},
    'Upload & Analyze': {'primary': '#10b981', 'pattern': 'linear-gradient(135deg, #10b981 0%, #059669 50%, #047857 100%)'},
    'Results': {'primary': '#8b5cf6', 'pattern': 'linear-gradient(135deg, #8b5cf6 0%, #7c3aed 50%, #6d28d9 100%)'},
    'Downloads': {'primary': '#0ea5e9', 'pattern': 'linear-gradient(135deg, #0ea5e9 0%, #0284c7 50%, #0369a1 100%)'},
    'Documentation': {'primary': '#f59e0b', 'pattern': 'linear-gradient(135deg, #f59e0b 0%, #d97706 50%, #b45309 100%)'},
}


def render_section_heading(title: str, page: str = None):
    """Render a uniform section heading with page-specific background colors.
    
    This is the canonical heading component for NonBDNAFinder. All section
    headings should use this function to maintain visual consistency.
    
    Features:
    - Thin box with page-specific gradient background color patterns
    - Increased font size for better readability
    - White glowing text with text-shadow effect
    - Compact, polished appearance
    
    Args:
        title: The heading text to display
        page: Optional page name for page-specific colors (Home, Upload & Analyze, Results, Downloads, Documentation)
    """
    # Escape HTML to prevent XSS
    safe_title = html.escape(title)
    
    # Determine colors based on page or use default
    if page and page in PAGE_HEADING_COLORS:
        colors = PAGE_HEADING_COLORS[page]
        background = colors['pattern']
        border_color = colors['primary']
    else:
        # Default to blue gradient
        primary_color = HOME_COLORS['primary']  # #0091FF
        secondary_color = HOME_COLORS['secondary']  # #00B4FF
        background = f"linear-gradient(135deg, {primary_color} 0%, {secondary_color} 100%)"
        border_color = primary_color
    
    white = GLOBAL_COLORS['white']
    
    st.markdown(f"""
    <div style="
        background: {background};
        padding: 0.6rem 1.5rem;
        border-radius: 6px;
        margin: 0.8rem 0 1rem 0;
        border: 1px solid {border_color};
        text-align: center;
        box-shadow: 0 0 12px rgba(0, 145, 255, 0.3);
        width: 100%;
    ">
        <h2 style="
            margin: 0;
            font-size: 1.5rem;
            font-weight: 700;
            color: {white};
            text-shadow: 0 0 10px rgba(255,255,255,0.8), 0 0 20px rgba(255,255,255,0.6);
            letter-spacing: -0.01em;
        ">
            {safe_title}
        </h2>
    </div>
    """, unsafe_allow_html=True)
