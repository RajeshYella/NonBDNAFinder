"""
Uniform section heading component for NonBDNAFinder.

This module provides a single, consistent heading style used across all pages
to ensure visual consistency following scientific dashboard design principles.
"""

import streamlit as st
import html

# Import centralized color configuration
from Utilities.config.colors import HOME_COLORS, GLOBAL_COLORS


def render_section_heading(title: str):
    """Render a uniform section heading as a thin blue box with white glowing text.
    
    This is the canonical heading component for NonBDNAFinder. All section
    headings should use this function to maintain visual consistency.
    
    Features:
    - Thin blue gradient box using centralized color tokens
    - White glowing text with text-shadow effect
    - Compact, polished appearance
    
    Args:
        title: The heading text to display
    """
    # Escape HTML to prevent XSS
    safe_title = html.escape(title)
    
    # Use centralized colors from HOME_COLORS for consistency
    primary_color = HOME_COLORS['primary']  # #0091FF
    secondary_color = HOME_COLORS['secondary']  # #00B4FF
    white = GLOBAL_COLORS['white']
    
    st.markdown(f"""
    <div style="
        background: linear-gradient(135deg, {primary_color} 0%, {secondary_color} 100%);
        padding: 0.75rem 1.5rem;
        border-radius: 8px;
        margin: 1rem 0 1.2rem 0;
        border: 2px solid {primary_color};
        text-align: center;
        box-shadow: 0 0 15px rgba(0, 145, 255, 0.4);
    ">
        <h2 style="
            margin: 0;
            font-size: 1.3rem;
            font-weight: 600;
            color: {white};
            text-shadow: 0 0 10px rgba(255,255,255,0.8), 0 0 20px rgba(255,255,255,0.6);
        ">
            {safe_title}
        </h2>
    </div>
    """, unsafe_allow_html=True)
