"""
Uniform section heading component for NonBDNAFinder.

This module provides a single, consistent heading style used across all pages
to ensure visual consistency following scientific dashboard design principles.
"""

import streamlit as st
import html


def render_section_heading(title: str):
    """Render a uniform section heading as a thin blue box with white glowing text.
    
    This is the canonical heading component for NonBDNAFinder. All section
    headings should use this function to maintain visual consistency.
    
    Features:
    - Thin blue gradient box
    - White glowing text with text-shadow effect
    - Compact, polished appearance
    
    Args:
        title: The heading text to display
    """
    # Escape HTML to prevent XSS
    safe_title = html.escape(title)
    
    st.markdown(f"""
    <div style="
        background: linear-gradient(135deg, #1e40af 0%, #3b82f6 100%);
        padding: 0.75rem 1.5rem;
        border-radius: 8px;
        margin: 1rem 0 1.2rem 0;
        border: 2px solid #1e40af;
        text-align: center;
        box-shadow: 0 0 15px rgba(30, 64, 175, 0.4);
    ">
        <h2 style="
            margin: 0;
            font-size: 1.3rem;
            font-weight: 600;
            color: #ffffff;
            text-shadow: 0 0 10px rgba(255,255,255,0.8), 0 0 20px rgba(255,255,255,0.6);
        ">
            {safe_title}
        </h2>
    </div>
    """, unsafe_allow_html=True)
