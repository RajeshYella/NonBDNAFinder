"""
Uniform section heading component for NonBDNAFinder.

This module provides a single, consistent heading style used across all pages
to ensure visual consistency following scientific dashboard design principles.
"""

import streamlit as st
import html


def render_section_heading(title: str):
    """Render a uniform section heading used across the app.
    
    This is the canonical heading component for NonBDNAFinder. All section
    headings should use this function to maintain visual consistency.
    
    Features:
    - Gradient accent bar (cyan to green)
    - Consistent typography (1.6rem, weight 700)
    - Uniform spacing (1.8rem top, 1.2rem bottom margin)
    - No captions or subtitles
    
    Args:
        title: The heading text to display
    """
    # Escape HTML to prevent XSS
    safe_title = html.escape(title)
    
    st.markdown(f"""
    <div style="
        display: flex;
        align-items: center;
        gap: 10px;
        margin: 1.8rem 0 1.2rem 0;
    ">
        <div style="
            width: 4px;
            height: 28px;
            background: linear-gradient(180deg, #0ea5e9, #22c55e);
            border-radius: 4px;
        "></div>
        <h2 style="
            margin: 0;
            font-size: 1.6rem;
            font-weight: 700;
            color: #1f2937;
        ">
            {safe_title}
        </h2>
    </div>
    """, unsafe_allow_html=True)
