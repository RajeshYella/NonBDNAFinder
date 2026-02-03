"""
Uniform section heading component for NonBDNAFinder.

This module provides a single, consistent heading style used across all pages
to ensure visual consistency following scientific dashboard design principles.
"""

import streamlit as st
import html


def render_section_heading(title: str, subtitle: str = None, icon: str = None):
    """Render a uniform section heading used across the app.
    
    This is the canonical heading component for NonBDNAFinder. All section
    headings should use this function to maintain visual consistency.
    
    Features:
    - Elegant gradient accent bar (purple to pink to blue)
    - Consistent typography (1.6rem, weight 700)
    - Optional subtitle with muted styling
    - Optional emoji icon
    - Uniform spacing (1.8rem top, 1.2rem bottom margin)
    
    Args:
        title: The heading text to display
        subtitle: Optional subtitle text (appears below title in muted style)
        icon: Optional emoji icon to display before title
    """
    # Escape HTML to prevent XSS
    safe_title = html.escape(title)
    safe_subtitle = html.escape(subtitle) if subtitle else None
    safe_icon = html.escape(icon) if icon else None
    
    icon_html = f'<span style="margin-right: 8px; font-size: 1.4rem;">{safe_icon}</span>' if safe_icon else ''
    subtitle_html = f'''
        <p style="
            margin: 4px 0 0 0;
            font-size: 0.9rem;
            font-weight: 500;
            color: #64748b;
            letter-spacing: 0.01em;
        ">{safe_subtitle}</p>
    ''' if safe_subtitle else ''
    
    st.markdown(f"""
    <div style="
        display: flex;
        align-items: flex-start;
        gap: 12px;
        margin: 1.8rem 0 1.2rem 0;
        padding: 16px 0;
        border-bottom: 1px solid #e2e8f0;
    ">
        <div style="
            width: 5px;
            min-height: 32px;
            height: 100%;
            background: linear-gradient(180deg, #6366f1 0%, #a855f7 50%, #ec4899 100%);
            border-radius: 4px;
            box-shadow: 0 0 10px rgba(99, 102, 241, 0.3);
        "></div>
        <div style="flex: 1;">
            <h2 style="
                margin: 0;
                font-size: 1.6rem;
                font-weight: 700;
                color: #1e293b;
                display: flex;
                align-items: center;
                letter-spacing: -0.02em;
            ">
                {icon_html}{safe_title}
            </h2>
            {subtitle_html}
        </div>
    </div>
    """, unsafe_allow_html=True)


def render_subsection_heading(title: str, icon: str = None):
    """Render a smaller subsection heading for nested content.
    
    Args:
        title: The heading text to display
        icon: Optional emoji icon
    """
    safe_title = html.escape(title)
    safe_icon = html.escape(icon) if icon else None
    icon_html = f'<span style="margin-right: 6px;">{safe_icon}</span>' if safe_icon else ''
    
    st.markdown(f"""
    <div style="
        display: flex;
        align-items: center;
        gap: 8px;
        margin: 1rem 0 0.75rem 0;
        padding-bottom: 8px;
        border-bottom: 2px solid transparent;
        border-image: linear-gradient(90deg, #6366f1 0%, #a855f7 50%, transparent 100%) 1;
    ">
        <h3 style="
            margin: 0;
            font-size: 1.15rem;
            font-weight: 600;
            color: #334155;
            letter-spacing: -0.01em;
        ">
            {icon_html}{safe_title}
        </h3>
    </div>
    """, unsafe_allow_html=True)
