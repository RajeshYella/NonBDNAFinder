"""
Typography configuration for NBDScanner.

This module contains font settings for the application including:
- Font families (primary and monospace)
- Font sizes for different heading levels and body text
- Font weights from light to extrabold
"""

# ==================== TYPOGRAPHY & FONTS ====================
# Control all font settings for the application - MODERN & READABLE
# Optimized for modern high-resolution displays with excellent readability
FONT_CONFIG = {
    # Primary font families (in order of preference)
    # The browser will use the first available font in the list
    'primary_font': "'Inter', 'IBM Plex Sans', 'Segoe UI', system-ui, -apple-system, sans-serif",
    'monospace_font': "'JetBrains Mono', 'Fira Code', 'Consolas', monospace",
    
    # Font sizes (in rem units, where 1rem ≈ 16px in most browsers)
    # Enhanced sizes for modern, bold, research-quality appearance
    # UPDATED: Larger and bolder for more vibrant and impactful display
    'h1_size': '3.0rem',      # Main page headers - bold, impactful, enlarged
    'h2_size': '2.25rem',     # Section headers - clear hierarchy, prominent
    'h3_size': '1.75rem',     # Subsection headers - organized structure, visible
    'h4_size': '1.35rem',     # Small headers - subtle distinction, readable
    'body_size': '1.0rem',    # Body text, paragraphs - optimal readability
    'small_size': '0.9rem',   # Small text, notes - clear but compact
    'caption_size': '0.8rem', # Captions, footnotes - supporting information
    
    # Font weights (100-900, where 400 is normal and 700 is bold)
    # UPDATED: Heavier weights for bold and vibrant appearance
    'light_weight': 400,       # Light weight (previously 300, now heavier)
    'normal_weight': 600,      # Normal weight (previously 400, now bolder)
    'medium_weight': 700,      # Medium weight (previously 500, now bold)
    'semibold_weight': 800,    # Semibold weight (previously 600, now extra bold)
    'bold_weight': 900,        # Bold weight (previously 700, now maximum)
    'extrabold_weight': 900,   # Extra bold weight (maximum for emphasis)
}
