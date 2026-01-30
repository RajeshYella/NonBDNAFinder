"""
╔══════════════════════════════════════════════════════════════════════════════╗
║              NATURE-READY VISUALIZATION STANDARDS MODULE                      ║
║          Publication-First Visualization & Redundancy Elimination             ║
╚══════════════════════════════════════════════════════════════════════════════╝

MODULE: visualization_standards.py
AUTHOR: Dr. Venkata Rajesh Yella
VERSION: 2025.1
LICENSE: MIT

NOTE: This module is maintained for backward compatibility.
      The canonical implementation has been moved to visualization/standards.py
      All imports should preferably use 'from visualization import ...'
      
DESCRIPTION:
    Centralized configuration for Nature-ready visualization standards.
    Implements automatic redundancy elimination, plot dominance rules,
    and journal-ready figure panel layouts.

CORE PRINCIPLES:
    1. Eliminate redundant plots automatically
    2. Expose only biologically interpretable metrics
    3. Guarantee non-overlapping, uncluttered figures
    4. Produce journal-ready panels with minimal user tuning
    5. Preserve all advanced outputs for export/supplementary use

REFERENCE:
    Design follows standards from:
    - Nature journal figure guidelines
    - Broad Institute visualization best practices
    - EMBL-EBI data visualization standards
"""

# Re-export everything from the new canonical location
from visualization.standards import (
    # Color schemes
    NATURE_MOTIF_COLORS,
    SUBCLASS_TONE_ADJUSTMENTS,
    
    # Plot management classes
    PlotConcept,
    PlotDominance,
    FigurePanel,
    MetricFilter,
    LabelPolicy,
    UILayout,
    ValidationThresholds,
    
    # Transparency notes
    TRANSPARENCY_NOTE,
    SUPPLEMENTARY_NOTE,
    
    # Helper functions
    get_plot_mapping,
    should_show_plot,
    get_nature_style_params,
)

__all__ = [
    'NATURE_MOTIF_COLORS',
    'SUBCLASS_TONE_ADJUSTMENTS',
    'PlotConcept',
    'PlotDominance',
    'FigurePanel',
    'MetricFilter',
    'LabelPolicy',
    'UILayout',
    'ValidationThresholds',
    'TRANSPARENCY_NOTE',
    'SUPPLEMENTARY_NOTE',
    'get_plot_mapping',
    'should_show_plot',
    'get_nature_style_params',
]
