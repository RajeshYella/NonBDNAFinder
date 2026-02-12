"""Download tab content for NonBDNAFinder"""

import streamlit as st
import pandas as pd
import numpy as np
import re
import io
import traceback
from collections import Counter

from Utilities.config.text import UI_TEXT
from Utilities.config.themes import TAB_THEMES
from UI.css import load_css
from UI.headers import render_section_heading
from UI.guards import generate_excel_bytes
from Utilities.utilities import (
    export_to_csv,
    export_to_json,
    export_to_excel,
    export_to_pdf,
    export_to_bed
)
from Utilities.export.export_validator import validate_export_data


def render():
    """Render the Download tab content"""
    # Apply Download tab theme based on configuration
    load_css(TAB_THEMES.get('Download', 'clinical_teal'))
    
    # Uniform section heading with page-specific color
    render_section_heading("Download & Export Results", page="Downloads")
    
    if not st.session_state.results:
        st.info(UI_TEXT['download_no_results'])
        
        # Show placeholder content explaining available export formats
        st.markdown("""
        <div style='background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%); 
                    padding: 1.2rem; border-radius: 12px; margin-top: 0.8rem;
                    border: 1px solid #bae6fd; text-align: center;'>
            <h3 style='color: #0284c7; margin: 0 0 0.6rem 0;'>Export Formats Available</h3>
            <p style='color: #6b7280; margin: 0 0 0.6rem 0;'>
                Once you upload and analyze a sequence, you can export results in the following formats:
            </p>
            <div style='display: flex; justify-content: center; gap: 0.6rem; flex-wrap: wrap;'>
                <span style='background: #0ea5e9; color: white; padding: 0.35rem 0.8rem; border-radius: 8px; font-weight: 600;'>CSV</span>
                <span style='background: #0ea5e9; color: white; padding: 0.35rem 0.8rem; border-radius: 8px; font-weight: 600;'>Excel</span>
                <span style='background: #0ea5e9; color: white; padding: 0.35rem 0.8rem; border-radius: 8px; font-weight: 600;'>JSON</span>
                <span style='background: #0ea5e9; color: white; padding: 0.35rem 0.8rem; border-radius: 8px; font-weight: 600;'>BED</span>
                <span style='background: #0ea5e9; color: white; padding: 0.35rem 0.8rem; border-radius: 8px; font-weight: 600;'>PDF</span>
            </div>
        </div>
        """, unsafe_allow_html=True)
        return
    
    primary_sequence_name = st.session_state.names[0] if st.session_state.names else "Unknown Sequence"
    analysis_time = st.session_state.get('analysis_time', 0)
    
    # Sanitize sequence name for safe filenames
    safe_filename = re.sub(r'[^\w\-]', '_', primary_sequence_name)[:50].strip('_')
    
    # Show selected class/subclass filter information (compact)
    selected_classes_used = st.session_state.get('selected_classes_used', [])
    selected_subclasses_used = st.session_state.get('selected_subclasses_used', [])
    analysis_mode_used = st.session_state.get('analysis_mode_used', 'Motif Level')
    
    if selected_classes_used:
        st.markdown(f"""
        <div style='background: linear-gradient(135deg, #fefce8 0%, #fef9c3 100%); 
                    padding: 0.4rem; border-radius: 6px; margin-bottom: 0.5rem;
                    border-left: 3px solid #ca8a04;'>
            <p style='color: #713f12; margin: 0; font-size: 0.8rem;'>
                <strong>Data Filter:</strong> Exports contain only selected classes ({len(selected_classes_used)}) and subclasses ({len(selected_subclasses_used)}) 
                from {analysis_mode_used} analysis
            </p>
        </div>
        """, unsafe_allow_html=True)
    
    # Prepare motif data and validate
    all_motifs = []
    for i, motifs in enumerate(st.session_state.results):
        for m in motifs:
            export_motif = m.copy()
            if 'Sequence_Name' not in export_motif:
                export_motif['Sequence_Name'] = st.session_state.names[i]
            all_motifs.append(export_motif)
    
    # Validate motifs before export (auto-normalize if needed)
    try:
        all_motifs = validate_export_data(all_motifs, auto_normalize=True, strict=False)
    except Exception as e:
        st.warning(f"Some motifs had invalid class/subclass data and were normalized: {e}")
    
    # Individual file downloads as main option (compact)
    st.markdown("### Export Options")
    st.markdown("""
    <div style='background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%); 
                padding: 0.6rem; border-radius: 10px; margin-bottom: 0.8rem;
                border-left: 4px solid #0ea5e9;'>
        <p style='color: #0c4a6e; margin: 0; font-size: 0.8rem;'>
            <strong>Quick Export:</strong> Choose your preferred format below. All exports include complete motif data.
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Enhanced download section with better visual hierarchy
    st.markdown("#### Standard Formats")
    col1, col2, col3, col4, col5 = st.columns(5)
    
    with col1:
        # CSV Export (All Motifs)
        if all_motifs:
            csv_data = export_to_csv(all_motifs, non_overlapping_only=False)
            st.download_button(
                "CSV (All Motifs)", 
                data=csv_data.encode('utf-8'), 
                file_name=f"{safe_filename}_all_motifs.csv", 
                mime="text/csv",
                use_container_width=True,
                type="primary",
                help="Download CSV with all detected motifs including Hybrid and Clusters"
            )
    
    with col2:
        # Excel Export (2-tab format)
        if all_motifs:
            try:
                excel_bytes = generate_excel_bytes(all_motifs, simple_format=True)
                st.download_button(
                    "Excel (2 tabs)", 
                    data=excel_bytes, 
                    file_name=f"{safe_filename}_results.xlsx", 
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    use_container_width=True,
                    type="primary",
                    help="Download Excel with 2 tabs: NonOverlappingConsolidated, OverlappingAll"
                )
            except Exception as e:
                st.error(f"Excel export error: {str(e)}")
    
    with col3:
        # JSON Export  
        if all_motifs:
            json_data = export_to_json(all_motifs, pretty=True)
            st.download_button(
                "JSON", 
                data=json_data.encode('utf-8'), 
                file_name=f"{safe_filename}_results.json", 
                mime="application/json",
                use_container_width=True,
                type="primary",
                help="Download results in JSON format"
            )
    
    with col4:
        # BED Export
        if all_motifs and st.session_state.names:
            bed_data = export_to_bed(all_motifs, st.session_state.names[0])
            st.download_button(
                "BED", 
                data=bed_data.encode('utf-8'), 
                file_name=f"{safe_filename}_results.bed", 
                mime="text/plain",
                use_container_width=True,
                type="primary",
                help="Download results in BED format for genome browsers"
            )
    
    with col5:
        # PDF Export for visual summaries
        if all_motifs and st.session_state.seqs:
            try:
                # Get sequence length from the first sequence (with safe None check)
                sequence_length = len(st.session_state.seqs[0]) if st.session_state.seqs and st.session_state.seqs[0] else 0
                if sequence_length > 0:
                    pdf_data = export_to_pdf(all_motifs, sequence_length, primary_sequence_name)
                    st.download_button(
                        "PDF (Visualizations)",
                        data=pdf_data,
                        file_name=f"{safe_filename}_visualizations.pdf",
                        mime="application/pdf",
                        use_container_width=True,
                        type="primary",
                        help="Download PDF containing graphical summarizations of sequence analyses"
                    )
                else:
                    st.warning("No sequence data available for PDF generation")
            except Exception as e:
                st.error(f"PDF export error: {str(e)}")
    
    # Add Distribution & Statistics Tables Download Section (compact)
    st.markdown("---")
    st.markdown("### Statistical Analysis Tables")
    st.markdown("""
    <div style='background: linear-gradient(135deg, #f0fdf4 0%, #dcfce7 100%); 
                padding: 0.6rem; border-radius: 10px; margin-bottom: 0.8rem;
                border-left: 4px solid #22c55e;'>
        <p style='color: #14532d; margin: 0; font-size: 0.8rem;'>
            <strong>Advanced Analytics:</strong> Detailed distribution and density statistics 
            for publication-quality analysis and reporting
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Calculate distribution and density statistics for all sequences
    if all_motifs and st.session_state.seqs:
        try:
            # Prepare distribution statistics
            distribution_data = []
            for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, st.session_state.results)):
                sequence_length = len(seq)
                
                # Show all motifs including hybrid/cluster motifs
                # No filtering is applied - all results are included in statistics
                filtered_motifs = motifs
                
                # Calculate class-level statistics
                class_counts = Counter(m.get('Class', 'Unknown') for m in filtered_motifs)
                for class_name, count in class_counts.items():
                    genomic_density = (sum(m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name) / sequence_length * 100) if sequence_length > 0 else 0
                    motifs_per_kbp = (count / sequence_length * 1000) if sequence_length > 0 else 0
                    avg_length = np.mean([m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name])
                    
                    distribution_data.append({
                        'Sequence Name': name,
                        'Motif Class': class_name.replace('_', ' '),
                        'Count': count,
                        'Genomic Density (%)': f"{genomic_density:.4f}",
                        'Motifs per kbp': f"{motifs_per_kbp:.2f}",
                        'Average Length (bp)': f"{avg_length:.1f}",
                        'Total Coverage (bp)': sum(m.get('Length', 0) for m in filtered_motifs if m.get('Class') == class_name)
                    })
            
            distribution_df = pd.DataFrame(distribution_data)
            
            # Prepare subclass-level statistics
            subclass_data = []
            for seq_idx, (seq, name, motifs) in enumerate(zip(st.session_state.seqs, st.session_state.names, st.session_state.results)):
                sequence_length = len(seq)
                
                # Show all motifs including hybrid/cluster motifs
                # No filtering is applied - all results are included in statistics
                filtered_motifs = motifs
                
                # Calculate subclass-level statistics
                subclass_counts = Counter(m.get('Subclass', 'Unknown') for m in filtered_motifs)
                for subclass_name, count in subclass_counts.items():
                    parent_class = next((m.get('Class') for m in filtered_motifs if m.get('Subclass') == subclass_name), 'Unknown')
                    genomic_density = (sum(m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name) / sequence_length * 100) if sequence_length > 0 else 0
                    motifs_per_kbp = (count / sequence_length * 1000) if sequence_length > 0 else 0
                    avg_length = np.mean([m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name])
                    
                    subclass_data.append({
                        'Sequence Name': name,
                        'Motif Class': parent_class.replace('_', ' '),
                        'Motif Subclass': subclass_name.replace('_', ' '),
                        'Count': count,
                        'Genomic Density (%)': f"{genomic_density:.4f}",
                        'Motifs per kbp': f"{motifs_per_kbp:.2f}",
                        'Average Length (bp)': f"{avg_length:.1f}",
                        'Total Coverage (bp)': sum(m.get('Length', 0) for m in filtered_motifs if m.get('Subclass') == subclass_name)
                    })
            
            subclass_df = pd.DataFrame(subclass_data)
            
            # Display preview of tables
            st.markdown("#### Class-Level Distribution Statistics")
            st.dataframe(distribution_df.head(10), use_container_width=True, height=300)
            st.caption(f"Showing first 10 of {len(distribution_df)} total records")
            
            st.markdown("#### Subclass-Level Distribution Statistics")
            st.dataframe(subclass_df.head(10), use_container_width=True, height=300)
            st.caption(f"Showing first 10 of {len(subclass_df)} total records")
            
            # Download buttons for statistics tables
            col_stat1, col_stat2, col_stat3 = st.columns(3)
            
            with col_stat1:
                # Class-level CSV
                class_csv = distribution_df.to_csv(index=False)
                st.download_button(
                    "Class Statistics (CSV)",
                    data=class_csv.encode('utf-8'),
                    file_name=f"{safe_filename}_class_statistics.csv",
                    mime="text/csv",
                    use_container_width=True,
                    help="Download class-level distribution statistics"
                )
            
            with col_stat2:
                # Subclass-level CSV
                subclass_csv = subclass_df.to_csv(index=False)
                st.download_button(
                    "Subclass Statistics (CSV)",
                    data=subclass_csv.encode('utf-8'),
                    file_name=f"{safe_filename}_subclass_statistics.csv",
                    mime="text/csv",
                    use_container_width=True,
                    help="Download subclass-level distribution statistics"
                )
            
            with col_stat3:
                # Combined Excel with both sheets
                try:
                    output = io.BytesIO()
                    with pd.ExcelWriter(output, engine='openpyxl') as writer:
                        distribution_df.to_excel(writer, sheet_name='Class Statistics', index=False)
                        subclass_df.to_excel(writer, sheet_name='Subclass Statistics', index=False)
                    
                    output.seek(0)
                    st.download_button(
                        "All Statistics (Excel)",
                        data=output.getvalue(),
                        file_name=f"{safe_filename}_all_statistics.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                        use_container_width=True,
                        help="Download all statistics in Excel format with separate sheets"
                    )
                except Exception as e:
                    st.error(f"Excel generation error: {str(e)}")
            
        except Exception as e:
            st.error(f"Error generating distribution statistics: {str(e)}")
            st.code(traceback.format_exc(), language="python")
