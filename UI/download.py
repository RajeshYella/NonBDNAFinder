"""
┌──────────────────────────────────────────────────────────────────────────────┐
│ Download Page - Export Data and Results                                      │
├──────────────────────────────────────────────────────────────────────────────┤
│ Author: Dr. Venkata Rajesh Yella | License: MIT | Version: 2024.1            │
└──────────────────────────────────────────────────────────────────────────────┘
"""
# ═══════════════════════════════════════════════════════════════════════════════
# IMPORTS
# ═══════════════════════════════════════════════════════════════════════════════
import streamlit as st; import pandas as pd; import numpy as np; import re; import io; import traceback
from collections import Counter
from Utilities.config.text import UI_TEXT; from Utilities.config.themes import TAB_THEMES
from UI.css import load_css; from UI.headers import render_section_heading; from UI.guards import generate_excel_bytes
from Utilities.utilities import export_to_csv, export_to_json, export_to_excel, export_to_pdf, export_to_bed
from Utilities.export.export_validator import validate_export_data

# ═══════════════════════════════════════════════════════════════════════════════
# TUNABLE PARAMETERS
# ═══════════════════════════════════════════════════════════════════════════════
FILENAME_MAX_LENGTH = 50
# ═══════════════════════════════════════════════════════════════════════════════

def render():
    load_css(TAB_THEMES.get('Download', 'clinical_teal')); render_section_heading("Download & Export Results", page="Downloads")
    if not st.session_state.results: st.info(UI_TEXT['download_no_results']); st.markdown("<div style='background:linear-gradient(135deg,#f0f9ff 0%,#e0f2fe 100%);padding:1.2rem;border-radius:12px;margin-top:0.8rem;border:1px solid #bae6fd;text-align:center;'><h3 style='color:#0284c7;margin:0 0 0.6rem 0;'>Export Formats Available</h3><p style='color:#6b7280;margin:0 0 0.6rem 0;'>Once you analyze a sequence, export results in:</p><div style='display:flex;justify-content:center;gap:0.6rem;flex-wrap:wrap;'><span style='background:#0ea5e9;color:white;padding:0.35rem 0.8rem;border-radius:8px;font-weight:600;'>CSV</span><span style='background:#0ea5e9;color:white;padding:0.35rem 0.8rem;border-radius:8px;font-weight:600;'>Excel</span><span style='background:#0ea5e9;color:white;padding:0.35rem 0.8rem;border-radius:8px;font-weight:600;'>JSON</span><span style='background:#0ea5e9;color:white;padding:0.35rem 0.8rem;border-radius:8px;font-weight:600;'>BED</span><span style='background:#0ea5e9;color:white;padding:0.35rem 0.8rem;border-radius:8px;font-weight:600;'>PDF</span></div></div>", unsafe_allow_html=True); return
    seq_name = st.session_state.names[0] if st.session_state.names else "Unknown"; safe_fn = re.sub(r'[^\w\-]', '_', seq_name)[:FILENAME_MAX_LENGTH].strip('_')
    cls_used = st.session_state.get('selected_classes_used', []); sub_used = st.session_state.get('selected_subclasses_used', []); mode_used = st.session_state.get('analysis_mode_used', 'Motif Level')
    if cls_used: st.markdown(f"<div style='background:linear-gradient(135deg,#fefce8 0%,#fef9c3 100%);padding:0.4rem;border-radius:6px;margin-bottom:0.5rem;border-left:3px solid #ca8a04;'><p style='color:#713f12;margin:0;font-size:0.8rem;'><strong>Data Filter:</strong> Exports contain selected classes ({len(cls_used)}) and subclasses ({len(sub_used)}) from {mode_used} analysis</p></div>", unsafe_allow_html=True)
    all_motifs = []
    for i, motifs in enumerate(st.session_state.results):
        for m in motifs: em = m.copy(); em.setdefault('Sequence_Name', st.session_state.names[i]); all_motifs.append(em)
    try: all_motifs = validate_export_data(all_motifs, auto_normalize=True, strict=False)
    except Exception as e: st.warning(f"Some motifs normalized: {e}")
    st.markdown("### Export Options"); st.markdown("<div style='background:linear-gradient(135deg,#f0f9ff 0%,#e0f2fe 100%);padding:0.6rem;border-radius:10px;margin-bottom:0.8rem;border-left:4px solid #0ea5e9;'><p style='color:#0c4a6e;margin:0;font-size:0.8rem;'><strong>Quick Export:</strong> Choose your preferred format.</p></div>", unsafe_allow_html=True)
    st.markdown("#### Standard Formats"); c1, c2, c3, c4, c5 = st.columns(5)
    with c1:
        if all_motifs: st.download_button("CSV (All Motifs)", data=export_to_csv(all_motifs, non_overlapping_only=False).encode('utf-8'), file_name=f"{safe_fn}_all_motifs.csv", mime="text/csv", use_container_width=True, type="primary", help="CSV with all motifs")
    with c2:
        if all_motifs:
            try: st.download_button("Excel (2 tabs)", data=generate_excel_bytes(all_motifs, simple_format=True), file_name=f"{safe_fn}_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True, type="primary", help="Excel workbook")
            except Exception as e: st.error(f"Excel error: {e}")
    with c3:
        if all_motifs: st.download_button("JSON", data=export_to_json(all_motifs, pretty=True).encode('utf-8'), file_name=f"{safe_fn}_results.json", mime="application/json", use_container_width=True, type="primary")
    with c4:
        if all_motifs and st.session_state.names: st.download_button("BED", data=export_to_bed(all_motifs, st.session_state.names[0]).encode('utf-8'), file_name=f"{safe_fn}_results.bed", mime="text/plain", use_container_width=True, type="primary")
    with c5:
        if all_motifs and st.session_state.seqs:
            try:
                slen = len(st.session_state.seqs[0]) if st.session_state.seqs and st.session_state.seqs[0] else 0
                if slen > 0: st.download_button("PDF", data=export_to_pdf(all_motifs, slen, seq_name), file_name=f"{safe_fn}_viz.pdf", mime="application/pdf", use_container_width=True, type="primary")
                else: st.warning("No sequence for PDF")
            except Exception as e: st.error(f"PDF error: {e}")
    st.markdown("---"); st.markdown("### Statistical Analysis Tables"); st.markdown("<div style='background:linear-gradient(135deg,#f0fdf4 0%,#dcfce7 100%);padding:0.6rem;border-radius:10px;margin-bottom:0.8rem;border-left:4px solid #22c55e;'><p style='color:#14532d;margin:0;font-size:0.8rem;'><strong>Advanced Analytics:</strong> Detailed distribution and density statistics</p></div>", unsafe_allow_html=True)
    if all_motifs and st.session_state.seqs:
        try:
            dist_data = []
            for seq, name, motifs in zip(st.session_state.seqs, st.session_state.names, st.session_state.results):
                slen = len(seq); cls_cnt = Counter(m.get('Class', 'Unknown') for m in motifs)
                for cname, cnt in cls_cnt.items():
                    gd = (sum(m.get('Length', 0) for m in motifs if m.get('Class') == cname) / slen * 100) if slen > 0 else 0
                    mkb = (cnt / slen * 1000) if slen > 0 else 0; avl = np.mean([m.get('Length', 0) for m in motifs if m.get('Class') == cname])
                    dist_data.append({'Sequence Name': name, 'Motif Class': cname.replace('_', ' '), 'Count': cnt, 'Genomic Density (%)': f"{gd:.4f}", 'Motifs per kbp': f"{mkb:.2f}", 'Average Length (bp)': f"{avl:.1f}", 'Total Coverage (bp)': sum(m.get('Length', 0) for m in motifs if m.get('Class') == cname)})
            dist_df = pd.DataFrame(dist_data)
            sub_data = []
            for seq, name, motifs in zip(st.session_state.seqs, st.session_state.names, st.session_state.results):
                slen = len(seq); sub_cnt = Counter(m.get('Subclass', 'Unknown') for m in motifs)
                for sname, cnt in sub_cnt.items():
                    pcls = next((m.get('Class') for m in motifs if m.get('Subclass') == sname), 'Unknown')
                    gd = (sum(m.get('Length', 0) for m in motifs if m.get('Subclass') == sname) / slen * 100) if slen > 0 else 0
                    mkb = (cnt / slen * 1000) if slen > 0 else 0; avl = np.mean([m.get('Length', 0) for m in motifs if m.get('Subclass') == sname])
                    sub_data.append({'Sequence Name': name, 'Motif Class': pcls.replace('_', ' '), 'Motif Subclass': sname.replace('_', ' '), 'Count': cnt, 'Genomic Density (%)': f"{gd:.4f}", 'Motifs per kbp': f"{mkb:.2f}", 'Average Length (bp)': f"{avl:.1f}", 'Total Coverage (bp)': sum(m.get('Length', 0) for m in motifs if m.get('Subclass') == sname)})
            sub_df = pd.DataFrame(sub_data)
            st.markdown("#### Class-Level Distribution Statistics"); st.dataframe(dist_df.head(10), use_container_width=True, height=300); st.caption(f"Showing first 10 of {len(dist_df)} records")
            st.markdown("#### Subclass-Level Distribution Statistics"); st.dataframe(sub_df.head(10), use_container_width=True, height=300); st.caption(f"Showing first 10 of {len(sub_df)} records")
            s1, s2, s3 = st.columns(3)
            with s1: st.download_button("Class Statistics (CSV)", data=dist_df.to_csv(index=False).encode('utf-8'), file_name=f"{safe_fn}_class_statistics.csv", mime="text/csv", use_container_width=True)
            with s2: st.download_button("Subclass Statistics (CSV)", data=sub_df.to_csv(index=False).encode('utf-8'), file_name=f"{safe_fn}_subclass_statistics.csv", mime="text/csv", use_container_width=True)
            with s3:
                try:
                    out = io.BytesIO()
                    with pd.ExcelWriter(out, engine='openpyxl') as w: dist_df.to_excel(w, sheet_name='Class Statistics', index=False); sub_df.to_excel(w, sheet_name='Subclass Statistics', index=False)
                    out.seek(0); st.download_button("All Statistics (Excel)", data=out.getvalue(), file_name=f"{safe_fn}_all_statistics.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", use_container_width=True)
                except Exception as e: st.error(f"Excel error: {e}")
        except Exception as e: st.error(f"Error generating statistics: {e}"); st.code(traceback.format_exc(), language="python")
