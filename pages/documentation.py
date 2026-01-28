"""
Upload & Analyze Page
=====================

Single-column, high-density Upload & Analyze UI for NonBDNAFinder.
- Dense table motif selector (deselect-only)
- Inline toggle bar
- Horizontal metrics strip
- Deterministic, high-performance execution
"""

# ─────────────────────────────────────────────────────────────
# Standard library
# ─────────────────────────────────────────────────────────────
import time
import gc
import logging
from typing import List, Dict, Any

# ─────────────────────────────────────────────────────────────
# Third-party
# ─────────────────────────────────────────────────────────────
import streamlit as st
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────
# Project imports (EXPLICIT – no magic)
# ─────────────────────────────────────────────────────────────
from config.text import UI_TEXT
from config.themes import TAB_THEMES

from config.motif_taxonomy import (
    VALID_CLASSES,
    CLASS_TO_SUBCLASSES,
    build_motif_selector_data,
    get_enabled_from_selector_data
)

from ui.css import load_css
from ui.formatters import format_time_scientific

from utilities import (
    parse_fasta,
    parse_fasta_chunked,
    get_file_preview,
    get_basic_stats,
    trigger_garbage_collection,
    get_memory_usage_mb
)

from nonbscanner import analyze_sequence
from job_manager import generate_job_id, save_job_results

# ─────────────────────────────────────────────────────────────
# Logger
# ─────────────────────────────────────────────────────────────
logger = logging.getLogger(__name__)

# ─────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────
EXAMPLE_FASTA = """>Example
GGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
"""

# ─────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────
def ensure_subclass(motif: Dict[str, Any]) -> Dict[str, Any]:
    """Guarantee every motif has Subclass"""
    if "Subclass" not in motif or motif["Subclass"] is None:
        motif["Subclass"] = motif.get("Subtype", "Other")
    return motif


def format_time_compact(seconds: float) -> str:
    m, s = int(seconds // 60), int(seconds % 60)
    return f"{m:02d}:{s:02d}"


# ─────────────────────────────────────────────────────────────
# Main Render Function
# ─────────────────────────────────────────────────────────────
def render() -> None:
    # ---------------------------------------------------------
    # Theme / CSS
    # ---------------------------------------------------------
    load_css(TAB_THEMES.get("Upload & Analyze", "nature_green"))

    # ---------------------------------------------------------
    # Header
    # ---------------------------------------------------------
    st.markdown(
        """
        <div style="text-align:center;margin-bottom:1.5rem">
          <h2 style="
            font-size:2rem;
            font-weight:700;
            background:linear-gradient(135deg,#10b981,#059669);
            -webkit-background-clip:text;
            -webkit-text-fill-color:transparent;">
            Upload & Analyze Sequences
          </h2>
          <p style="color:#64748b">
            High-performance Non-B DNA motif detection
          </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ---------------------------------------------------------
    # Sequence Input
    # ---------------------------------------------------------
    input_method = st.radio(
        UI_TEXT["upload_input_method_prompt"],
        [
            UI_TEXT["upload_method_file"],
            UI_TEXT["upload_method_paste"],
            UI_TEXT["upload_method_example"],
        ],
        horizontal=True,
        label_visibility="collapsed",
    )

    seqs: List[str] = []
    names: List[str] = []

    if input_method == UI_TEXT["upload_method_file"]:
        fasta = st.file_uploader(
            UI_TEXT["upload_file_prompt"],
            type=["fa", "fasta", "fna", "txt"],
        )
        if fasta:
            preview = get_file_preview(fasta, max_sequences=3)
            st.info(
                f"{preview['num_sequences']} sequences · "
                f"{preview['total_bp']:,} bp"
            )
            for name, seq in parse_fasta_chunked(fasta):
                names.append(name)
                seqs.append(seq)

    elif input_method == UI_TEXT["upload_method_paste"]:
        raw = st.text_area("Paste FASTA", height=160)
        if raw:
            parsed = parse_fasta(raw)
            names = list(parsed.keys())
            seqs = list(parsed.values())

    elif input_method == UI_TEXT["upload_method_example"]:
        parsed = parse_fasta(EXAMPLE_FASTA)
        names = list(parsed.keys())
        seqs = list(parsed.values())
        st.success("Loaded example sequence")

    if seqs:
        st.session_state.seqs = seqs
        st.session_state.names = names

    # ---------------------------------------------------------
    # Motif Selector — Dense Table (OPTION 4)
    # ---------------------------------------------------------
    st.markdown("### Motif & Submotif Selection")

    if "motif_selector_data" not in st.session_state:
        st.session_state.motif_selector_data = build_motif_selector_data()

    selector_df = pd.DataFrame(st.session_state.motif_selector_data)

    edited_df = st.data_editor(
        selector_df,
        hide_index=True,
        num_rows="fixed",
        height=280,
        use_container_width=True,
        column_config={
            "Enabled": st.column_config.CheckboxColumn("✓", width="small"),
            "Motif Class": st.column_config.TextColumn("Class", disabled=True),
            "Submotif": st.column_config.TextColumn("Submotif", disabled=True),
        },
    )

    st.session_state.motif_selector_data = edited_df.to_dict("records")

    enabled_classes, enabled_subclasses = get_enabled_from_selector_data(
        st.session_state.motif_selector_data
    )

    st.session_state.selected_classes = enabled_classes
    st.session_state.selected_subclasses = enabled_subclasses

    st.markdown(
        f"""
        <div style="background:#ecfeff;padding:8px 12px;border-radius:6px">
          <strong>{len(enabled_classes)}</strong> classes ·
          <strong>{len(enabled_subclasses)}</strong> submotifs enabled
        </div>
        """,
        unsafe_allow_html=True,
    )

    # ---------------------------------------------------------
    # Inline Toggle Bar
    # ---------------------------------------------------------
    st.markdown("### Analysis Options")
    c1, c2, c3, c4, c5 = st.columns(5)

    with c1:
        detailed = st.checkbox("📊 Detailed", value=True)
    with c2:
        validation = st.checkbox("✅ Validation", value=True)
    with c3:
        parallel = st.checkbox("⚡ Parallel", value=True)
    with c4:
        chunk_view = st.checkbox("📦 Chunks", value=False)
    with c5:
        memory_view = st.checkbox("💾 Memory", value=False)

    # ---------------------------------------------------------
    # Run Analysis
    # ---------------------------------------------------------
    st.markdown("---")

    can_run = bool(seqs and enabled_classes)
    run = st.button(
        "Run Analysis",
        type="primary",
        use_container_width=True,
        disabled=not can_run,
    )

    if not run:
        return

    # ---------------------------------------------------------
    # Analysis Execution (DETERMINISTIC)
    # ---------------------------------------------------------
    job_id = generate_job_id()
    start_time = time.time()

    all_results: List[List[Dict[str, Any]]] = []
    total_bp = 0

    progress = st.progress(0)
    status = st.empty()

    for i, (seq, name) in enumerate(zip(seqs, names)):
        status.info(f"Processing {name} ({len(seq):,} bp)")
        results = analyze_sequence(seq, name)
        results = [ensure_subclass(m) for m in results]

        filtered = [
            m for m in results
            if m.get("Class") in enabled_classes
            and m.get("Subclass") in enabled_subclasses
        ]

        all_results.append(filtered)
        total_bp += len(seq)

        if len(seq) > 1_000_000:
            trigger_garbage_collection()

        progress.progress((i + 1) / len(seqs))

    elapsed = time.time() - start_time
    speed = total_bp / elapsed if elapsed else 0

    st.session_state.results = all_results

    # ---------------------------------------------------------
    # Horizontal Metrics Strip
    # ---------------------------------------------------------
    st.markdown(
        f"""
        <div class="metrics-strip metrics-strip--success">
          <div class="metric-card">
            <span class="metric-card__icon">⏱</span>
            <span class="metric-card__value">{format_time_scientific(elapsed)}</span>
            <span class="metric-card__label">Analysis</span>
          </div>
          <div class="metric-card">
            <span class="metric-card__icon">🧬</span>
            <span class="metric-card__value">{total_bp:,}</span>
            <span class="metric-card__label">Base Pairs</span>
          </div>
          <div class="metric-card">
            <span class="metric-card__icon">⚡</span>
            <span class="metric-card__value">{speed:,.0f}</span>
            <span class="metric-card__label">bp/sec</span>
          </div>
          <div class="metric-card">
            <span class="metric-card__icon">🧪</span>
            <span class="metric-card__value">{sum(len(r) for r in all_results)}</span>
            <span class="metric-card__label">Motifs</span>
          </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

    save_job_results(
        job_id=job_id,
        results=all_results,
        sequences=seqs,
        names=names,
        metadata={
            "analysis_time": elapsed,
            "bp_processed": total_bp,
            "speed_bp_per_sec": speed,
        },
    )

    st.success(f"Analysis complete · Job ID **{job_id}**")
