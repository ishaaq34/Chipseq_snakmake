import streamlit as st
import yaml
import pandas as pd
import os
import subprocess

st.set_page_config(page_title="ChIP-seq Manager", layout="wide", page_icon="🧬")

# --- Constants & Helpers ---
CONFIG_FILE = "config.yaml"

def load_yaml():
    if os.path.exists(CONFIG_FILE):
        with open(CONFIG_FILE, 'r') as f:
            return yaml.safe_load(f)
    return {}

def save_yaml(data):
    with open(CONFIG_FILE, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, sort_keys=False)

def get_metadata_path(config):
    return config.get('macs3', {}).get('metadata_tsv', 'samples_metadata.tsv')

def load_metadata(path):
    if os.path.exists(path):
        return pd.read_csv(path, sep='\t')
    return pd.DataFrame(columns=["raw_id", "condition", "replicate", "assay", "mark_type"])

def save_metadata(path, df):
    df.to_csv(path, sep='\t', index=False)

# Load config
config = load_yaml()

# --- Title ---
st.title("🧬 ChIP-seq Pipeline Manager")
st.markdown("---")

# --- Sidebar ---
with st.sidebar:
    st.header("🚀 Execution")
    if st.button("🧪 Dry Run", use_container_width=True):
        with st.spinner("Checking workflow..."):
            try:
                res = subprocess.run(["snakemake", "-n"], capture_output=True, text=True)
                st.code(res.stdout, language="text")
                if res.stderr: st.error(res.stderr)
            except Exception as e:
                st.error(str(e))
                
    if st.button("▶️ Run Pipeline", type="primary", use_container_width=True):
        threads = config.get('threads', 4)
        try:
            subprocess.Popen(["snakemake", "--cores", str(threads), "--use-conda"])
            st.success("Pipeline started in background!")
        except Exception as e:
            st.error(str(e))
            
    if st.button("📊 Generate Report", use_container_width=True):
        try:
            subprocess.run(["snakemake", "--report", "chipseq_report.zip"])
            st.success("Report created: chipseq_report.zip")
        except Exception as e:
            st.error(str(e))
    
    st.markdown("---")
    st.markdown("[Open Report](chipseq_report.zip)")

# --- Tabs ---
tabs = st.tabs(["🌍 Global", "✂️ Trim/Align", "🏔 Peaks", "📈 Viz", "⚙️ Advanced", "📋 Metadata"])

with tabs[0]:
    st.header("Global Configuration")
    col1, col2 = st.columns(2)
    with col1:
        config['threads'] = st.number_input("Max Threads", 1, 128, config.get('threads', 6))
        config['sequencing_type'] = st.selectbox("Seq Type", ["single", "paired"], index=0 if config.get('sequencing_type')=='single' else 1)
        config['start_from'] = st.selectbox("Start From", ["raw", "fastq", "bam"], index=["raw", "fastq", "bam"].index(config.get('start_from', 'raw')))
    with col2:
        genome = config.setdefault('genome', {})
        genome['fasta'] = st.text_input("Genome FASTA", genome.get('fasta', ''))
        genome['gtf'] = st.text_input("Genome GTF", genome.get('gtf', ''))
        genome['index_prefix'] = st.text_input("Index Prefix", genome.get('index_prefix', 'trim/ssindex'))
        genome['effective_size'] = st.text_input("Genome Eff. Size", str(genome.get('effective_size', '')))
        with st.expander("Auto-Download URLs"):
            genome['fasta_url'] = st.text_input("FASTA URL", genome.get('fasta_url', ''))
            genome['gtf_url'] = st.text_input("GTF URL", genome.get('gtf_url', ''))

with tabs[1]:
    st.header("Trimming & Alignment")
    c1, c2, c3 = st.columns(3)
    with c1:
        st.subheader("Trim Galore")
        tg = config.setdefault('trim_galore', {})
        tg['quality'] = st.slider("Quality", 10, 40, int(tg.get('quality', 20)))
        tg['min_length'] = st.number_input("Min Length", value=int(tg.get('min_length', 20)))
        tg['adapter'] = st.text_input("Adapter", tg.get('adapter', 'auto'))
    with c2:
        st.subheader("Bowtie2")
        bt = config.setdefault('bowtie2', {})
        bt['sensitivity'] = st.selectbox("Sensitivity", ["--very-fast", "--fast", "--sensitive", "--very-sensitive"], index=3)
        bt['no_unal'] = st.checkbox("No Unaligned", bt.get('no_unal', True))
        bt['no_mixed'] = st.checkbox("No Mixed", bt.get('no_mixed', True))
        bt['no_discordant'] = st.checkbox("No Discordant", bt.get('no_discordant', True))
    with c3:
        st.subheader("Samtools")
        st_conf = config.setdefault('samtools', {})
        st_conf['mapq_threshold'] = st.number_input("MapQ", 0, 60, int(st_conf.get('mapq_threshold', 30)))
        st_conf['sam_flags_exclude'] = st.number_input("Exclude Flags", value=int(st_conf.get('sam_flags_exclude', 768)))

with tabs[2]:
    st.header("Peaks & Analysis")
    c1, c2 = st.columns(2)
    with c1:
        st.subheader("MACS3")
        macs = config.setdefault('macs3', {'params': {}})
        mp = macs.get('params', {})
        mp['genome_size'] = st.text_input("Genome Size Code", mp.get('genome_size', 'hs'))
        
        # P/Q Value Toggle
        use_q = st.toggle("Use Q-value instead of P-value", value='qvalue_cutoff' in mp)
        if use_q:
            if 'pvalue_cutoff' in mp: del mp['pvalue_cutoff']
            mp['qvalue_cutoff'] = st.number_input("Q-value", value=float(mp.get('qvalue_cutoff', 0.05)), format="%.4f")
        else:
            if 'qvalue_cutoff' in mp: del mp['qvalue_cutoff']
            mp['pvalue_cutoff'] = st.number_input("P-value", value=float(mp.get('pvalue_cutoff', 0.01)), format="%.4f")
            
        mp['keep_dup'] = st.text_input("Keep Dup", mp.get('keep_dup', 'all'))
        mp['broad_cutoff'] = st.number_input("Broad Cutoff", value=float(mp.get('broad_cutoff', 0.1)))
        mp['format'] = st.selectbox("Format", ["BAM", "BAMPE"], index=0 if mp.get('format')=='BAM' else 1)
        
        with st.expander("Advanced"):
            mp['generate_bedgraph'] = st.checkbox("BedGraph (-B)", mp.get('generate_bedgraph', False))
            mp['bedgraph_spmr'] = st.checkbox("SPMR", mp.get('bedgraph_spmr', False))
            mp['no_model'] = st.checkbox("No Model", mp.get('no_model', False))
            mp['extsize'] = st.number_input("ExtSize", value=int(mp.get('extsize', 200)))
        config['macs3']['params'] = mp

    with c2:
        st.subheader("Optional")
        opt = config.setdefault('optional_tools', {})
        opt['run_spp'] = st.checkbox("Run SPP", opt.get('run_spp', False))
        opt['spp_script'] = st.text_input("SPP Path", opt.get('spp_script', ''))
        opt['run_homer'] = st.checkbox("Run HOMER", opt.get('run_homer', False))
        h = opt.setdefault('homer_params', {})
        h['size'] = st.number_input("Homer Size", value=int(h.get('size', 200)))
        h['length'] = st.text_input("Homer Lengths", h.get('length', '8,10,12'))
        h['mask'] = st.checkbox("Mask Repeats", h.get('mask', True))

with tabs[3]:
    st.header("Visualization (DeepTools)")
    dt = config.setdefault('deeptools', {})
    c1, c2 = st.columns(2)
    with c1:
        dt['normalization_method'] = st.selectbox("Norm Method", ["RPGC", "CPM", "BPM", "RPKM", "None"], index=0)
        dt['bin_size'] = st.number_input("Bin Size", value=int(dt.get('bin_size', 1000)))
        plots = dt.setdefault('plots', {})
        plots['dpi'] = st.number_input("DPI", value=int(plots.get('dpi', 600)))
        
        with st.expander("PCA Styling"):
            pca = plots.setdefault('pca', {})
            pca['plot_width'] = st.number_input("PCA Width", value=float(pca.get('plot_width', 14)))
            pca['plot_height'] = st.number_input("PCA Height", value=float(pca.get('plot_height', 12)))
            # Markers/Colors
            m_str = st.text_input("Markers (comma sep)", value=",".join(pca.get('markers', [])))
            pca['markers'] = [x.strip() for x in m_str.split(',')] if m_str else []
            c_str = st.text_input("Colors (comma sep)", value=",".join(pca.get('colors', [])))
            pca['colors'] = [x.strip() for x in c_str.split(',')] if c_str else []
            
    with c2:
        st.subheader("Region Definitions")
        dt['tss_upstream'] = st.number_input("TSS Up", value=int(dt.get('tss_upstream', 3000)))
        dt['tss_downstream'] = st.number_input("TSS Down", value=int(dt.get('tss_downstream', 10000)))
        dt['gene_body_length'] = st.number_input("Body Len", value=int(dt.get('gene_body_length', 5000)))

with tabs[4]:
    st.header("Advanced")
    dirs = config.setdefault('directories', {})
    for k,v in dirs.items():
        dirs[k] = st.text_input(f"Dir: {k}", v)
    
    st.markdown("---")
    res = config.setdefault('resources', {})
    c1, c2 = st.columns(2)
    with c1:
        res['align_mem_mb'] = st.number_input("Align Mem", value=int(res.get('align_mem_mb', 8000)))
        res['markdup_mem_mb'] = st.number_input("MarkDup Mem", value=int(res.get('markdup_mem_mb', 4000)))
    with c2:
        res['deeptools_mem_mb'] = st.number_input("DeepTools Mem", value=int(res.get('deeptools_mem_mb', 16000)))
        res['default_mem_mb'] = st.number_input("Default Mem", value=int(res.get('default_mem_mb', 4000)))

with tabs[5]:
    st.header("Sample Metadata")
    mp = get_metadata_path(config)
    try:
        df = load_metadata(mp)
        edf = st.data_editor(df, num_rows="dynamic", use_container_width=True)
        if st.button("💾 Save Metadata", type="primary"):
            save_metadata(mp, edf)
            st.success("Metadata saved.")
    except Exception as e:
        st.error(str(e))

st.markdown("---")
if st.button("💾 Save All Config Changes", type="primary", use_container_width=True):
    save_yaml(config)
    st.success("Config saved!")
