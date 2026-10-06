"""Small helpers that were duplicated across the analysis scripts."""
import re

from matplotlib.backends.backend_pdf import PdfPages


def read_fasta(fa):
    """Read a FASTA file and return the sequence only, headers dropped."""
    seq = ''
    with open(fa, "r") as fa_file:
        for line in fa_file:
            if not line.startswith('>'):
                seq += line.strip()
    return seq


def parse_attributes(attr_str: str) -> dict:
    """Parse a GTF attributes column into a dictionary."""
    pairs = re.findall(r'(\S+)\s+"?([^";]+)"?;', attr_str)
    return dict(pairs)


def save_figures_to_pdf(output_pdf_path, figures):
    """Write a list of matplotlib figures to one multi-page PDF."""
    with PdfPages(output_pdf_path) as pdf:
        for fig in figures:
            pdf.savefig(fig)
