"""
Count the spectroscopic redshift sources per facility and write a LaTeX table.
"""
from pathlib import Path

from astropy.table import Table
# -------------------------------------------------------------------------------------------- #
specz_dir = Path(__file__).resolve().parent.parent
specz = Table.read(specz_dir / '2026-11-02.MOO_1142+1527.speczs.all.fits')
source = [str(s).strip() for s in specz['Source']]


def count(*names):
    return sum(s in names for s in source)


n_jwst = sum(s.startswith('J') for s in source)
n_internal = sum(not s.startswith('J') for s in source)

# Survey row and N_members are not filled in yet
rows = [
    ('JWST', n_jwst, f'JP ({count("JP")})', f'JG ({count("JG")})', ''),
    ('Internal', n_internal, f'KL/r ({count("KL", "KLr")})', f'KG ({count("KG")})', ''),
    ('Survey', '', '', '', ''),
]

lines = [r'\begin{tabular}{lcccc}',
         r'\hline',
         r' & $N_{\mathrm{tot}}$ & Facility 1 & Facility 2 & $N_{\mathrm{members}}$ \\',
         r'\hline']
lines += [' & '.join(str(c) for c in row) + r' \\' for row in rows]
lines += [r'\hline', r'\end{tabular}']
table = '\n'.join(lines) + '\n'

out_dir = specz_dir / 'latex_tables'
out_dir.mkdir(exist_ok=True)
(out_dir / 'specz_counts.tex').write_text(table)
print(table)
