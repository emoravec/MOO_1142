"""
Count the spectroscopic redshift sources per facility and write a LaTeX table.
"""
from pathlib import Path

from astropy.table import Table
# -------------------------------------------------------------------------------------------- #
specz_dir = Path(__file__).resolve().parent.parent
specz = Table.read(specz_dir / '2026-11-12.MOO_1142+1527.speczs.all.fits')
source = [str(s).strip() for s in specz['Source']]


def count(*names):
    return sum(s in names for s in source)


jwst_indices = [i for i, s in enumerate(source) if s.startswith('J')]
survey_indices = [i for i, s in enumerate(source) if s in ('DS', 'SD')]
internal_indices = [
    i for i, s in enumerate(source)
    if not s.startswith('J') and s not in ('DS', 'SD')
]


def is_q_gt_one(value):
    value = value.decode().strip() if isinstance(value, bytes) else str(value).strip()
    try:
        return float(value) > 1
    except ValueError:
        return False


def count_q_gt_one(indices):
    return sum(is_q_gt_one(specz['Q'][i]) for i in indices)


def count_in_cluster(indices):
    return sum(0.9 <= specz['Spec-z'][i] <= 1.3 for i in indices)


def count_good_in_cluster(indices):
    return sum(
        0.9 <= specz['Spec-z'][i] <= 1.3 and is_q_gt_one(specz['Q'][i])
        for i in indices
    )


n_jwst = len(jwst_indices)
n_internal = len(internal_indices)
n_survey = len(survey_indices)
n_jwst_clu = count_in_cluster(jwst_indices)
n_internal_q_gt_one = count_q_gt_one(internal_indices)
n_internal_clu = count_in_cluster(internal_indices)
n_internal_good_clu = count_good_in_cluster(internal_indices)
n_survey_clu = count_in_cluster(survey_indices)

rows = [
    ('JWST', n_jwst, f'JP ({count("JP")})', f'JG ({count("JG")})',
     n_jwst, n_jwst_clu, n_jwst_clu, ''),
    ('Internal', n_internal, f'KL/r ({count("KL", "KLr")})', f'KG ({count("KG")})',
     n_internal_q_gt_one, n_internal_clu, n_internal_good_clu, ''),
    ('Survey', n_survey, f'DS ({count("DS")})', f'SD ({count("SD")})',
     n_survey, n_survey_clu, n_survey_clu, ''),
    ('Total', len(source), '', '',
     n_internal_q_gt_one + n_jwst + n_survey,
     n_internal_clu + n_jwst_clu + n_survey_clu,
     n_internal_good_clu + n_jwst_clu + n_survey_clu, ''),
]

lines = [r'\begin{tabular}{lccccccc}',
         r'\hline',
         r' & $N_{\mathrm{tot}}$ & Facility 1 & Facility 2 & $N_{\mathrm{good}}$ & $N_{\mathrm{clu}}$ & $N_{\mathrm{clu},Q>1}$ & $N_{\mathrm{members}}$ \\',
         r'\hline']
lines += [' & '.join(str(c) for c in row) + r' \\' for row in rows]
lines += [r'\hline', r'\end{tabular}']
table = '\n'.join(lines) + '\n'

out_dir = specz_dir / 'latex_tables'
out_dir.mkdir(exist_ok=True)
(out_dir / 'specz_counts.tex').write_text(table)
print(table)
