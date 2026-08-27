# Public acceptance sequences

These complete plasmid sequences are reproducible acceptance fixtures fetched from
NCBI Nucleotide by `scripts/fetch_acceptance_fixtures.py` on 2026-08-26.

| File | Accession | Purpose |
|---|---|---|
| `pKpQIL.fasta` | NC_014016.1 | Tn4401 reference detection |
| `pOXA-48a.fasta` | JN626286.1 | Negative Tn4401 control and Tn1999 boundary test |
| `pEK499.fasta` | EU935739.1 | ISEcp1 and false-positive acceptance tests |

`SHA256SUMS` pins the downloaded sequence after normalising the final newline. NCBI sequence records are public data;
the accession remains the authoritative source and attribution. Re-run the fetcher to
verify or deliberately refresh these fixtures.
