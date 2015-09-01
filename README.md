# msgea
MSGEA (Multi-Species Gene Enrichment Analysis)

Package for enrichment analysis across 2 or 3 species.

Authors: Philipp Schiffer (EMBL), Luca Ferretti (The Pirbright Institute)


Standard single-species Fisher’s exact tests for enrichment cannot capture coherent signals of moderate enrichment across different species. MSGEA is sensitive to ontology terms enriched in several species at the same time. These represent a considerable fraction (30-70%) of the significant terms detected by our method in test data.
As far as we know, MSGEA is the first method for enrichment analysis that is able to detect these signals. Our method computes an exact p-value, grounded in the same statistics as the usual GO enrichment analyses, and in fact it reduces to the Fisher’s exact test once applied to a single species, thereby ensuring consistency with existing approaches.
