"""Read and filter HMMER domain tables.

`hmmer_import` parses an **hmmsearch** `--domtblout` file and applies the Methods' filters:
full-sequence E-value < 1e-4, independent domain E-value < 2, and a length check
`tlen > 0.9 * qlen`.

The search direction matters and is easy to get wrong. In hmmsearch output the target is the
*sequence* and the query is the *HMM*, so `tlen` is the ectodomain length and `qlen` the model
length -- which is what makes that third filter mean "the ECD is at least 90% as long as the
domain model", as the Methods describe. Running hmmscan instead swaps the two and silently
inverts the filter.
"""
import pandas as pd

def hmmer_import(filepath, filter = True, Evalue = 0.0001, i_Evalue = 2):
    hmmer_dom = pd.read_csv(filepath, skiprows = 3, header = None,
                    sep=r'\s+', skipfooter = 10, engine='python',
                   names = ['target_name','accession', 'tlen', 'query_name','Pfam', 'qlen', 'E-value_full',
                            'score_full', 'bias_full', 'dom_num', 'type_count', 'c-Evalue_dom','i-Evalue_dom',
                            'score_dom', 'bias_dom', 'hmm_from', 'hmm_to', 'align_from', 'align_to',
                            'env_from', 'env_to', 'accuracy','description', 'nothing'])
    hmmer_dom.drop(columns = ['accession', 'description', 'nothing'], inplace = True)

    if filter:
        #filtering for specific domain assignments and then lest specific individual domains seems to give best agreement with notch1
        filtered_dom = hmmer_dom[(hmmer_dom['E-value_full']  < Evalue) & (hmmer_dom['i-Evalue_dom']  < i_Evalue)]
        filtered_dom = filtered_dom.assign(Pfam = filtered_dom['Pfam'].str.split('.').str[0])

        #length of domain assignement must cover 90% of the domain hmm profile
        filtered = filtered_dom[filtered_dom['tlen'] > 0.9*filtered_dom['qlen']]

        return filtered

    else:
        return hmmer_dom
