import streamlit as st
import pandas as pd
import re
from collections import Counter
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from goatools.base import get_godag
import plotly.express as px
import io



class GOTerm:
    def __init__(self, id, name, namespace, level, depth):
        self.id = id
        self.name = name
        self.namespace = namespace
        self.level = level
        self.depth = depth

    def get_all_upper(self):
        # Placeholder: Return a list of ancestor GO terms
        return []


godag = get_godag('go-basic.obo', optional_attrs='relationship')


def extract_GO_terms(row):
    row_str = ' '.join(row.astype(str))
    go_terms = re.findall(r'GO:\d+', row_str)
    unique_go_terms = list(set(go_terms))
    return unique_go_terms, row['protein_id']


def extract_ancestors(go_id):
    if go_id not in godag:
        return []
    term = godag[go_id]
    ancestors = term.get_all_upper()
    return [go_id] + list(ancestors)


def GO_enrichments(df, pvalue, refined=True, l=0, d=0, namespaces=None,filtered=["regulation"]):
    df['GOs_and_protein'] = df.apply(extract_GO_terms, axis=1)

    go_to_proteins = {}

    newGOs_background = []
    for i, row in df.iterrows():
        for go in row['GOs_and_protein'][0]:
            ancestors = extract_ancestors(go)
            for ancestor in ancestors:
                GO = godag[ancestor].namespace + " - " + godag[ancestor].name + " L" + str(
                    godag[ancestor].level) + " D" + str(godag[ancestor].depth)
                newGOs_background.append(GO)
                if GO not in go_to_proteins:
                    go_to_proteins[GO] = set(row['GOs_and_protein'][1])
                go_to_proteins[GO].add(row['GOs_and_protein'][1])

    background = newGOs_background
    if namespaces:
        background = [x for x in background if any(ns in x for ns in namespaces)]
    if l != 0:
        background = [x for x in background if int(x.split(" ")[-2][1:]) >= l]
    if d != 0:
        background = [x for x in background if int(x.split(" ")[-1][1:]) >= d]

    col = "pvalue_refined" if refined else "pvalue"
    df_cut = df[df[col] < pvalue]
    GOs = df_cut["GOs_and_protein"].apply(lambda x: x[0]).tolist()
    GOs = [item for sublist in GOs for item in sublist]
    GOs = [item for sublist in (extract_ancestors(x) for x in GOs) for item in sublist]
    # create new go_to_proteins that has only pids that are in df_cut['protein_id'].tolist()
    new_go_to_proteins = {}
    for go in go_to_proteins:
        new_go_to_proteins[go] = go_to_proteins[go].intersection(
            df_cut['GOs_and_protein'].apply(lambda x: x[1]).tolist())
        new_go_to_proteins[go] = [pid + " (" + df.loc[df['protein_id'] == pid, 'product'].values[0] + ")" for pid in
                                  new_go_to_proteins[go]]
    newGOs = []
    for x in GOs:
        GO = godag[x].namespace + " - " + godag[x].name + " L" + str(godag[x].level) + " D" + str(godag[x].depth)
        newGOs.append(GO)

    if namespaces:
        newGOs = [x for x in newGOs if any(ns in x for ns in namespaces)]
    if l != 0:
        newGOs = [x for x in newGOs if int(x.split(" ")[-2][1:]) >= l]
    if d != 0:
        newGOs = [x for x in newGOs if int(x.split(" ")[-1][1:]) >= d]
    newGOs = [x for x in newGOs if not any(f in x for f in filtered)]

    GOs = Counter(newGOs)
    GOs_background = Counter(background)
    M = len(df)
    n = len(df_cut)

    if n == 0 or len(GOs) == 0:
        print("No GO terms meet the specified criteria. Try adjusting your filters.")
        return pd.DataFrame()

    p_values = {}
    enrichments = {}

    for go in GOs:
        N = GOs_background[go]
        x = GOs[go]
        if N > 0 and M > 0:
            enrichments[go] = (x / n) / (N / M)
            observed_frequency = x / n
            expected_frequency = N / M
            enrichment_score = observed_frequency / expected_frequency
            p_values[go] = hypergeom.sf(x - 1, M, N, n)
        else:
            enrichments[go] = 0
            p_values[go] = 1

    if len(p_values) == 0:
        print("No valid GO terms for enrichment analysis. Try adjusting your filters.")
        return pd.DataFrame()

    pvals = list(p_values.values())
    corrected_pvals_bonferroni = multipletests(pvals, method='bonferroni')[1]
    corrected_pvals_benjamini_hochberg = multipletests(pvals, method='fdr_bh')[1]

    data = {
        'GO term': list(p_values.keys()),
        'p-value': list(p_values.values()),
        'corrected p-value Bonferroni': corrected_pvals_bonferroni,
        'corrected p-value Benjamini-Hochberg': corrected_pvals_benjamini_hochberg,
        'associated_proteins': [', '.join(new_go_to_proteins.get(go, [])) for go in p_values.keys()]
    }
    dfGO = pd.DataFrame(data)
    dfGO['enrichment'] = dfGO['GO term'].map(enrichments)
    dfGO = dfGO.sort_values(by='p-value', ascending=True)
    return dfGO


def main():
    st.set_page_config(layout="wide")
    st.title("Applicazione per Gene Enrichment Analysis di dati tabulari")
    st.markdown("""
    Questa app esegue l'analisi di arricchimento genico basata sui termini Gene Ontology (GO) usando come set i termini GO associati a un p-value sotto una certa soglia.
    Carica il tuo file di dati genici (.csv o .xlsx) e regola i parametri per eseguire l'analisi. Ci può essere un numero variabile di colonne con termini GO, lo script estrae tutti i termini GO per ogni riga. Ci deve essere una colonna chiamata pvalue e opzionale una chiamata pvalue_refined
    """)
    use_default = st.checkbox("Usa la tabella dati di default", value=True)


    if use_default:
        #load file github https://github.com/DavideGotta/TESI/blob/master/operoniscorestream.csv
        uploaded_file = 'https://raw.githubusercontent.com/DavideGotta/TESI/master/operoniscorestream.csv'
        st.info("Deseleziona la casella qui sopra per caricare il tuo file.")
    else:
        uploaded_file = st.file_uploader("Choose a CSV or XLSX file", type=["csv", "xlsx"])

    if uploaded_file is not None:
        try:
            if isinstance(uploaded_file, str) or uploaded_file.name.endswith('.csv'):
                df = pd.read_csv(uploaded_file)
            else:
                df = pd.read_excel(uploaded_file)

            st.subheader("Data Preview")
            st.dataframe(df, use_container_width=True)

            col1, col2 = st.columns(2)

            with col1:
                pvalue_threshold = st.slider("Seleziona la soglia per il p-value", 0.0, 1.0, 0.01, 0.001)
                use_refined = st.checkbox("Usa il p-value dello score raffinato (che usa l'informazione degli ortologhi)", value=True)

                st.markdown("""
                **Proprietà dei termini GO:**
                I termini GO troppo in alto nella gerarchia sono troppo generici e potrebbero risultare poco utili al fine dell'analisi, per questo si può filtrare in base alla posizione nella gerarchia di Gene Ontology.
                - **Level**: La lunghezza del percorso più breve dal termine radice. I valori più bassi includono termini più generali.
                - **Depth**: La lunghezza del percorso più lungo dal termine radice. I valori più alti includono termini più specifici.
                """)

                min_level = st.number_input("Min level per i termini GO", min_value=0, value=4, step=1)
                min_depth = st.number_input("Min depth per i termini GO", min_value=0, value=0, step=1)

            with col2:
                namespace_options = ["biological_process", "cellular_component", "molecular_function"]
                selected_namespaces = st.multiselect(
                    "Select GO namespaces to include",
                    options=namespace_options,
                    default=["biological_process"]
                )
                st.markdown("**Filtra termini GO per parole chiave:**")
                default_keywords = ["regulation"]
                additional_keywords = st.text_input("Aggiungi parole chiave separate da virgola", "")
                all_keywords = default_keywords + [kw.strip() for kw in additional_keywords.split(",") if kw.strip()]
                selected_keywords = st.multiselect("Seleziona le parole chiave da filtrare", options=all_keywords, default=default_keywords)


            if st.button("Run Gene Enrichment Analysis"):
                with st.spinner('Running analysis...'):
                    results = GO_enrichments(df, pvalue_threshold, refined=use_refined,
                                             l=min_level, d=min_depth, namespaces=selected_namespaces, filtered=selected_keywords)

                if not results.empty:
                    st.subheader("Gene Enrichment Results")

                    # Display results in an interactive table
                    st.dataframe(results.style.format({
                        'p-value': "{:.2e}",
                        'corrected p-value Bonferroni': "{:.2e}",
                        'corrected p-value Benjamini-Hochberg': "{:.2e}",
                        'enrichment': "{:.2f}"
                    }), use_container_width=True)

                    # Create a bar plot of top enriched GO terms
                    top_terms = results.sort_values('enrichment', ascending=False).head(10)
                    fig = px.bar(top_terms, x='enrichment', y='GO term', orientation='h',
                                 title='Top 10 Enriched GO Terms',
                                 labels={'enrichment': 'Enrichment Score', 'GO term': 'GO Term'},
                                 color='p-value', color_continuous_scale='Viridis')
                    st.plotly_chart(fig, use_container_width=True)

                    # Download button
                    csv = results.to_csv(index=False)
                    st.download_button(
                        label="Download results as CSV",
                        data=csv,
                        file_name="gene_enrichment_results.csv",
                        mime="text/csv",
                    )
                else:
                    st.warning("No results to display. Try adjusting your filters.")
        except Exception as e:
            st.error(f"An error occurred: {str(e)}")


if __name__ == "__main__":
    main()
