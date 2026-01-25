for f in /Users/noah/Desktop/ufc_repository/results/analysis_4e_results/4e_prs_results/*.tsv; do
    echo "Running on $f"
    python3 create_figures.4e.v3.py "$f"
done
