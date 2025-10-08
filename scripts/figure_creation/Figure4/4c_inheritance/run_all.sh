for f in /Users/noah/Desktop/ufc_repository/results/analysis_4c_results/*.tsv; do
    echo "Running on $f"
    python3 create_figures.4c.py "$f"
done
