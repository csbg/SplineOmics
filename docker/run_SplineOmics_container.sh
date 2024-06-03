# Run the analysis script
docker run -it --rm \
    -v $(pwd)/data.xlsx:/workspace/data.xlsx \
    -v $(pwd)/meta.xlsx:/workspace/meta.xlsx \
    -v $(pwd)/params.json:/workspace/params.json \
    -v $(pwd)/output:/output \
    splineomics-env:0.1.0 Rscript /workspace/run_analysis.R
