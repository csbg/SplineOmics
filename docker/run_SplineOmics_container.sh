# Run the analysis script
docker run -it --rm \
    -v $(pwd)/data.xlsx:/input/data.xlsx \
    -v $(pwd)/meta.xlsx:/input/meta.xlsx \
    -v $(pwd)/params.json:/input/params.json \
    -v $(pwd)/output:/output \
    splineomics:0.1.0 Rscript /app/run_analysis.R
