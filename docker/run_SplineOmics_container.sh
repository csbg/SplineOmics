docker run -it -d \
    -v $(pwd)/data.xlsx:/input/data.xlsx \
    -v $(pwd)/meta.xlsx:/input/meta.xlsx \
    -v $(pwd)/output:/output \
    -p 8888:8787 \
    -e PASSWORD=password \
    thomasrauter/splineomics:0.1.0



