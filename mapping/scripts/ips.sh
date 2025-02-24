docker run --rm \
    -v $PWD/interproscan-5.72-103.0/data:/opt/interproscan/data \
    -v $PWD:/work \
    --memory=60g \
    interpro/interproscan:5.72-103.0 \
    --input data/top.faa \
    --output-dir data/jobs \
    --cpu 16 \
    --goterms \
    --iprlookup \
    --formats TSV