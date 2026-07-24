FROM pytorch/pytorch:2.2.1-cuda12.1-cudnn8-runtime

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    libz-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    numpy==1.23.5 \
    matplotlib==3.5.1 \
    pysam==0.19.1

WORKDIR /app

COPY src/swave /app/swave

ENV PYTHONPATH="/app/swave"

CMD ["python3"]