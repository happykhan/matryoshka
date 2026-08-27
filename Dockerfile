FROM python:3.12-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends ncbi-blast+ \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app
RUN python -m pip install --no-cache-dir .

ENTRYPOINT ["matryoshka"]
CMD ["--help"]
