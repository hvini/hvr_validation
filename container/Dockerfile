# Dockerfile for CNN-based neuroimaging pipeline
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Install dependencies via Conda
COPY env.yaml /tmp/env.yaml
RUN conda env create -f /tmp/env.yaml
RUN conda clean --all --yes

# Install ImageMagick
RUN apt-get update && apt-get install -y imagemagick

# Copy files
COPY . /app

# Set environment variable for Python module path
ENV PYTHONPATH=/app/..

# Entrypoint
ENTRYPOINT ["conda", "run", "-n", "hvr_cnn", "python", "main.py"]
