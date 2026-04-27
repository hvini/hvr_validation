# Dockerfile for CNN-based neuroimaging pipeline
FROM continuumio/miniconda3:latest

# Set working directory
WORKDIR /app

# Install dependencies via Conda
COPY env.yaml /tmp/env.yaml
RUN conda env create -f /tmp/env.yaml
RUN conda clean --all --yes

RUN /bin/bash -c "source /opt/conda/bin/activate hvr_cnn && \
    pip uninstall -y torch torchvision torchaudio && \
    pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu128"

# Install ImageMagick
RUN apt-get update && apt-get install -y imagemagick

# Copy files
COPY . /app

RUN sed -i 's/cpu=True/cpu=False/g' /app/cnn/cnn_hcvc.py

RUN sed -i 's/torch.load(params.model)/torch.load(params.model, weights_only=False)/g' /app/cnn/apply_segm_model.py

# Set environment variable for Python module path
ENV PYTHONPATH=/app/..

# Entrypoint
ENTRYPOINT ["conda", "run", "-n", "hvr_cnn", "python", "main.py"]
