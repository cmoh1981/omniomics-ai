# OmniOmics AI Dockerfile
# Multi-stage build for production

FROM python:3.11-slim as base

# Set environment variables
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    curl \
    git \
    wget \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install R and required packages for bioinformatic tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    r-base \
    r-base-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Install R packages
RUN R -e "install.packages(c('MOFA2', 'Seurat', 'tidyverse', ' patchwork', 'ggplot2'), repos='https://cloud.r-project.org/')"

# Copy application code
COPY . .

# Create non-root user for security
RUN useradd -m -u 1000 omniomics && \
    chown -R omniomics:omniomics /app
USER omniomics

# Expose ports
EXPOSE 8000 8501

# Default command
CMD ["python", "-m", "uvicorn", "omniomics_ai.api.main:app", "--host", "0.0.0.0", "--port", "8000"]
