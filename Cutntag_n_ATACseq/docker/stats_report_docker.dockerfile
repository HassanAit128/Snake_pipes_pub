FROM python:3.9-slim

# Update and install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    libssl-dev \
    libffi-dev \
    python3-dev \
    python3-pip \
    samtools \ 
    && pip install openpyxl numpy pandas jinja2 matplotlib \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set environment variables
ENV LC_ALL=C
ENV LANG=C
